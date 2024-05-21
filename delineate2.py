"""
Fast delineation of waterhshed subbasins, using data from
MERIT-Basins and MERIT-Hydro

Creates watersheds for outlet points in an input CSV file.

Usage: Carefully edit the file config2.py, then run delineator2.py
See README for more detailed instructions.

"""

import warnings
import numpy as np
import pickle
import pandas as pd
import os
import geopandas as gpd
import re
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
import shapely.ops
from shapely.wkt import loads
import sigfig  # for formatting numbers to significant digits
from py.fast_dissolve import dissolve_geopandas, fill_geopandas
import pyproj
from functools import partial
from config2 import *
from py.mapper import make_map, create_folder_if_not_exists
from py.merit_detailed import split_catchment
import matplotlib.pyplot as plt

from plot_network import draw_graph
from util.graph_tools import calculate_strahler_stream_order, calculate_shreve_stream_order, calculate_num_incoming, \
    make_river_network, insert_node

warnings.simplefilter(action='ignore', category=FutureWarning)

# The WGS84 projection string, used in a few places
PROJ_WGS84 = 'EPSG:4326'

# Regular expression used to find numbers so I can round lat, lng coordinates in GeoJSON files to make them smaller
simpledec = re.compile(r"\d*\.\d+")


def has_unique_elements(lst: list) -> bool:
    """
    Determine whether all the items in a list are unique
    :param lst: A Python LIST
    :return: boolean,
        True if all items in the list are unique
        False if there are any duplicated items
    """
    return len(lst) == len(set(lst))


def find_repeated_elements(lst: list) -> list:
    """
    Finds any repeated (or duplicate) items in a Python list.
    Input argument: a List
    Outputs: a list containing any duplicate values.
      The duplicates are not repeated, we just id the values that were repeated.
      If there were no dupes, this function returns an empty list []
    """
    seen = set()
    duplicates = set()
    for elem in lst:
        if elem in seen:
            duplicates.add(elem)
        else:
            seen.add(elem)
    return list(duplicates)


def mround(match):
    # Utility function for rounding the coordinates in GeoJSON files to make them smaller
    return "{:.5f}".format(float(match.group()))


def validate(gages_df: pd.DataFrame) -> bool:
    """
    After we have read in the user's input CSV file with their desired delineation locations
    (I refer to these as gages as that was my original use case), check whether the input
    data is valid.
    (1) required columns are present -- at a minimum id, lat, lng)
    (2) data types are correct
    (3) input values are in the appropriate range (i.e. lat between -90 and 90).

    returns: True, if inputs are valid
             throws an Exception if inputs are not valid.

    """
    cols = gages_df.columns
    required_cols = ['id', 'lat', 'lng']
    for col in required_cols:
        if col not in cols:
            raise Exception(f"Missing column in CSV file: {col}")

    # Check that the ids are all unique
    if len(gages_df['id'].unique()) != len(gages_df):
        raise Exception("Each id in your CSV file must be unique.")

    # Check that lat, lng are numeric
    fields = ['lat', 'lng']
    for field in fields:
        if gages_df[field].dtype != 'float64':
            raise Exception(f"In outlets CSV, the column {field} is not numeric.")

    # Check that all the lats are in the right range
    lats = gages_df["lat"].tolist()
    lngs = gages_df["lng"].tolist()

    if not all(lat > -60 for lat in lats):
        raise Exception("All latitudes must be greater than -60°")

    if not all(lat < 85 for lat in lats):
        raise Exception("All latitudes must be less than 85°")

    if not all(lng > -180 for lng in lngs):
        raise Exception("All longitudes must be greater than -180°")

    if not all(lng < 180 for lng in lngs):
        raise Exception("All longitudes must be less than 180°")

    # Check that every row has an id
    ids = gages_df["id"].tolist()

    if not all(len(str(wid)) > 0 for wid in ids):
        raise Exception("Every watershed outlet must have an id in the CSV file")

    # Check that the ids are unique. We cannot have duplicate ids, because they are used as the index in DataFrames
    if not has_unique_elements(ids):
        raise Exception("Outlet ids must be unique. No duplicates are allowed!")

    return True


def calc_area(poly: Polygon) -> float:
    """
    Calculates the approximate area of a Shapely polygon in raw lat, lng coordinates (CRS=4326)
    First projects it into the Albers Equal Area projection to facilitate calculation.
    No
    Args:
        poly: Shapely polygon
    Returns:
         area of the polygon in km²
    """
    if poly.is_empty:
        return 0

    projected_poly = shapely.ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat_1=poly.bounds[1],
                lat_2=poly.bounds[3]
            )
        ),
        poly)

    # Get the area in m^2
    return projected_poly.area / 1e6


def calc_length(line: LineString) -> float:
    """
    Calculates the approximate length in km of a Shapely LineString in raw lat, lng coordinates (CRS=4326)
    First projects it into the Albers Equal Area projection to facilitate calculation.

    Args:
        line: Shapely LineString
    Returns:
         length of the LineString in kilometers.
    """
    if line.is_empty:
        return 0

    projected_line = shapely.ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat_1=line.bounds[1],
                lat_2=line.bounds[3]
            )
        ),
        line)

    # Get the area in m^2
    return projected_line.length / 1e3


def load_megabasins() -> gpd.GeoDataFrame:
    """
    Reads the "megabasin" data from disk and returns a GeoDataFrame.
    The program uses this data to determine what dataset is needed for analyses.
    I refer to the MERIT-Basins Pfafstetter Level 2 basins as megabasins.

    This function gets the data from a shapefile or a pickle file, if it exists.
    If the .pkl does not exist, create it for faster processing in the future,
    since reading shapefiles is slow.
    """

    # Check whether the pickle file exists
    pickle_fname = f"{PICKLE_DIR}/megabasins.pkl"
    if os.path.isfile(pickle_fname):
        gdf = pickle.load(open(pickle_fname, "rb"))
        return gdf

    else:
        # This file has the merged "megabasins_gdf" in it
        merit_basins_shp = 'data/shp/basins_level2/merit_hydro_vect_level2.shp'
        megabasins_gdf = gpd.read_file(merit_basins_shp)

        # The CRS string in the shapefile is EPSG 4326 but does not match verbatim, so set it here
        megabasins_gdf.to_crs(PROJ_WGS84, inplace=True)
        if not megabasins_gdf.loc[0].BASIN == 11:
            raise Exception("An error occurred loading the Level 2 basins shapefile")

        # Try saving to a pickle file for future speedups
        if len(PICKLE_DIR) > 0:
            make_folders()
            pickle.dump(megabasins_gdf, open(pickle_fname, "wb"))

    return megabasins_gdf


def get_megabasin(points_gdf) -> int:
    """
    Finds out what Pfafstetter Level 2 "megabasin" a point is in.
    Arguments:
        lat, lon of a point on the map
    Returns:
        the ID of the megabasin, an integer from 11 to 91
    """
    if VERBOSE: print("Finding out which Level 2 megabasin(s) your outlets are in")

    megabasins_gdf = load_megabasins()

    # Overlay the gage points on the Level 2 Basins polygons to find out which
    # PFAF_2 basin each point falls inside of, using a spatial join
    if SEARCH_DIST == 0:
        gages_basins_join = gpd.overlay(points_gdf, megabasins_gdf, how="intersection")

    # Needed to set this option in order to avoid a warning message in Geopandas.
    # https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
    pd.options.mode.chained_assignment = None  # default='warn'

    # Get a list of the DISTINCT Level 2 basins, and a count of how many gages in each.
    basins_df = gages_basins_join.groupby("BASIN").id.nunique()
    basins = basins_df.index.tolist()

    # If the points fall in more than one megabasin, not possible to process: Raise Error
    if len(basins) > 1:
        print(f"ERROR: Your watershed outlets are in {len(basins)} continental-scale megabasin(s). "
              f"Fix before continuing.")
        raise Exception

    # Function returns the megabasin of the first point (all are in the same megabasin)
    return basins[0]


def make_folders():
    """
    This function makes sure that the folders that the user specified
    in config.py are existant. If they are not, it tries to create them.
    If it cannot find the folder, and it cannot create the folder, the
    function will raise an error.

    :return: Nothing, but throws an error if fails.
    """
    # Check that the OUTPUT directories are there. If not, try to create them.
    folder_exists = create_folder_if_not_exists(OUTPUT_DIR)
    if not folder_exists:
        raise Exception(f"No folder for output. Stopping")
    # Check for the folder to put Python PICKLE files
    if PICKLE_DIR != "":
        folder_exists = create_folder_if_not_exists(OUTPUT_DIR)
        if not folder_exists:
            raise Exception(f"No folder for pickle files. Stopping")
    # Check if the MAP_FOLDER is there
    if MAKE_MAP:
        folder_exists = create_folder_if_not_exists(MAP_FOLDER)
        if not folder_exists:
            raise Exception(f"No folder for the map files. Stopping")


def get_pickle_filename(geotype: str, basin: int, high_resolution: bool) -> str:
    """Simple function to get the standard filename for the pickle files used by this project.
    The filenames look like this:
       PICKLE_DIR/catchments_##_hires.pkl
       PICKLE_DIR/catchments_##_lores.pkl

       PICKLE_DIR/rivers_##_hires.pkl
       PICKLE_DIR/rivers_##_lores.pkl

    where ## is the megabasin number (11-91)

    """

    if high_resolution:
        resolution_str = 'hires'
    else:
        resolution_str = 'lores'
    fname = f'{PICKLE_DIR}/{geotype}_{basin}_{resolution_str}.pkl'
    return fname


def load_gdf(geotype: str, basin: int, high_resolution: bool) -> gpd.GeoDataFrame:
    """
    Returns the unit catchments vector polygon dataset as a GeoDataFrame
    Gets the data from the MERIT-Basins shapefile the first time,
    and after that from a saved .pkl file on disk.
    Uses some global parameters from config.py

    :param geotype: either "catchments" or "rivers" depending on which one we want to open.
    :param basin: the Pfafstetter level 2 megabasin, an integer from 11 to 91
    :param high_resolution: True to load the standard (high-resolution) file,
      False to load the low-resolution version (for faster processing, slightly less accurate results)

    :return: a GeoPandas GeoDataFrame

    """

    # First, check for the presence of a pickle file
    if PICKLE_DIR != '':
        pickle_fname = get_pickle_filename(geotype, basin, high_resolution)
        if os.path.isfile(pickle_fname):
            if VERBOSE: print(f"Fetching BASIN # {basin} catchment data from pickle file.")
            gdf = pickle.load(open(pickle_fname, "rb"))
            return gdf

    # Open the shapefile for the basin
    if geotype == "catchments":
        if high_resolution:
            directory = HIGHRES_CATCHMENTS_DIR
        else:
            directory = LOWRES_CATCHMENTS_DIR
        shapefile = f"{directory}/cat_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"
    elif geotype == "rivers":
        shapefile = f"{RIVERS_DIR}/riv_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"

    if not os.path.isfile(shapefile):
        raise Exception(f"Could not find the file: {shapefile}")

    if VERBOSE: print(f"Reading geodata in {shapefile}")
    gdf = gpd.read_file(shapefile)
    #gdf.set_index('COMID', inplace=True)

    # This line is necessary because some of the shapefiles provided by reachhydro.com do not include .prj files
    gdf.set_crs(PROJ_WGS84, inplace=True, allow_override=True)

    # Before we exit, save the GeoDataFrame as a pickle file, for future speedups!
    save_pickle(geotype, gdf, basin, high_resolution)
    return gdf


def save_pickle(geotype: str, gdf: gpd.GeoDataFrame, basin: int, high_resolution: bool):
    # If we loaded the catchments from a shapefile, save the gdf to a pickle file for future speedup
    if PICKLE_DIR != '':

        # Check whether the GDF has a spatial index.
        # Note: I don't think this is ever necessary. Since version 0.7.0 (March 2020), GeoPandas
        # creates a spatial index by default.
        has_spatial_index = hasattr(gdf, 'sindex')
        if not has_spatial_index:
            gdf.sindex.create_index()

        # Get the standard project filename for the pickle files.
        pickle_fname = get_pickle_filename(geotype, basin, high_resolution)
        if not os.path.isfile(pickle_fname):
            if VERBOSE: print(f"Saving GeoDataFrame to pickle file: {pickle_fname}")
            try:
                pickle.dump(gdf, open(pickle_fname, "wb"))
            except:
                raise Warning("Could not save pickle file to: {pickle_fname}")


def fix_polygon(poly: Polygon or MultiPolygon) -> Polygon:
    """
    When we use the difference() method in Shapely to subtract one polygon from another,
    it's common to end up with small slivers or unwanted geometries around the edges of
    the input polygon due to precision issues or complex boundaries.
    To eliminate these small slivers, we can use a combination of simplification and filtering based on area.

    input: a Shapely Polygon or MultiPolygon
    output: a Shapely Polygon that has been fixed to remove small slivers.
    """
    simplify_tolerance = 0.0001
    simplified_poly = poly.simplify(tolerance=simplify_tolerance, preserve_topology=True)

    # Filter out small slivers by area
    min_area_threshold = 0.0001  # Define your own threshold

    # Handle MultiPolygon and Polygon cases
    if simplified_poly.geom_type == 'Polygon':
        cleaned_poly = simplified_poly if simplified_poly.area >= min_area_threshold else Polygon()
    elif simplified_poly.geom_type == 'MultiPolygon':
        cleaned_poly = MultiPolygon([poly for poly in simplified_poly.geoms if poly.area >= min_area_threshold])

    else:
        cleaned_poly = Polygon()  # If it's neither, return an empty Polygon

    # Optional: Simplify again if needed
    final_poly = cleaned_poly.simplify(tolerance=simplify_tolerance, preserve_topology=True)

    return final_poly


def write_geodata(gdf: gpd.GeoDataFrame, fname: str):
    """
    Write a GeoDataFrame to disk in the user's pre
    """
    if VERBOSE: print('Writing geodata to disk')

    # This line rounds all the vertices to fewer digits. For text-like formats GeoJSON or KML, makes smaller
    # files with minimal loss of precision. For other formats (shp, gpkg), it doesn't change file size, so don't bother.
    if OUTPUT_EXT.lower() in ['geojson', 'kml']:
        gdf.geometry = gdf.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning)
        gdf.to_file(fname)


def plot_basins(basins_gdf: gpd.GeoDataFrame, outlets_gdf: gpd.GeoDataFrame, fname: str):
    """
    Makes a plot of the unit catchments that are in the watershed

    """
    # subbasins_gdf.plot(column='area', edgecolor='gray', legend=True)
    [fig, ax] = plt.subplots(1, 1, figsize=(10, 8))

    # Plot each unit catchment with a different color
    for x in basins_gdf.index:
        color = np.random.rand(3, )
        basins_gdf.loc[[x]].plot(facecolor=color, edgecolor=color, alpha=0.5, ax=ax)

    # Plot the gage points
    outlets_gdf.plot(ax=ax, c='red', edgecolors='black')

    plt.savefig(f"plots/{fname}.png")
    plt.close(fig)


def delineate2():
    """
    MAIN watershed delineation routine
    Make sure to set the variables in `config2.py` before running.

    Version #2, where all the points are in a single watershed,
    and we are interested in returning a set of subwatersheds.
    Reads a list of outlet points from a .csv file,
    the first row should contain the main watershed outlet.
    Subsequent rows of the .csv should contain uptstream points
    within the watershed that we will use to create subwatershed outlets.

    Outputs geodata (.shp, .gpkg, etc.) and optionally a CSV file with a summary of results
    (including the watershed id, names, and areas).

    """

    def addnode(B, node):
        """"
        Recursive function to assemble the list of upstream unit catchments
        B is a Python List of the unit catchments that make up our watershed.
        B is for BASIN...
        List items are `comid`s, unique identifiers of each unit catchment.

        """
        # first, append the node to the list that is in the Basin, B
        B.append(node)

        # next, check whether the fields up1, up2, up3, and up4 contain a node ID
        up1 = rivers_gdf['up1'].loc[node]
        if up1 != 0:
            addnode(B, up1)

        up2 = rivers_gdf['up2'].loc[node]
        if up2 != 0:
            addnode(B, up2)

        up3 = rivers_gdf['up3'].loc[node]
        if up3 != 0:
            addnode(B, up3)

        up4 = rivers_gdf['up4'].loc[node]
        if up4 != 0:
            addnode(B, up4)

    # Make sure the user folders from `config.py` exist. If not create them.
    make_folders()

    # Check that the CSV file is there
    if not os.path.isfile(OUTLETS_CSV):
        raise Exception(f"Could not your outlets file at: {OUTLETS_CSV}")

    # Read the outlet points CSV file and put the data into a Pandas DataFrame
    # (I call the outlet points gages, because I usually in delineated watersheds at streamflow gages)
    if VERBOSE: print(f"Reading your outlets data in: {OUTLETS_CSV}")
    gages_df = pd.read_csv(OUTLETS_CSV, header=0, dtype={'id': 'str', 'lat': 'float', 'lng': 'float'})

    # Check that the CSV file includes at a minimum: id, lat, lng and that all values are appropriate
    validate(gages_df)

    # Get the number of points in the gages file
    n_gages = len(gages_df)

    # Boolean vars to track whether user's outlets file had fields `name`
    bNames = 'name' in gages_df

    # Create fields for the "snapped" points, and the distance from the input coordinates
    gages_df['lat_snap'] = np.nan
    gages_df['lng_snap'] = np.nan
    gages_df['snap_dist'] = 0

    # Convert gages_df to a GeoPandas GeoDataFrame (adds geography, lets us do geo. operations)
    coordinates = [Point(xy) for xy in zip(gages_df['lng'], gages_df['lat'])]
    points_gdf = gpd.GeoDataFrame(gages_df, crs=PROJ_WGS84, geometry=coordinates)

    # Get the megabasin(s) of the points
    megabasin = get_megabasin(points_gdf)

    # Now add fields to gages_df so we can reuse it to create a table to output to CSV
    #gages_df.set_index('id', inplace=True)
    gages_df['area_calc'] = 0
    gages_df['result'] = "failed"
    num_gages = len(gages_df)

    # Dict failed: key = id, value = string, explanation of failure
    failed = {}

    gages_counter = 0

    if VERBOSE: print('Reading data table for unit catchments in basin %s' % megabasin)
    catchments_gdf = load_gdf("catchments", megabasin, True)

    # The network data is in the RIVERS file rather than the CATCHMENTS file
    # (this is just how the MeritBASIS authors did it)
    if VERBOSE: print('Reading data table for rivers in basin %s' % megabasin)
    rivers_gdf = load_gdf("rivers", megabasin, True)
    rivers_gdf.set_index('COMID', inplace=True)

    # Perform a Spatial join on gages (points) and the unit catchments (polygons)
    # to find the corresponding unit catchment for each gage
    # Adds the fields COMID and unitarea
    if VERBOSE: print(f"Performing spatial join on {num_gages} outlet points in basin #{megabasin}")
    gages_joined = gpd.overlay(points_gdf, catchments_gdf, how="intersection")

    # For any gages for which we could not find a unit catchment, add them to Dict failed
    gages_matched = gages_joined['id'].tolist()
    gage_basin_ids = gages_df['id'].tolist()
    for wid in gage_basin_ids:
        if wid not in gages_matched:
            failed[wid] = f"Could not assign to a unit catchment in Level 2 basin #{megabasin}"

    # First, let us find the set of unit catchments upstream of the outlet.
    # Let wid be the watershed ID. Get the lat, lng coords of the gage.
    terminal_node_id = gages_joined['id'].iloc[0]

    # The terminal comid is the unit catchment that contains (overlaps) the outlet point
    terminal_comid = gages_joined['COMID'].iloc[0]

    # Let B be the list of unit catchments (and river reaches) that are in the basin
    B = []

    # Add the first node, and the rest will be added recursively
    addnode(B, terminal_comid)

    # Next, check that all the other points provided by the user are
    # in this set; otherwise, it means they are not upstream.
    n = len(gages_joined)
    problemFound = False
    for i in range(0, n):
        id = gages_joined.at[i, 'id']
        comid = gages_joined.at[i, 'COMID']
        if comid not in B:
            print(f"Problem: Point with id = {id} is in the watershed of the first point.")
            problemFound = True

    if problemFound:
        raise Exception("One or more points not in watershed. Please check your inputs and try again.")

    # Next, check that we are not trying to subdivide a single unit catchment more than once.
    comids = gages_joined['COMID'].tolist()

    if not has_unique_elements(comids):
        print("Problem: More than one point falls in each unit catchment. "
              "This script is not currently set up to handle this")

        # Get a list of the offending points
        repeated_comids = find_repeated_elements(comids)
        for comid in repeated_comids:
            ids_rows = gages_joined[gages_joined['comid'] == comid, 'id']
            ids_list = ids_rows.tolist()
            print(f"The following gages are in unit catchment {comid}: {', '.join(ids_list)}")

    # Next, we need to split every subcatchment that contains an outlet point.
    # First, we will create a DataFrame with the subcatchment data.
    catchments_gdf.set_index('COMID', inplace=True)
    subbasins_gdf = catchments_gdf.loc[B]

    plot_basins(subbasins_gdf, points_gdf, 'before')

    # We are going to add a new subbasin for each outlet point!
    # This info is in the rivers dataframe, geom field. It is the *START* point (!)
    rivers_gdf['end_point'] = rivers_gdf['geometry'].apply(lambda x: x.coords[0])
    rivers_gdf['lng'] = rivers_gdf['end_point'].apply(lambda x: x[0])
    rivers_gdf['lat'] = rivers_gdf['end_point'].apply(lambda x: x[1])

    # Drop the 'end_point' column as it's no longer needed
    rivers_gdf = rivers_gdf.drop(columns=['end_point'])

    # Join the GeoDataFrames on their indices
    subbasins_gdf = subbasins_gdf.join(rivers_gdf[['lat', 'lng', 'NextDownID']])

    # Re-name a field for my convenience, and make sure it is an integer
    subbasins_gdf.rename(columns={'NextDownID': 'nextdown'}, inplace=True)
    subbasins_gdf['nextdown'] = subbasins_gdf['nextdown'].astype(int)

    # The next steps are easier to figure out using graphs.
    # We are going to insert nodes.
    G = make_river_network(subbasins_gdf, terminal_comid)

    G = calculate_strahler_stream_order(G)
    if PLOTS: draw_graph(G, 'plots/before')

    # Create a dictionary of the new nodes to add : the comid that we're inserting them into!
    gages_joined.set_index('id', inplace=True)
    new_nodes = gages_joined['COMID'].to_dict()

    # Insert the new nodes into the flow network.
    for node, comid in new_nodes.items():
        G = insert_node(G, node, comid)

    if PLOTS: draw_graph(G, 'plots/after')

    # Now we have a new, accurate network topology. Now we only need to take care of the geography,
    # by splitting the polygons to find the geography (catchment polygon) for the
    # the newly inserted nodes in the river network.

    # First, we can insert the new nodes into our subbasins GeoDataFrame,
    # and update the topology data (column `nextdown`) based on the Graph.
    G = calculate_strahler_stream_order(G)
    G = calculate_shreve_stream_order(G)

    # First insert the new nodes
    for node, comid in new_nodes.items():
        lat = gages_joined.at[node, 'lat']
        lng = gages_joined.at[node, 'lng']
        subbasins_gdf.loc[node] = [None, None, lat, lng, int(comid)]

    # Kind of hackish, but coerce the field nextdown back to integer type
    # to fix the problem where the comids were being turned into floats
    subbasins_gdf['nextdown'] = subbasins_gdf['nextdown'].astype(int)

    # Append the column `type` (to hold 'leaf' or 'stem')
    # This will be important later when we squash and simplify the river network,
    # because we MUST NOT remove these user-defined nodes!
    # TODO: Consider making this just 'new' ?
    subbasins_gdf['type'] = ''

    # Let's update the field `nextdown` by iterating over the rows and getting the data from the GRAPH.
    for index in subbasins_gdf.index:
        try:
            nextdown = list(G.out_edges(index))[0][1]
            subbasins_gdf.at[index, 'nextdown'] = nextdown
        except IndexError:
            subbasins_gdf.at[index, 'nextdown'] = None

        # Add whether it is a leaf or a stem
        try:
            node_type = G.nodes[index]['type']
        except KeyError:
            node_type = ''
        subbasins_gdf.at[index, 'type'] = node_type

    # Next, run the split_catchment() routine to find the geometry of the new nodes!
    for node, comid in new_nodes.items():
        lat = gages_joined.at[node, 'lat']
        lng = gages_joined.at[node, 'lng']
        catchment_poly = subbasins_gdf.loc[comid, 'geometry']
        node_type = subbasins_gdf.loc[node, 'type']
        isLeaf = (node_type == 'leaf')

        # Split catchment routine
        new_poly, lat_snap, lng_snap = split_catchment(node, megabasin, lat, lng, catchment_poly, isLeaf)

        # Assign the split polygon as the geometry of the new node.
        #  i.e. this is the subcatchment area corresponding to a user-defined outlet.
        subbasins_gdf.at[node, 'geometry'] = new_poly

        # Calculate the area of the new unit catchment and insert value into GeoDataFrame
        area = calc_area(new_poly)
        subbasins_gdf.at[node, 'unitarea'] = round(area, 1)

        # Now SUBTRACT the new polygon from the target unit catchments polygon
        # and update the geometry of the unit catchment with comid, as we have removed some of its area.
        target_polygon = subbasins_gdf.at[comid, 'geometry']
        revised_poly = target_polygon.difference(new_poly)
        cleaned_poly = fix_polygon(revised_poly)
        subbasins_gdf.at[comid, 'geometry'] = cleaned_poly
        # Update the area of the clipped unit catchment and insert value into GeoDataFrame
        area = calc_area(cleaned_poly)
        subbasins_gdf.at[comid, 'unitarea'] = round(area, 1)

    if PLOTS:
        plot_basins(subbasins_gdf, points_gdf, 'after')

    # After splitting the unit catchments, we can remove the most downstream node
    terminal_comid = new_nodes[terminal_node_id]
    subbasins_gdf = subbasins_gdf.drop(terminal_comid)

    # Split the rivers data, based on the results above.
    #   FOr most unit catchments, the river will not change.
    #   But for the set of nodes and comids in our dictionary new_nodes, we need to do a clipping operation
    myrivers_gdf = rivers_gdf.loc[B]

    # Drop all columns except a couple that we want to keep. Why? Because after splitting the river reach segments,
    #  this data is not going to be accurate anymore.
    columns_to_keep = ['lengthkm', 'geometry']
    myrivers_gdf = myrivers_gdf[columns_to_keep]

    # First, insert a copy of the river reach into the new nodes we inserted.
    for node, comid in new_nodes.items():
        row_to_copy = myrivers_gdf.loc[comid]
        row_to_copy_df = row_to_copy.to_frame().T
        row_to_copy_df.index = [node]

        myrivers_gdf = pd.concat([myrivers_gdf, row_to_copy_df])

    # Now in both the new nodes and the unit catchment, clip the river using the polygon boundary.
    for node, comid in new_nodes.items():
        for id in [node, comid]:
            if id != terminal_comid:
                river_polyline = myrivers_gdf.at[id, 'geometry']
                basin_polygon = subbasins_gdf.at[id, 'geometry']
                clipped_line = river_polyline.intersection(basin_polygon)
                myrivers_gdf.at[id, 'geometry'] = clipped_line
                lengthkm = calc_length(clipped_line)
                myrivers_gdf.at[id, 'lengthkm'] = round(lengthkm, 1)

    myrivers_gdf = myrivers_gdf.drop(terminal_comid)

    # Export the geodata for (1) subbasins, (2) outlets, and (3) rivers
    fname = f"{OUTPUT_DIR}/subbasins.{OUTPUT_EXT}"
    write_geodata(subbasins_gdf, fname)

    # Extract the outlet data and write it to disk
    outlets_gdf = subbasins_gdf.copy()
    outlets_gdf['geometry'] = outlets_gdf.apply(lambda row: Point(row['lng'], row['lat']), axis=1)
    fname = f"{OUTPUT_DIR}/outlets.{OUTPUT_EXT}"
    write_geodata(outlets_gdf, fname)

    # Write the rivers data to disk.
    myrivers_gdf.reset_index(inplace=True)
    myrivers_gdf.rename(columns={'index': 'comid'}, inplace=True)
    fname = f"{OUTPUT_DIR}/rivers.{OUTPUT_EXT}"
    write_geodata(myrivers_gdf, fname)

    if VERBOSE: print("Ran successfully!")


if __name__ == "__main__":
    delineate2()
