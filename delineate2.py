"""
Fast watershed delineation for anywhere on the Earth's surface.
Creates watersheds for outlet points in an input CSV file.

Interactive demo: https://mghydro.com/watersheds

Usage: Carefully edit the file config2.py, then run delineator2.py
See README for more detailed instructions.

"""

import warnings
import time
import numpy as np
import pickle
import pandas as pd
import os
import geopandas as gpd
import re
from shapely.geometry import Point, Polygon, box
import shapely.ops
from shapely.wkt import loads
import sigfig  # for formatting numbers to significant digits
from py.fast_dissolve import dissolve_geopandas, fill_geopandas
import pyproj
from functools import partial
from config2 import *
from py.mapper import make_map, create_folder_if_not_exists
import py.merit_detailed
import matplotlib.pyplot as plt


warnings.simplefilter(action='ignore', category=FutureWarning)

#gpd.options.use_pygeos = True

# The WGS84 projection string, used in a few places
PROJ_WGS84 = 'EPSG:4326'


def validate(gages_df: pd.DataFrame) -> bool:
    """
    After we have read in the user's input CSV file with their desired delineation locations
    (I refer to these as gages as that was my original use case), check whether the input
    data is valid -- required columns are there (at a minimum id, lat, lng) and that the
    data types are correct and values are in the appropriate range (i.e. lat between -90 and 90).

    returns: True (valid) or throws an Exception

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

    return True


def get_area(poly: Polygon) -> float:
    """
    Projects a Shapely polygon in raw lat, lng coordinates, and calculates its area
    Args:
        poly: Shapely polygon
    :return: area in km²
    """
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

def load_megabasins() -> gpd.GeoDataFrame:
    """
    Reads the megabasin data from disk and returns a GeoDataFrame.
    Will get the data from a shapefile or a pickle file, if it exists.
    If the .pkl does not exist, create it for faster processing in the future.
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
    Finds out what Pfafstetter Level 2 megabasin a point is in.
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

def delineate2():
    """
    *** MAIN watershed delineation routine ***
    Make sure to set the variables in `config.py` before running.

    Reads a list of outlet points from a .csv file,
    then finds their watersheds or drainage basins,
    using hybrid of vector- and raster-based methods.

    Outputs geodata (.shp, .gpkg, etc.) and optionally a CSV file with a summary of results
    (including the watershed id, names, and areas).

    Optionally creates an HTML page with a handy map viewer to review the results.
    """

    # Regular expression used to find numbers so I can round lat, lng coordinates in GeoJSON files to make them smaller
    simpledec = re.compile(r"\d*\.\d+")

    def mround(match):
        # For rounding the coordinates in GeoJSON files to make them smaller
        return "{:.5f}".format(float(match.group()))

    def addnode(B, node):
        """"
        Recursive function to assemble the list of upstream unit catchments
        B is a list of the unit catchments that make up our watershed.
        List items are `comid`s, unique identifiers of each unit catchment.
        """
        # first, append the node to the basin
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

    def plot_basins(suffix: str):
        """
        Makes a plot of the unit catchments that are in the watershed

        It is all the upstream unit catchments, and the *split* terminal unit catchment.

        """
        # subbasins_gdf.plot(column='area', edgecolor='gray', legend=True)
        [fig, ax] = plt.subplots(1, 1, figsize=(10, 8))

        # Plot each unit catchment with a different color
        for x in subbasins_gdf.index:
            color = np.random.rand(3, )
            subbasins_gdf.loc[[x]].plot(facecolor=color, edgecolor=color, alpha=0.5, ax=ax)

        # Plot the gage point
        gages_joined.iloc[[i]].plot(ax=ax, c='red', edgecolors='black')

        if suffix == "post" and HIGH_RES:
            plt.scatter(x=lng_snap, y=lat_snap, c='cyan', edgecolors='black')
            plt.title(f"Showing the {len(subbasins_gdf)-1} upstream unit catchments and split terminal unit catchment")
        else:
            plt.title(f"Found {len(subbasins_gdf)} unit catchments for watershed id = {wid}")

        plt.savefig(f"plots/{wid}_vector_unit_catchments_{suffix}.png")
        plt.close(fig)

    def write_geodata():
        # SAVE the WATERSHED to disk as a GeoJSON file or a shapefile
        if VERBOSE: print(f' Writing WATERSHED geodata for watershed {wid}')
        outfile = f"{OUTPUT_DIR}/wshed_{wid}.{OUTPUT_EXT}"

        # This line rounds all the vertices to fewer digits. For text-like formats GeoJSON or KML, makes smaller
        # files with minimal loss of precision. For other formats (shp, gpkg), doesn't change file size
        if OUTPUT_EXT.lower() in ['geojson', 'kml']:
            mybasin_gdf.geometry = mybasin_gdf.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=UserWarning)
            mybasin_gdf.to_file(outfile)

        # RIVERS
        if VERBOSE: print(f' Writing RIVERS geodata for watershed {wid}')
        outfile = f"{OUTPUT_DIR}/rivers_{wid}.{OUTPUT_EXT}"

        # TODO: Clip the little downstream portion of the last river reach? 
        myrivers_gdf = rivers_gdf.loc[B]

        # This line rounds all the vertices to fewer digits. For text-like formats GeoJSON or KML, makes smaller
        # files with minimal loss of precision. For other formats (shp, gpkg), doesn't change file size
        if OUTPUT_EXT.lower() in ['geojson', 'kml']:
            myrivers_gdf.geometry = rivers_gdf.geometry.apply(
                lambda x: loads(re.sub(simpledec, mround, x.wkt)))

        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=UserWarning)
            myrivers_gdf.to_file(outfile)

    def write_map_data():
        """
        Creates the .js files for the interactive .html map.
        We need to write a second, slightly different version of the GeoJSON files,
        because we need it in a .js file assigned to a variable, to avoid cross-origin restrictions
        of modern web browsers. The workaround for this issue would be to use a simple webserver
        on the local machine instead of just opening the .html file, but that seemed too complicated.
        """

        # WATERSHED Boundary Polygon
        watershed_js = f"{MAP_FOLDER}/{wid}.js"
        with open(watershed_js, 'w') as f:
            s = f"gage_coords = [{gages_df.loc[wid, 'lat']}, {gages_df.loc[wid, 'lng']}];\n"
            f.write(s)
            s = f"snapped_coords = [{gages_df.loc[wid, 'lat_snap']}, {gages_df.loc[wid, 'lng_snap']}];\n"
            f.write(s)

            f.write("basin = ")
            f.write(mybasin_gdf.to_json())

        # RIVERS polylines
        if MAP_RIVERS:
            myrivers_gdf = rivers_gdf.loc[B]

            # Keep only the fields lengthkm and order
            myrivers_gdf = myrivers_gdf[['lengthkm', 'order', 'geometry']]

            # Filter out the little headwater streams in large watersheds.
            max_order = myrivers_gdf.order.max()
            min_order = max_order - NUM_STREAM_ORDERS
            # Drop rows where order < min_order
            myrivers_gdf = myrivers_gdf[myrivers_gdf.order >= min_order]
            myrivers_gdf = myrivers_gdf.round(1)

            myrivers_gdf.geometry = myrivers_gdf.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

            if SIMPLIFY:
                myrivers_gdf.geometry = myrivers_gdf.geometry.simplify(tolerance=SIMPLIFY_TOLERANCE)

            rivers_js = f"{MAP_FOLDER}/{wid}_rivers.js"
            with open(rivers_js, 'w') as f:
                f.write("rivers = ")
                f.write(myrivers_gdf.to_json())

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
    gages_df.set_index('id', inplace=True)
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

    # Perform a Spatial join on gages (points) and the unit catchments (polygons)
    # to find the corresponding unit catchment for each gage
    # Adds the fields COMID and unitarea
    if VERBOSE: print(f"Performing spatial join on {num_gages} outlet points in basin #{megabasin}")
    gages_joined = gpd.overlay(gages_df, catchments_gdf, how="intersection")
    gages_joined.rename(columns={"index_right": "COMID"}, inplace=True)

    # For any gages for which we could not find a unit catchment, add them to failed
    gages_matched = gages_joined['id'].tolist()
    gage_basin_ids = gages_df['id'].tolist()
    for wid in gage_basin_ids:
        if wid not in gages_matched:
            failed[wid] = f"Could not assign to a unit catchment in Level 2 basin #{megabasin}"

        # Revise the number of gages in the basin, based on those which have a matching COMID
        num_gages = len(gages_joined)

        # Iterate over the gages and assemble the watershed
        for i in range(0, num_gages):

            gages_counter += 1
            if VERBOSE: print(f"\n* Delineating watershed {gages_counter} of {n_gages}, with outlet id = {wid}")

            # Let wid be the watershed ID. Get the lat, lng coords of the gage.
            wid = gages_joined['id'].iloc[i]
            lat = gages_joined['lat'].iloc[i]
            lng = gages_joined['lng'].iloc[i]

            # The terminal comid is the unit catchment that contains (overlaps) the outlet point
            terminal_comid = gages_joined['COMID'].iloc[i]

            # Get the upstream area of the unit catchment we found, according to MERIT-Basins
            up_area = rivers_gdf.loc[terminal_comid].uparea

            # Let B be the list of unit catchments (and river reaches) that are in the basin
            B = []

            # Add the first node, and the rest will be added recursively
            addnode(B, terminal_comid)
            if VERBOSE: print(f"  found {len(B)} unit catchments in the watershed")

            subbasins_gdf = catchments_gdf.loc[B]


            if VERBOSE: print("Performing detailed raster-based delineation for "
                              "the downstream portion of the watershed")
            # Let split_catchment_poly be the polygon of the terminal unit catchment
            assert terminal_comid == B[0]
            catchment_poly = catchments_hires_gdf.loc[terminal_comid].geometry
            bSingleCatchment = len(B) == 1
            split_catchment_poly, lat_snap, lng_snap = py.merit_detailed.split_catchment(wid, megabasin, lat, lng,
                                                                                         catchment_poly,
                                                                                         bSingleCatchment)
            if split_catchment_poly is None:
                failed[wid] = "An error occured in pysheds detailed delineation."
                continue
            else:
                # Create a temporary GeoDataFrame to create the geometry of the split catchment polygon
                # This is just a shortcut method to transfer it to our watershed's subbasins GDF
                split_gdf = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[split_catchment_poly])
                split_geom = split_gdf.loc[0, 'geometry']
                subbasins_gdf.loc[terminal_comid, 'geometry'] = split_geom

            # Make a plot of the watersheds unit catchments AFTER we have split the terminal unit catchment
            if PLOTS:
                plot_basins("post")

            # New feature 2024-04-21: Return the unit catchments *without* merging and dissolving
            if DISSOLVE:
                if VERBOSE: print("Dissolving...")
                # mybasin_gs is a GeoPandas GeoSeries

                mybasin_gs = dissolve_geopandas(subbasins_gdf)

                if TIMER:
                    # End timer (do not include writing output)
                    end_time = time.time()
                    gages_df.at[wid, 'time'] = sigfig.round(end_time - start_time, 4)

                if FILL:
                    # Fill donut holes in the watershed polygon
                    # Recall we asked the user for the fill threshold in terms of number of pixels
                    PIXEL_AREA = 0.000000695  # Constant for the area of a single pixel in MERIT-Hydro, in decimal degrees
                    area_max = FILL_THRESHOLD * PIXEL_AREA
                    mybasin_gs = fill_geopandas(mybasin_gs, area_max=area_max)

                # Let `mybasin_gdf` be a GeoPandas DataFrame with the geometry, and the id and area of our watershed
                mybasin_gdf = gpd.GeoDataFrame(geometry=mybasin_gs)
                mybasin_gdf['id'] = wid
                basin_poly = mybasin_gdf.geometry.values[0]
                up_area = get_area(basin_poly)
                if bNames:
                    mybasin_gdf['name'] = gages_df.loc[wid, 'name']

                # Add the upstream area of the delineated watershed to the DataFrame
                mybasin_gdf['area_calc'] = up_area
                gages_df.at[wid, 'area_calc'] = up_area

                # If the user provided an a basin area, calculate the percent difference in areas
                if bAreas:
                    area_reported = gages_df.loc[wid].area_reported
                    mybasin_gdf['area_reported'] = area_reported
                    perc_diff = sigfig.round((up_area - area_reported) / area_reported * 100, 2)
                    gages_df.at[wid, 'perc_diff'] = perc_diff


            else:
                mybasin_gdf = subbasins_gdf.copy()


                snapped_outlet = rivers_gdf.loc[terminal_comid].geometry.coords[0]
                lat_snap = snapped_outlet[1]
                lng_snap = snapped_outlet[0]

            # Get the (approx.) snap distance
            geod = pyproj.Geod(ellps='WGS84')
            snap_dist = geod.inv(lng, lat, lng_snap, lat_snap)[2]
            gages_df.at[wid, 'snap_dist'] = sigfig.round(snap_dist, 2)
            gages_df.at[wid, 'lat_snap'] = round(lat_snap, 3)
            gages_df.at[wid, 'lng_snap'] = round(lng_snap, 3)

            if OUTPUT_EXT != "":
                write_geodata()

            # Create the HTML Viewer Map?
            if MAKE_MAP:
                write_map_data()

    # CREATE OUTPUT.CSV, a data table of the outputs
    # id, status (hi, low, failed), name, area_reported, area_calculated
    if OUTPUT_CSV:
        if OUTPUT_FNAME == '':
            output_fname = f"{OUTPUT_DIR}/OUTPUT.csv"
        else:
            output_fname = OUTPUT_FNAME

        gages_df.to_csv(output_fname)

    # FAILED.csv: If there were any failures, write this to a separate CSV file
    if len(failed) > 0:
        print(f"### FAILED to find watersheds for {len(failed)} locations. Check FAILED.csv for info.")

        failfile = f"{OUTPUT_DIR}/FAILED.csv"
        with open(failfile, 'w') as f:
            f.write("ID, EXPLANATION\n")
            for k, v in failed.items():
                f.write(f'{k},"{v}"\n')

    # If the user wants the browser map, make it
    if MAKE_MAP:
        if VERBOSE: print("* Creating viewer.html *")
        make_map(gages_df)

    # Finished, print a little status message
    if VERBOSE: print(f"It's over! See results in {output_fname}")


def make_folders():
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
    gdf.set_index('COMID', inplace=True)

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


if __name__ == "__main__":
    delineate2()
