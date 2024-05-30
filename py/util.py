# Functions used by the watershed subbasins delineator script. Matt Heberger, May 2022
import json
import os
from functools import partial
from collections import Counter
import networkx
import pyproj
import shapely
from shapely.geometry import MultiPolygon, Polygon, LineString, MultiLineString
import geopandas as gpd
import pandas as pd
import re
import pickle
import warnings
import matplotlib.pyplot as plt
from subbasins_config import PICKLE_DIR, OUTPUT_DIR, VERBOSE, OUTPUT_EXT, RIVERS_DIR, CATCHMENTS_DIR
from numpy import random

# The WGS84 projection string, used in a few places
PROJ_WGS84 = 'EPSG:4326'

# Regular expression used to find numbers so I can round lat, lng coordinates in GeoJSON files to make them smaller
simpledec = re.compile(r"\d*\.\d+")


def get_largest(input_poly: MultiPolygon or Polygon) -> Polygon:
    """
    Converts a Shapely MultiPolygon to a Shapely Polygon
    For multipart polygons, will only keep the largest polygon
    in terms of area. In my testing, this was usually good enough

    Note: can also do this via PostGIS query... see myqueries_merit.py, query19a and query19b
          Not sure one approach is better than the other. They both seem to work well.
    Args:
        input_poly: A Shapely Polygon or MultiPolygon

    Returns:
        a shapely Polygon
    """
    if input_poly.geom_type == "MultiPolygon":
        areas = []
        polygons = list(input_poly.geoms)

        for poly in polygons:
            areas.append(poly.area)

        max_index = areas.index(max(areas))

        return polygons[max_index]
    else:
        return input_poly


def create_folder_if_not_exists(folder_path: str) -> bool:
    """
    Check if a folder exists at the specified path. If it does not, create it.

    :param folder_path: The path of the folder to check/create.
    :return: True if the folder exists or is created, False otherwise.
    """
    try:
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Created folder: {folder_path}")
        else:
            pass  # Folder is already there

        return True

    except Exception as e:
        print(f"Error creating folder: {e}")
        return False


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


def find_repeats_with_frequency(lst: list) -> dict:
    r"""
    Finds repeated elements in a list with their frequency
    > lst = list('collections')
    ['c', 'o', 'l', 'l', 'e', 'c', 't', 'i', 'o', 'n', 's']
    > find_repeats_with_frequency(lst)
    {'c': 2, 'o': 2, 'l': 2}
    """
    # Count the elements in the list
    element_count = Counter(lst)

    # Filter elements that are repeated
    repeated_elements = {element: count for element, count in element_count.items() if count > 1}
    return repeated_elements


def decrement_repeats_dict(dictionary: dict, key) -> dict:
    """
    Works with the dictionary of repeats with frequency created above.
    Decrements an entry based on its key. Each time it is called,
    the value decreases by one. When the value gets to one, the entry is removed
    (because it is no longer a duplicate or repeated entry!)
    """
    if key in dictionary:
        dictionary[key] -= 1
        if dictionary[key] == 1:
            del dictionary[key]
    return dictionary


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
    return 27
    if VERBOSE: print("Finding out which Pfafstetter Level 2 'megabasin' your outlets are in")
    megabasins_gdf = load_megabasins()

    # Overlay the gage points on the Level 2 Basins polygons to find out which
    # PFAF_2 basin each point falls inside of, using a spatial join
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
        directory = CATCHMENTS_DIR
        shapefile = f"{directory}/cat_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"
    elif geotype == "rivers":
        shapefile = f"{RIVERS_DIR}/riv_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp"

    if not os.path.isfile(shapefile):
        raise Exception(f"Could not find the file: {shapefile}")

    if VERBOSE: print(f"Reading geodata in {shapefile}")
    gdf = gpd.read_file(shapefile)

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
        gdf.geometry = gdf.geometry.apply(lambda x: shapely.wkt.loads(re.sub(simpledec, mround, x.wkt)))

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
        color = random.rand(3, )
        basins_gdf.loc[[x]].plot(facecolor=color, edgecolor=color, alpha=0.5, ax=ax)

    # Plot the gage points
    outlets_gdf.plot(ax=ax, c='red', edgecolors='black')

    plt.savefig(f"plots/{fname}.png")
    plt.close(fig)


def save_network(G: networkx.Graph, prefix: str, file_ext: str):
    """
    Saves the NetworkX graph to disk
    :param G: the graph object
    :param prefix: a string to prepend to the output filename
    :param file_ext: format to save, choose from among pkl, gml, xml, json
    :return:
    """

    # Here are 4 different options for how to save the graph; other options are possible,
    #  especially if NetworkX has a `write_##()` method built-in. See:
    #    https://networkx.org/documentation/stable/reference/readwrite/index.html

    allowed_formats = ['pkl', 'gml', 'xml', 'json']
    if file_ext not in allowed_formats:
        print(f"River network graph not saved. Please choose one of the following formats: "
              f"{', '.join(allowed_formats)}")
        raise Warning("Did not save graph data.")

    filename = f"{OUTPUT_DIR}/{prefix}_graph.{file_ext}"

    match file_ext:
        case 'pkl':
            # (1) Python pickle file
            pickle.dump(G, open(filename, "wb"))
        case 'json':
            # (2) JSON file
            data = networkx.node_link_data(G)
            with open(filename, "w") as f:
                json.dump(data, f)
        case 'gml':
            # (3) GML (Graph Modeling Language), a common graph file format.
            networkx.write_gml(G, filename)
        case 'xml':
            # (4) GraphML is an XML-based file format for graphs.
            networkx.write_graphml(G, filename)


def multilinestring_to_linestring(multi_line):
    """
    Coerce a MultiLineString into a single LineString by connecting the different pieces.

    Parameters:
    multi_line (MultiLineString): The input MultiLineString to be converted.

    Returns:
    LineString: The resulting connected LineString.
    """
    if isinstance(multi_line, LineString):
        return multi_line

    if not isinstance(multi_line, MultiLineString):
        raise ValueError("Input must be a MultiLineString")

    # Extract the individual LineStrings
    line_strings = list(multi_line.geoms)

    # MODIFY to just return the longest one...
    # Nevermind, gave poor results.
    #lengths = [line.length for line in line_strings]
    #max_index = lengths.index(max(lengths))
    #return line_strings[max_index]
    # END

    if len(line_strings) == 0:
        return LineString()

    # Function to sort LineStrings so that they connect end-to-end
    def sort_lines(lines):
        sorted_lines = [lines.pop(0)]
        while lines:
            last_point = sorted_lines[-1].coords[-1]
            found_next = False
            for i, line in enumerate(lines):
                if line.coords[0] == last_point:
                    sorted_lines.append(lines.pop(i))
                    found_next = True
                    break
                elif line.coords[-1] == last_point:
                    sorted_lines.append(LineString(line.coords[::-1]))
                    lines.pop(i)
                    found_next = True
                    break
            if not found_next:
                # If no next line is found, break the loop
                break
        return sorted_lines

    # Find out the number of vertices in each linestring, and if it is 1, remove it
    for line in line_strings:
        if len(line.coords) < 3:
            line_strings.remove(line)

    if len(line_strings) == 1:
        return line_strings[0]


    # Sort the line segments
    sorted_lines = sort_lines(line_strings)

    # Merge the sorted LineStrings into a single LineString
    coords = []
    for line in sorted_lines:
        coords.extend(line.coords)

    # Remove duplicate coordinates that occur at the junctions
    coords = [coords[0]] + [coord for i, coord in enumerate(coords[1:]) if coord != coords[i]]

    # Create the final LineString
    final_line = LineString(coords)

    return final_line
