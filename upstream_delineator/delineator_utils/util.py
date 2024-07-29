# Functions used by the watershed subbasins delineator script. Matt Heberger, May 2022
from typing import Union
import json
import os
from functools import partial
import networkx
import pyproj
import shapely
from shapely.ops import unary_union
from geopandas import GeoDataFrame
from geopandas.io.file import _read_file_pyogrio
from shapely.geometry import MultiPolygon, Polygon, LineString
import geopandas as gpd
import pandas as pd
import re
import pickle
import warnings
import matplotlib.pyplot as plt
from numpy import random

from upstream_delineator import config

# The WGS84 projection string, used in a few places
PROJ_WGS84 = 'EPSG:4326'
CATCHMENT_PATH = os.getenv("CATCHMENT_PATH")
assert CATCHMENT_PATH
RIVER_PATH = os.getenv("RIVER_PATH")
assert RIVER_PATH
MEGABASINS_PATH = os.getenv("MEGABASINS_PATH")
assert MEGABASINS_PATH


# Regular expression used to find numbers so I can round lat, lng coordinates in GeoJSON files to make them smaller
simpledec = re.compile(r"\d*\.\d+")


def get_largest(input_poly: Union[MultiPolygon, Polygon]) -> Polygon:
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
    required_cols = ['id', 'lat', 'lng', 'is_outlet']
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

    # Check that `is_outlet` is boolean
    is_outlet_type = gages_df['is_outlet'][0]
    if is_outlet_type != 'bool':
        if is_outlet_type == 'O':
            vals = gages_df['id'].unique()
            if True not in vals:
                raise Exception("The field `is_outlet` must be boolean (true/false)")
    else:
        raise Exception("The field `is_outlet` must be boolean (true/false)")

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
            pyproj.Proj(init=PROJ_WGS84),
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
            pyproj.Proj(init=PROJ_WGS84),
            pyproj.Proj(
                proj='aea',
                lat_1=line.bounds[1],
                lat_2=line.bounds[3]
            )
        ),
        line)

    # Get the area in m^2
    return projected_line.length / 1e3


def load_megabasins(bounds: tuple[float]) -> gpd.GeoDataFrame:
    """
    Reads the "megabasin" data from disk and returns a GeoDataFrame.
    The program uses this data to determine what dataset is needed for analyses.
    I refer to the MERIT-Basins Pfafstetter Level 2 basins as megabasins.

    This function gets the data from a flatgeobuf file or a pickle file, if it exists.
    If the .pkl does not exist, create it for faster processing in the future,
    since reading flatgeobuf files is slow.
    """
    megabasins_gdf = load_gdf("megabasins", True, bounds)

    # The CRS string in the flatgeobuf file is EPSG 4326 but does not match verbatim, so set it here
    megabasins_gdf.to_crs(PROJ_WGS84, inplace=True)

    # Check that data is well-formed 
    ((11 <= megabasins_gdf.BASIN) & (megabasins_gdf.BASIN <= 91)).all()

    return megabasins_gdf


def get_megabasins(points_gdf: GeoDataFrame) -> dict:
    """
    Finds out what Pfafstetter Level 2 "megabasins" the outlet points are in
    Arguments:
        GeoDataFrame containing a set of points
    Returns:
        A dictionary, keys are the unique megabasins: integers from 11 to 91.
        Values are lists of outlet points (type is variable, whatever the user entered in the input CSV file)
    """
    all_points = unary_union(points_gdf['geometry'])
    megabasins_gdf = load_megabasins(all_points.bounds)
    if config.get("VERBOSE"): print("Finding out which Pfafstetter Level 2 'megabasin' your outlets are in")

    # Overlay the gage points on the Level 2 Basins polygons to find out which
    # PFAF_2 basin each point falls inside of, using a spatial join
    gages_basins_join = gpd.overlay(points_gdf, megabasins_gdf, how="intersection")

    excluded_points = points_gdf.loc[~points_gdf['id'].isin(gages_basins_join['id'])]
    if not excluded_points.empty:
        raise ValueError(f'These points do not fall inside any of the continental-scale megabasins.\n{excluded_points}')
    
    # Needed to set this option in order to avoid a warning message in Geopandas.
    # https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
    pd.options.mode.chained_assignment = None

    # Get a list of the DISTINCT Level 2 basins, and a count of how many gages in each.
    basins_dict = gages_basins_join.groupby('BASIN')['id'].apply(list).to_dict()

    return basins_dict, megabasins_gdf.set_index('BASIN')


def make_folders():
    """
    This function makes sure that the folders that the user specified
    in config.py are existant. If they are not, it tries to create them.
    If it cannot find the folder, and it cannot create the folder, the
    function will raise an error.

    :return: Nothing, but throws an error if it fails.
    """
    # Check that the OUTPUT directories are there. If not, try to create them.
    folders = [config.get("OUTPUT_DIR"), config.get("PLOTS_DIR"), config.get("PICKLE_DIR")]
    for folder in folders:
        if folder == '':
            continue
        else:
            folder_exists = create_folder_if_not_exists(folder)
        if not folder_exists:
            raise Exception(f"Could not create folder `{folder}`. Stopping")


def get_pickle_filename(geotype: str, bounds: tuple, high_resolution: bool) -> str:
    """Simple function to get the standard filename for the pickle files used by this project.
    The filenames look like this:
       config.get("PICKLE_DIR")/catchments_##_##_##_##_hires.pkl
       config.get("PICKLE_DIR")/catchments_##_##_##_##_lores.pkl

       config.get("PICKLE_DIR")/rivers_##_##_##_##_hires.pkl
       config.get("PICKLE_DIR")/rivers_##_##_##_##_lores.pkl

    where ## is a coord in bounds

    """

    if high_resolution:
        resolution_str = 'hires'
    else:
        resolution_str = 'lores'
    bounds_str = '_'.join(map(str, bounds))
    fname = f'{config.get("PICKLE_DIR")}/{geotype}_{bounds_str}_{resolution_str}.pkl'
    return fname


def load_gdf(geotype: str, high_resolution: bool, bounds: tuple[float]) -> gpd.GeoDataFrame:
    """
    Returns the unit catchments vector polygon dataset as a GeoDataFrame
    Gets the data from the MERIT-Basins flatgeobuf file the first time,
    and after that from a saved .pkl file on disk.
    Uses some global parameters from config.py

    :param geotype: either "catchments" or "rivers" depending on which one we want to open.
    :param high_resolution: True to load the standard (high-resolution) file,
      False to load the low-resolution version (for faster processing, slightly less accurate results)
    :param bounds: tuple of floats representing a bounding box.
      
    :return: a GeoPandas GeoDataFrame

    """

    # First, check for the presence of a pickle file
    if config.get("PICKLE_DIR") != '':
        pickle_fname = get_pickle_filename(geotype, bounds, high_resolution)
        if os.path.isfile(pickle_fname):
            if config.get("VERBOSE"): print(f"Loading BASIN catchment data from pickle file.")
            gdf = pickle.load(open(pickle_fname, "rb"))
            return gdf

    if geotype == "catchments":
        gis_path = CATCHMENT_PATH
    elif geotype == "rivers":
        gis_path = RIVER_PATH
    elif geotype == "megabasins":
        gis_path = MEGABASINS_PATH

    if config.get("VERBOSE"): print(f"Reading geodata in {gis_path}")
    # use _read_file_pygrio instead of gpd.read_file b/c it's performing an unneeded check that causes a 403 error 
    gdf = _read_file_pyogrio(gis_path, bbox=bounds)
    # This line is necessary because some of the gis_paths provided by reachhydro.com do not include .prj files
    gdf.set_crs(PROJ_WGS84, inplace=True, allow_override=True)

    # Before we exit, save the GeoDataFrame as a pickle file, for future speedups!
    save_pickle(geotype, gdf, high_resolution, bounds)

    return gdf


def save_pickle(geotype: str, gdf: gpd.GeoDataFrame, high_resolution: bool, bounds: tuple):
    # If we loaded the catchments from a flatgeobuf, save the gdf to a pickle file for future speedup
    if config.get("PICKLE_DIR") != '':

        # Check whether the GDF has a spatial index.
        # Note: I don't think this is ever necessary. Since version 0.7.0 (March 2020), GeoPandas
        # creates a spatial index by default.
        has_spatial_index = hasattr(gdf, 'sindex')
        if not has_spatial_index:
            gdf.sindex.create_index()

        # Get the standard project filename for the pickle files.
        pickle_fname = get_pickle_filename(geotype, bounds, high_resolution)
        if not os.path.isfile(pickle_fname):
            if config.get("VERBOSE"): print(f"Saving GeoDataFrame to pickle file: {pickle_fname}")
            try:
                pickle.dump(gdf, open(pickle_fname, "wb"))
            except Exception:
                raise Warning("Could not save pickle file to: {pickle_fname}")


def fix_polygon(poly: Union[Polygon, MultiPolygon]):
    """
    When we use the difference() method in Shapely to subtract one polygon from another,
    it's common to end up with small slivers or unwanted geometries around the edges of
    the input polygon due to precision issues or complex boundaries.
    To eliminate these small slivers, we can use a combination of simplification and filtering based on area.

    input: a Shapely Polygon or MultiPolygon
    output: a Shapely Polygon or MultiPolygon that has been fixed to remove small slivers.
    """
    simplify_tolerance = 0.00001
    simplified_poly = poly.simplify(tolerance=simplify_tolerance, preserve_topology=True)

    # Try buffering?
    # buffered = simplified_poly.buffer(simplify_tolerance)

    # Merge the buffered polygons
    # merged = shapely.ops.unary_union(buffered)

    # Remove the buffer by shrinking back
    # final_merged = merged.buffer(-simplify_tolerance)

    # Filter out small slivers by area
    min_area_threshold = 0.00001  # Define your own threshold

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
    if config.get("VERBOSE"): print('Writing geodata to disk')

    # This line rounds all the vertices to fewer digits. For text-like formats GeoJSON or KML, makes smaller
    # files with minimal loss of precision. For other formats (shp, gpkg), it doesn't change file size, so don't bother.
    if config.get("OUTPUT_EXT").lower() in ['geojson', 'kml']:
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

    plt.savefig(f'{config.get("PLOTS_DIR")}/{fname}.png')
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

    filename = f'{config.get("OUTPUT_DIR")}/{prefix}_graph.{file_ext}'

    if file_ext == 'pkl':
        # (1) Python pickle file
        pickle.dump(G, open(filename, "wb"))
    elif file_ext == 'json':
        # (2) JSON file
        data = networkx.node_link_data(G)
        with open(filename, "w") as f:
            json.dump(data, f)
    elif file_ext == 'gml':
        # (3) GML (Graph Modeling Language), a common graph file format.
        networkx.write_gml(G, filename)
    elif file_ext == 'xml':
        # (4) GraphML is an XML-based file format for graphs.
        networkx.write_graphml(G, filename)
    else:
        raise ValueError(f'Unhandled file extension {file_ext}')
