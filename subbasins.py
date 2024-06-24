r"""
Delineation of watershed subbasins, using data from
MERIT-Basins and MERIT-Hydro.
Created by Matthew Heberger, May 2024.

See README for more detailed instructions.

Quick Start:

First, set parameters in the file config.py.

Run this script from the command line with two required arguments:
$ python subbasins.py outlets.csv testrun

or with a full file path as follows on Windows or Linux:
$ python subbasins.py C:\Users\matt\Desktop\outlets.csv test
$ python subbasins.py /home/files/outlets.csv test

or in Python as follows:
>> from subbasins import delineate
>> delineate('outlets.csv', 'testrun')

"""

# Standard Python libraries. See requirements.txt for recommended versions.
import argparse
import sys
from shapely.geometry import Point
import topojson
import warnings

# My stuff
from config import *  # This file contains a bunch of variables to be set before running this script
from py.consolidate import consolidate_network, show_area_stats
from py.util import make_folders, get_megabasins, load_gdf, plot_basins, calc_area, \
    find_repeated_elements, fix_polygon, calc_length, validate, save_network, \
    write_geodata, PROJ_WGS84  # Contains a bunch of functions
from py.graph_tools import make_river_network, calculate_strahler_stream_order, calculate_shreve_stream_order, \
    prune_node, upstream_nodes  # Functions for working with river network information as a Python NetworkX graph
from py.merit_detailed import split_catchment
from py.plot_network import draw_graph
from py.fast_dissolve import dissolve_geopandas, buffer, close_holes

from os.path import isfile
import geopandas as gpd
import pandas as pd
import networkx as nx

# Shapely throws a bunch of FutureWarnings. Safe to ignore for now, as long as we
# are using a virtual environment, and use the library versions in requirements.txt.
warnings.simplefilter(action='ignore', category=FutureWarning)

PIXEL_AREA = 0.000000695  # Constant for the area of a single pixel in MERIT-Hydro, in decimal degrees
FILL_AREA_MAX = FILL_THRESHOLD * PIXEL_AREA


def get_wshed_rows(df: gpd.GeoDataFrame, outlet):
    """
    Extracts rows of the gages GeoDataFrame for an outlet and any upstream points.

    """
    start_index = df[df['id'] == outlet].index[0]

    # Slice the DataFrame from the start_index to the end
    subset_df = df.iloc[start_index:]

    # Find the index of the first row where is_outlet is True
    indices = list(subset_df[subset_df['is_outlet'] == True].index)

    if len(indices) > 1:
        end_index = indices[1] - 1
        # Select all rows from start_index to end_index (inclusive)
        result_df = subset_df.loc[0:end_index]
    else:
        result_df = subset_df

    return result_df


def delineate(input_csv: str, output_prefix: str):
    """
    Finds the watershed for a set of outlets.
    Make sure to set the variables in `config.py` before running.

    Version #2, where all the points are in a single watershed,
    and we are interested in returning a set of subbasins.

    Reads a list of outlet points from a .csv file.
    The first row should contain the main watershed outlet.
    Subsequent rows of the .csv should contain uptstream points
    within the watershed that we will use to create *subbasin* outlets.

    Outputs geodata (.shp, .gpkg, etc.) and optionally a CSV file with a summary of results
    (including the watershed id, names, and areas).

    """

    # Make sure the user folders from `config.py` exist (for output & plots). If not create them.
    make_folders()

    # Read the outlet points CSV file and put the data into a Pandas DataFrame
    # (I call the outlet points gages, because I usually in delineated watersheds at streamflow gages)
    gages_gdf = make_gages_gdf(input_csv)

    # Create a filtered version with only the *outlets*
    outlets_gdf = gages_gdf[gages_gdf['is_outlet'] == True]

    # Get the megabasin(s) in which the points are located
    # This returns a dictionary. Key: megabasin, Value: list of outlets that are in the megabasin
    # This way, we can process the gages one megabasin at a time, so we only have to read geodata files once.
    gage_basins_dict = get_megabasins(outlets_gdf)

    G = None
    subbasins_gdf = None
    myrivers_gdf = None

    # Iterate over the megabasins
    for megabasin in gage_basins_dict.keys():
        # Iterate over the outlets:
        outlets = gage_basins_dict[megabasin]

        if VERBOSE: print('Reading geodata for unit catchments in megabasin %s' % megabasin)
        catchments_gdf = load_gdf("catchments", megabasin, True)

        # The _network_ data is in the RIVERS file rather than the CATCHMENTS file
        # (this is just how the MERIT-Basins authors did it)
        if VERBOSE: print('Reading geodata for rivers in megabasin %s' % megabasin)
        rivers_gdf = load_gdf("rivers", megabasin, True)
        rivers_gdf.set_index('COMID', inplace=True)
        # We wish to report the outlet point for each subbasin.
        # We can get this information from end point of the river polylines.
        rivers_gdf['end_point'] = rivers_gdf['geometry'].apply(lambda x: x.coords[0])
        rivers_gdf['lng'] = rivers_gdf['end_point'].apply(lambda x: x[0])
        rivers_gdf['lat'] = rivers_gdf['end_point'].apply(lambda x: x[1])

        for outlet in outlets:
            wshed_gages_gdf = get_wshed_rows(gages_gdf, outlet)
            wshed_G, wshed_subbasins_gdf, wshed_rivers_gdf = get_watershed(wshed_gages_gdf, megabasin,
                                                                           catchments_gdf, rivers_gdf)
            # Merge the results with the master
            if G is not None:
                G = nx.compose(G, wshed_G)
                subbasins_gdf = pd.concat([subbasins_gdf, wshed_subbasins_gdf])
                myrivers_gdf = pd.concat([myrivers_gdf, wshed_rivers_gdf])
            else:
                G = wshed_G
                subbasins_gdf = wshed_subbasins_gdf
                myrivers_gdf = wshed_rivers_gdf

    # Finally, write the results to disk
    gages_list = gages_gdf['id'].tolist()
    write_outputs(G, myrivers_gdf, subbasins_gdf, gages_list, output_prefix)

    if NETWORK_DIAGRAMS:
        draw_graph(G, f'plots/{output_prefix}_network_final')

    if VERBOSE:
        print("Ran succesfully!")


def get_watershed(gages_gdf: gpd.GeoDataFrame, megabasin: int, catchments_gdf, rivers_gdf):
    """
    Finds the watershed and subbasins upstream of an outlet, including any intermediate locations
    (such as gages) where the subbasins should be split.

    Inputs:
        gages_gdf: GeoDataFrame with the gages to delineate subbasins for


    :return:
    """

    def addnode(B: list, node_id):
        """"
        Recursive function to assemble the list of upstream unit catchments
        B is a Python List of the unit catchments that make up our watershed.
        B is for BASIN...
        List items are `comid`s, unique identifiers of each unit catchment.

        """
        # first, append the node to the list that is in the Basin, upstream_comids
        B.append(node_id)

        # next, check whether the fields up1, up2, up3, and up4 contain a node ID
        up1 = rivers_gdf['up1'].loc[node_id]
        if up1 != 0:
            addnode(B, up1)

        up2 = rivers_gdf['up2'].loc[node_id]
        if up2 != 0:
            addnode(B, up2)

        up3 = rivers_gdf['up3'].loc[node_id]
        if up3 != 0:
            addnode(B, up3)

        up4 = rivers_gdf['up4'].loc[node_id]
        if up4 != 0:
            addnode(B, up4)

    # Perform an overlay analysis on gages (points) and the unit catchments (polygons)
    # to find the corresponding unit catchment in which each gages is located.
    # Adds the fields COMID and unitarea to `gages_gdf`

    gages_list = gages_gdf['id'].tolist()  # Get the list before doing the join.
    num_gages = len(gages_list)
    if VERBOSE:
        print(f"Performing overlay analysis on {num_gages} outlet points in basin #{megabasin}")
    catchments_gdf.reset_index(inplace=True)
    gages_gdf = gpd.overlay(gages_gdf, catchments_gdf, how="intersection", make_valid=True)
    gages_gdf.set_index('id', inplace=True)
    gages_gdf.set_crs(crs=PROJ_WGS84)

    # For any gages for which we could not find a unit catchment, add issue a warning
    # Basically checking which rows do not appear after doing the overlay
    gages_matched_list = gages_gdf.index.tolist()
    for id in gages_list:
        if id not in gages_matched_list:
            raise Warning(f"Could not assign to a unit catchment to gage with id {id}")

    # First, let us find the set of unit catchments upstream of the outlet.
    terminal_node_id = gages_gdf.index[0]

    # The terminal comid is the unit catchment that contains (overlaps) the outlet point
    terminal_comid = gages_gdf['COMID'].iloc[0]

    # Let upstream_comids be the list of unit catchments (and river reaches) that are in the basin
    upstream_comids = []

    # Add the first node, and the rest will be added recursively

    addnode(upstream_comids, terminal_comid)

    # Next, check that all the other points provided by the user are
    # fall in unit catchments that are in this set. Otherwise, it means they are not upstream
    # of the outlet, and therefore we cannot process them to get expected results.
    gage_list = gages_gdf.index.tolist()
    for id in gage_list:
        comid = gages_gdf.at[id, 'COMID']
        if comid not in upstream_comids:
            gages_gdf.drop(id, inplace=True)
            raise Warning(f"The point with id = {id} is not contained in the watershed of the first point.")

    # subbasins_gdf is the set of unit catchments in our watershed. This will ultimately become our output
    catchments_gdf.set_index('COMID', inplace=True)
    subbasins_gdf = catchments_gdf.loc[upstream_comids]
    # Add lat, lng, and NextDownID to subbasins_gdf.
    subbasins_gdf = subbasins_gdf.join(rivers_gdf[['lat', 'lng', 'NextDownID']])
    # Re-name the NextDownID field, and make sure it is an integer
    subbasins_gdf.rename(columns={'NextDownID': 'nextdown'}, inplace=True)
    subbasins_gdf['nextdown'] = subbasins_gdf['nextdown'].astype(int)
    subbasins_gdf['custom'] = False  # Adds a column that shows whether a subbasin is connected to a custom pour point

    if PLOTS: plot_basins(subbasins_gdf, gages_gdf, 'before')

    # For debugging mostly
    # G = make_river_network(subbasins_gdf, terminal_comid)
    # if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{OUTPUT_PREFIX}_network_before')

    # Create two copies of river network data!
    # `allrivers_gdf` will contain all the available polylines in the watershed, a nice
    # visual representation of the river network.
    # `myrivers_gdf` will contain the topologically correct network where there is exactly
    # one polyline per unit catchment.
    # TODO: Will need to update this feature if we want to keep it
    # if OUTPUT_ALL_RIVERS:
    #    allrivers_gdf = rivers_gdf.loc[upstream_comids]
    #    fname = f"{OUTPUT_DIR}/{output_prefix}_allrivers.{OUTPUT_EXT}"
    #    write_geodata(allrivers_gdf, fname)
    #    del allrivers_gdf

    # With this version, we will try to create a topologically correct visual representation
    # of the river network, where there is a 1:1 mapping of river reaches to subbasins.
    myrivers_gdf = rivers_gdf.loc[upstream_comids]

    # For now, put the split polygon geometry into a field in `gages_gdf`
    gages_gdf['polygon'] = None
    gages_gdf['polygon_area'] = 0

    # Iterate over the gages, and run `split_catchment()` for every gage
    for gage_id in gages_gdf.index:
        comid = gages_gdf.at[gage_id, 'COMID']
        lat = gages_gdf.at[gage_id, 'lat']
        lng = gages_gdf.at[gage_id, 'lng']
        catchment_poly = subbasins_gdf.loc[comid, 'geometry']

        is_leaf = rivers_gdf.at[comid, 'up1'] == 0  # A leaf is an unit catchment with no upstream neighbor.
        # SPLIT
        node_poly, lat_snap, lng_snap = split_catchment(gage_id, megabasin, lat, lng, catchment_poly, is_leaf)
        gages_gdf.at[gage_id, 'polygon'] = node_poly
        gages_gdf.at[gage_id, 'lat_snap'] = lat_snap
        gages_gdf.at[gage_id, 'lng_snap'] = lng_snap

        # Find the area of the clipped unit catchment and insert value into GeoDataFrame
        area = calc_area(node_poly)
        gages_gdf.at[gage_id, 'polygon_area'] = round(area, 1)

    subbasins_gdf = update_split_catchment_geo(gage_id, gages_gdf, myrivers_gdf, rivers_gdf, subbasins_gdf)

    # Now, we no longer need the downstream portion of the terminal unit catchment
    # so remove its row from the subbasins GeoDataFrame
    subbasins_gdf = subbasins_gdf.drop(terminal_comid)
    subbasins_gdf.at[terminal_node_id, 'nextdown'] = 0

    # Create a NETWORK GRAPH of the river basin.
    if VERBOSE: print("Creating Network GRAPH")
    G = make_river_network(subbasins_gdf, terminal_comid)
    if VERBOSE: show_area_stats(G)
    G = calculate_shreve_stream_order(G)
    G = calculate_strahler_stream_order(G)

    # Add an attribute to the Graph, identifying the 'custom' nodes added by the script
    for gage in gage_list:
        G.nodes[gage]['custom'] = True

    # Draw the network before consolidating? Mostly useful for debugging.
    # if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{OUTPUT_PREFIX}_premerge')

    # CHECK FOR Null Geometries.
    # If the user has placed one of their points very close to an existing basin outlet,
    # we can get a null geometry. This then creates problems. So remove it.

    null_nodes = subbasins_gdf.index[subbasins_gdf['geometry'].is_empty].to_list()
    subbasins_gdf.drop(null_nodes, inplace=True)

    # Remove the null nodes from the graph
    for node in null_nodes:
        prune_node(G, node)

    if len(null_nodes) > 0:
        if VERBOSE: print(f"Pruned {len(null_nodes)} empty unit catchments from the river network.")
        print(null_nodes)
        # if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{OUTPUT_PREFIX}_network_pruned')

    null_rivers = myrivers_gdf.index[myrivers_gdf['geometry'].is_empty].tolist()
    myrivers_gdf.drop(null_rivers, inplace=True)

    # SIMPLIFY the geodata?
    # Can speed up subsequent geo operations, but sometimes also causes problems!
    if SIMPLIFY:
        topo = topojson.Topology(subbasins_gdf, prequantize=False)
        subbasins_gdf = topo.toposimplify(SIMPLIFY_TOLERANCE).to_gdf()
        # The simplification often creates topology errors (ironic!), so try my trick to fix.
        subbasins_gdf.geometry = subbasins_gdf.geometry.apply(lambda p: buffer(p))
        topo = topojson.Topology(myrivers_gdf, prequantize=False)
        myrivers_gdf = topo.toposimplify(SIMPLIFY_TOLERANCE).to_gdf()

    # Add attributes for subbasin area and river reach length to the Graph
    # We will keep track of these quantities as we delete and merge nodes
    for node in list(G.nodes):
        area = subbasins_gdf.at[node, 'unitarea']
        G.nodes[node]['area'] = round(area, 1)
        try:
            length = myrivers_gdf.at[node, 'lengthkm']
        except Exception:
            length = 0
        G.nodes[node]['length'] = round(length, 1)

    # If the user wants larger subbasins, we can merge adjacent unit catchments.
    # I call this "consolidating the river network."
    if CONSOLIDATE:
        G, MERGES, rivers2merge, rivers2delete = consolidate_network(G, threshold_area=MAX_AREA)
        # if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{OUTPUT_PREFIX}_network_final')
        subbasins_gdf['target'] = subbasins_gdf.index

        # Dissolve unit catchments based on information in MERGES.
        for node, target in MERGES.items():
            subbasins_gdf.at[node, 'target'] = target

        agg = {
            'unitarea': 'sum',
            'lat': 'last',
            'lng': 'last',
        }

        if VERBOSE: print("Dissolving geometries. Can be slow. Please wait...")
        subbasins_gdf = subbasins_gdf.dissolve(by="target", aggfunc=agg)
        subbasins_gdf.geometry = subbasins_gdf.geometry.apply(lambda p: buffer(p))
        subbasins_gdf.geometry = subbasins_gdf.geometry.apply(lambda p: close_holes(p, 0))
        subbasins_gdf.reset_index(inplace=True)
        subbasins_gdf.rename(columns={'target': 'comid'}, inplace=True)
        subbasins_gdf.set_index('comid', inplace=True)

        # After the dissolve operation, no way to preserve correct information in column 'nextdown'
        # But it is present in the Graph, so update the GeoDataFrame subbasins_gdf with that information
        # Add column `nextdownid` and the stream orders based on data in the graph
        subbasins_gdf['nextdown'] = 0

        for idx in subbasins_gdf.index:
            try:
                nextdown = list(G.successors(idx))[0]
            except Exception:
                nextdown = 0
            subbasins_gdf.at[idx, 'nextdown'] = nextdown

    # Round all the areas to four decimals (just for appearances)
    subbasins_gdf['unitarea'] = subbasins_gdf['unitarea'].round(1)
    subbasins_gdf['lat'] = subbasins_gdf['lat'].round(4)
    subbasins_gdf['lng'] = subbasins_gdf['lng'].round(4)

    if CONSOLIDATE:
        # Now handle the river reaches, deleting some rows, dissolving others.
        myrivers_gdf.drop(rivers2delete, inplace=True)

        if len(rivers2merge) > 0:
            myrivers_gdf['target'] = myrivers_gdf.index
            for target, node_list in rivers2merge.items():
                for node in node_list:
                    myrivers_gdf.at[node, 'target'] = target

            agg = {'lengthkm': 'sum'}
            myrivers_gdf = myrivers_gdf.dissolve(by="target", aggfunc=agg)
            myrivers_gdf.reset_index(inplace=True)
            myrivers_gdf.rename(columns={'target': 'comid'}, inplace=True)
            myrivers_gdf.set_index('comid', inplace=True)

    # We can now delete the river reach segment that belonged to the downstream portion of the
    # terminal unit catchment, as we no longer need it.
    try:
        myrivers_gdf.drop(terminal_comid, inplace=True)
    except Exception:
        pass

    # Add the fields `nextdown` and the stream orders to the rivers.
    for idx in myrivers_gdf.index:
        try:
            nextdown = list(G.successors(idx))[0]
        except Exception:
            nextdown = 0
        myrivers_gdf.at[idx, 'nextdown'] = nextdown
        myrivers_gdf.at[idx, 'strahler_order'] = G.nodes[idx]['strahler_order']
        myrivers_gdf.at[idx, 'shreve_order'] = G.nodes[idx]['shreve_order']

    # Update the lat/lng coordinates of the subbasin outlets
    subbasins_gdf['custom'] = False
    subbasins_gdf['strahler_order'] = 0
    subbasins_gdf['shreve_order'] = 0

    for idx in subbasins_gdf.index:
        subbasins_gdf.at[idx, 'strahler_order'] = G.nodes[idx]['strahler_order']
        subbasins_gdf.at[idx, 'shreve_order'] = G.nodes[idx]['shreve_order']

        if idx in rivers_gdf.index:
            subbasins_gdf.at[idx, 'lat'] = rivers_gdf.at[idx, 'lat']
            subbasins_gdf.at[idx, 'lng'] = rivers_gdf.at[idx, 'lng']

    for gage in gage_list:
        subbasins_gdf.at[gage, 'custom'] = True

    # Before exporting geodata, make 'comid' a regular column
    try:
        subbasins_gdf.drop(columns=['COMID'], inplace=True)
    except Exception:
        pass

    subbasins_gdf.reset_index(inplace=True)
    subbasins_gdf.rename(columns={'index': 'comid'}, inplace=True)
    myrivers_gdf.reset_index(inplace=True)
    myrivers_gdf.rename(columns={'index': 'comid'}, inplace=True)

    # The rivers data will no longer be accurate, so drop these columns
    cols = ["lengthdir", "sinuosity", "slope", "uparea", "order", "strmDrop_t", "slope_taud", "NextDownID",
            "maxup", "up1", "up2", "up3", "up4", "end_point", "lat", "lng"]
    try:
        myrivers_gdf.drop(columns=cols, inplace=True)
    except Exception:
        pass

    return G, subbasins_gdf, myrivers_gdf


def update_split_catchment_geo(gage_id, gages_gdf, myrivers_gdf, rivers_gdf, subbasins_gdf):
    # Handle the case where we had more than one gage point in a single unit catchment.
    #   (1) figure out their order from upstream to downstream based on split catchment area
    #   (2) clip the polygons as appropriate
    #   (3) insert this data into subbasins_gdf, including polygon geometry and `nextdown` gage id
    # Get a list of unit catchments that contains outlet points (may contain duplicate items)
    comids = gages_gdf['COMID'].tolist()
    # Get a list of unit catchments which contains more than one outlet point
    repeats = find_repeated_elements(comids)
    # Get a list of unit catchments that contain a single outlet point
    singles = [item for item in comids if item not in repeats]
    # Drop columns that we no longer need.
    gages_gdf.drop(columns=['lat', 'lng', 'geometry', 'unitarea'], inplace=True)
    # Rename certain columns to make it identical to `subbasins_gdf` so we can concatenate the rows.
    rnmap = {'lat_snap': 'lat',
             'lng_snap': 'lng',
             'polygon_area': 'unitarea',
             'polygon': 'geometry'
             }
    gages_gdf.rename(columns=rnmap, inplace=True)
    gages_gdf.set_crs(crs=PROJ_WGS84)
    gages_gdf.set_geometry(col="geometry")
    # First, handle the gages where there is only one gage in a unit catchment (standard treatment)
    # The new unit catchments (or nodes in the network) will always be upstream of the unit catchment
    # that we are inserting it into
    #  insert these rows into `subbasins_gdf`
    if len(singles) > 0:
        selected_rows = gages_gdf[gages_gdf['COMID'].isin(singles)]
        selected_rows.set_crs(crs=PROJ_WGS84)  # Just needed to eliminate an annoying warning

        # This creates the dictionary `new_nodes` that maps gage id : unit catchment comid
        new_nodes = selected_rows['COMID'].to_dict()
        selected_rows.drop(columns=['COMID'], inplace=True)

        # Add the column `nextdown` to our temporary GeoDataFrame selected_rows to prevent type problems below
        selected_rows['nextdown'] = -999

        # Here is where we add the newly-created split catchments corresponding to our outlet points.
        subbasins_gdf = pd.concat([subbasins_gdf, selected_rows])

        # Find and fix the 'nextdown' attribute for any nodes that were u/s of terminal catchments
        # into which we just inserted our nodes.
        # Add the column `nextdown`; we will populate it below

        for node, comid in new_nodes.items():
            rows = subbasins_gdf['nextdown'] == comid
            subbasins_gdf.loc[rows, 'nextdown'] = node
            subbasins_gdf.at[node, 'nextdown'] = comid

            # Subtract the polygon geometry of the split catchment from its parent unit catchment
            comid_poly = subbasins_gdf.loc[comid, 'geometry']
            node_poly = subbasins_gdf.loc[node, 'geometry']
            updated_poly = comid_poly.difference(node_poly)
            updated_poly = fix_polygon(updated_poly)
            subbasins_gdf.loc[comid, 'geometry'] = updated_poly

            # Split the river polyline at the new unit catchment boundary
            # First, get the piece in what is left of the old unit catchment
            comid_line = rivers_gdf.at[comid, 'geometry']
            updated_line = comid_line.intersection(updated_poly)
            length = calc_length(updated_line)
            myrivers_gdf.at[comid, 'geometry'] = updated_line
            myrivers_gdf.at[comid, 'lengthkm'] = length

            # Now get the river reach polyline for the new node
            node_line = comid_line.intersection(node_poly)
            length = calc_length(node_line)
            myrivers_gdf.at[node, 'geometry'] = node_line
            myrivers_gdf.at[node, 'lengthkm'] = length
    # Next, handle the case where there are multiple gages in a single unit catchment (special treatment)
    # We need to handle these one unit catchment at a time.
    # I call these nested outlets for want of a better name.
    # (Note: It is not really possible to understand the following code without a picture of what it is doing!!!)
    for comid in repeats:
        # Find all the gages that fall in this unit catchment.
        gages_set = gages_gdf[gages_gdf['COMID'] == comid]
        gages_set.set_crs(crs=PROJ_WGS84)
        # We want to handle these in order from downstream to upstream.
        # This is the same as largest area to smallest area
        gages_set.sort_values(by='unitarea', inplace=True)
        # This one-liner fills in the correct network connection information for nested outlets
        gages_set['nextdown'] = gages_set.index.to_series().shift(-1).fillna(comid)

        gages_set.sort_values(by='unitarea', ascending=False, inplace=True)
        subbasins_gdf = pd.concat([subbasins_gdf, gages_set])

        # Get the (whole) river polyline in the original, unsplit unit catchment, so we can split into pieces
        comid_poly = subbasins_gdf.at[comid, 'geometry']
        comid_line = rivers_gdf.at[comid, 'geometry']

        # Get a list of the gages
        gages = gages_set.index

        # Update the catchment polygon for the base unit catchment by subtracting the *first* clipped catchment
        first_node = gages[0]
        poly1 = gages_set.at[first_node, 'geometry']
        updated_poly = comid_poly.difference(poly1)
        subbasins_gdf.at[comid, 'geometry'] = updated_poly
        area = calc_area(updated_poly)
        subbasins_gdf.at[comid, 'unitarea'] = area

        # Update the river polyline for the base unit catchment
        node_line = comid_line.intersection(updated_poly)
        myrivers_gdf.at[comid, 'geometry'] = node_line
        length = calc_length(node_line)
        myrivers_gdf.at[comid, 'lengthkm'] = length

        # Update the last (most upstream, smallest) unit catchment's polygon and line
        n = len(gages)
        last_node = gages[n - 1]
        last_node_poly = gages_set.at[last_node, 'geometry']
        subbasins_gdf.at[last_node, 'geometry'] = last_node_poly
        area = calc_area(last_node_poly)
        subbasins_gdf.at[last_node, 'unitarea'] = area

        node_line = comid_line.intersection(last_node_poly)
        myrivers_gdf.at[last_node, 'geometry'] = node_line
        length = calc_length(node_line)
        myrivers_gdf.at[last_node, 'lengthkm'] = length

        for i in range(0, n - 1):
            # Handle the unit catchments from largest to smallest
            node_poly = gages_set.at[gages[i], 'geometry']
            nextup_poly = gages_set.at[gages[i + 1], 'geometry']
            updated_poly = node_poly.difference(nextup_poly)
            updated_poly = fix_polygon(updated_poly)
            subbasins_gdf.at[gages[i], 'geometry'] = updated_poly
            area = calc_area(updated_poly)
            subbasins_gdf.at[gage_id, 'unitarea'] = area

            # Insert the new, clipped river reach polyline into our rivers GeoDataFrame
            node_line = comid_line.intersection(updated_poly)
            myrivers_gdf.at[gages[i], 'geometry'] = node_line
            length = calc_length(node_line)
            myrivers_gdf.at[gages[i], 'lengthkm'] = length

        # After the last nested gage, we need to fix the connection info for
        # any rows that previously had the unit catchment with comid as its 'nextdown'
        rows = subbasins_gdf['nextdown'] == comid
        indices = list(rows.index[rows])
        indices.remove(first_node)
        subbasins_gdf.loc[indices, 'nextdown'] = last_node
    return subbasins_gdf


def make_gages_gdf(input_csv: str) -> gpd.GeoDataFrame:
    """
    Reads user data from a CSV file containing information about the desired watershed outlet points.
    and returns a GeoPandas GeoDataFrame where the geometry field contains XY points in unprojected
    lat, lng (CRS 4326).
    """
    # Check that the CSV file is there
    if not isfile(input_csv):
        raise Exception(f"Could not find your outlets file at: {input_csv}")

    if VERBOSE: print(f"Reading your outlets data in: {input_csv}")
    gages_df = pd.read_csv(input_csv, header=0, skipinitialspace=True,
                           dtype={'id': 'str', 'lat': 'float', 'lng': 'float'})
    # Check that the CSV file includes at a minimum: id, lat, lng and that all values are appropriate
    validate(gages_df)
    # Convert gages_df to a GeoPandas GeoDataFrame (adds geography, lets us do geo. operations)
    coordinates = [Point(xy) for xy in zip(gages_df['lng'], gages_df['lat'])]
    gages_gdf = gpd.GeoDataFrame(gages_df, crs=PROJ_WGS84, geometry=coordinates)
    return gages_gdf


def write_outputs(G, myrivers_gdf, subbasins_gdf, gages_list, output_prefix):

    # (0) Save the river network GRAPH
    save_network(G, output_prefix, 'pkl')

    # Save the GEODATA for (1) subbasins, (2) outlets, and (3) rivers

    # (1) SUBBASINS
    fname = f"{OUTPUT_DIR}/{output_prefix}_subbasins.{OUTPUT_EXT}"
    write_geodata(subbasins_gdf, fname)

    # (2) Get the OUTLETS data and write it to disk
    outlets_gdf = subbasins_gdf.copy()
    outlets_gdf['geometry'] = outlets_gdf.apply(lambda row: Point(row['lng'], row['lat']), axis=1)
    fname = f"{OUTPUT_DIR}/{output_prefix}_outlets.{OUTPUT_EXT}"
    write_geodata(outlets_gdf, fname)

    # (3) Write the RIVERS data to disk.
    fname = f"{OUTPUT_DIR}/{output_prefix}_rivers.{OUTPUT_EXT}"
    write_geodata(myrivers_gdf, fname)

    # Finally, return the larger watersheds for each outlet point, if desired
    if WATERSHEDS:
        create_watersheds(G, gages_list, subbasins_gdf)


def create_watersheds(G, gages_list, subbasins_gdf):
    """
    Creates the watersheds (single polygon) upstream of each user outlet point
    and writes geodata to disk as watershed_##.gpkg

    """
    subbasins_gdf.set_index('comid', inplace=True)
    if VERBOSE: print(f"Creating MERGED watersheds geodata for {len(gages_list)} outlets.")
    # Iterate over each of the user's nodes
    for node in gages_list:
        if VERBOSE: print(f"Creating watershed for outlet: {node}")
        upnodes = upstream_nodes(G, node)
        upnodes.append(node)
        mysubs_gdf = subbasins_gdf.loc[upnodes]
        if len(upnodes) == 1:
            watershed_gdf = mysubs_gdf
        else:
            mysubs_gdf.geometry = mysubs_gdf.geometry.apply(lambda p: buffer(p))
            watershed_gs = dissolve_geopandas(mysubs_gdf)
            watershed_gdf = gpd.GeoDataFrame(watershed_gs, columns=['geometry'])
            watershed_gdf['id'] = node

        if FILL:
            watershed_gdf.geometry = watershed_gdf.geometry.apply(lambda p: close_holes(p, FILL_AREA_MAX))

        fname = f"{OUTPUT_DIR}/wshed_{node}.{OUTPUT_EXT}"
        write_geodata(watershed_gdf, fname)


def _run_from_terminal():
    """
    Routine which is run when you call subbasins.py from the command line.
    should be run like:
    python subbasins.py outlets.csv run_name
    """
    # Create the parser for command line inputs.
    description = "Delineate subbasins using data in and input CSV file. Writes " \
        "a set of output files beginning with the output prefix string."

    parser = argparse.ArgumentParser(description=description)

    # Add the arguments
    parser.add_argument('input_csv', help="Input CSV filename, for example 'gages.csv'")
    parser.add_argument('output_prefix', help="Output prefix, a string. The output files will start with this string")

    # Parse the arguments
    args = parser.parse_args()

    # Call the main function, passing the command line arguments
    delineate(args.input_csv, args.output_prefix)


def main():
    # Run directly, for convenience or during development and debugging
    input_csv = 'test_inputs/susquehanna.csv'
    out_prefix = 'susquehanna'
    delineate(input_csv, out_prefix)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Run with command-line arguments
        _run_from_terminal()
    else:
        main()
