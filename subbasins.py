r"""
Delineation of waterhshed subbasins, using data from
MERIT-Basins and MERIT-Hydro.
Created by Matthew Heberger, May 2024.

See README for detailed instructions.

Quick guide:

First, set parameters in the file config.py.

Run it from the command line like this.
$ python subbsins.py outlets.csv testrun

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
import numpy as np
from shapely.geometry import Point
import topojson

# My stuff
from subbasins_config import *  # This file contains a bunch of variables
from py.util import *  # This one contains a bunch of functions
from py.graph_tools import *  # Functions for working with river network information as a Python NetworkX graph
from py.merit_detailed import split_catchment
from py.plot_network import draw_graph

# Shapely throws a bunch of FutureWarnings. Safe to ignore for now, as long as we
# are using a virtual environment, and use the library versions in requirements.txt.
warnings.simplefilter(action='ignore', category=FutureWarning)


def delineate(input_csv: str, output_prefix: str):
    """
    Finds the watershed for a set of outlets.
    Make sure to set the variables in `subbasins_config.py` before running.

    Version #2, where all the points are in a single watershed,
    and we are interested in returning a set of subbasins.

    Reads a list of outlet points from a .csv file.
    The first row should contain the main watershed outlet.
    Subsequent rows of the .csv should contain uptstream points
    within the watershed that we will use to create *subbasin* outlets.

    Outputs geodata (.shp, .gpkg, etc.) and optionally a CSV file with a summary of results
    (including the watershed id, names, and areas).

    """

    def addnode(B: list, node):
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
    if not os.path.isfile(input_csv):
        raise Exception(f"Could not find your outlets file at: {input_csv}")

    # Read the outlet points CSV file and put the data into a Pandas DataFrame
    # (I call the outlet points gages, because I usually in delineated watersheds at streamflow gages)
    if VERBOSE: print(f"Reading your outlets data in: {input_csv}")
    gages_df = pd.read_csv(input_csv, header=0, dtype={'id': 'str', 'lat': 'float', 'lng': 'float'})

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
    gages_df['area_calc'] = 0
    gages_df['result'] = "failed"
    num_gages = len(gages_df)

    # Dict failed_dict: key = id, value = string, explanation of failure
    failed_dict = {}

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
    gages_joined_gdf = gpd.overlay(points_gdf, catchments_gdf, how="intersection")

    # For any gages for which we could not find a unit catchment, add them to Dict failed_dict
    gages_matched_list = gages_joined_gdf['id'].tolist()
    gage_basin_ids = gages_df['id'].tolist()
    for wid in gage_basin_ids:
        if wid not in gages_matched_list:
            failed_dict[wid] = f"Could not assign to a unit catchment in Level 2 basin #{megabasin}"

    # First, let us find the set of unit catchments upstream of the outlet.
    # Let wid be the watershed ID. Get the lat, lng coords of the gage.
    terminal_node_id = gages_joined_gdf['id'].iloc[0]

    # The terminal comid is the unit catchment that contains (overlaps) the outlet point
    terminal_comid = gages_joined_gdf['COMID'].iloc[0]

    # Let B be the list of unit catchments (and river reaches) that are in the basin
    B = []

    # Add the first node, and the rest will be added recursively
    addnode(B, terminal_comid)

    # Next, check that all the other points provided by the user are
    # in this set; otherwise, it means they are not upstream.
    n = len(gages_joined_gdf)
    problemFound = False
    for i in range(0, n):
        id = gages_joined_gdf.at[i, 'id']
        comid = gages_joined_gdf.at[i, 'COMID']
        if comid not in B:
            print(f"Problem: Point with id = {id} is in the watershed of the first point.")
            problemFound = True

    if problemFound:
        raise Exception("One or more points not in watershed. Please check your inputs and try again.")

    # Next, check that we are not trying to subdivide a single unit catchment more than once.
    comids = gages_joined_gdf['COMID'].tolist()

    if not has_unique_elements(comids):
        print("Problem: More than one point falls in each unit catchment. "
              "This script is not currently set up to handle this")

        # Get a list of the offending points
        repeated_comids = find_repeated_elements(comids)
        for comid in repeated_comids:
            ids_rows = gages_joined_gdf[gages_joined_gdf['comid'] == comid, 'id']
            ids_list = ids_rows.tolist()
            print(f"The following gages are in unit catchment {comid}: {', '.join(ids_list)}")

    # Next, we need to split every subcatchment that contains an outlet point.
    # First, we will create a DataFrame with the subcatchment data.
    catchments_gdf.set_index('COMID', inplace=True)
    subbasins_gdf = catchments_gdf.loc[B]

    if PLOTS: plot_basins(subbasins_gdf, points_gdf, 'before')

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
    if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{output_prefix}_network_before')

    # Create a dictionary of the new nodes to add : the comid that we're inserting them into!
    gages_joined_gdf.set_index('id', inplace=True)
    new_nodes = gages_joined_gdf['COMID'].to_dict()

    # Insert the new nodes into the flow network.
    for node, comid in new_nodes.items():
        G = insert_node(G, node, comid)

    if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{output_prefix}_network_after')

    # Now we have a new, accurate network topology. Now we only need to take care of the geography,
    # by splitting the polygons to find the geography (catchment polygon) for the
    # the newly inserted nodes in the river network.

    # First, we can insert the new nodes into our subbasins GeoDataFrame,
    # and update the topology data (column `nextdown`) based on the Graph.
    G = calculate_strahler_stream_order(G)
    G = calculate_shreve_stream_order(G)

    # First insert the new nodes
    for node, comid in new_nodes.items():
        lat = gages_joined_gdf.at[node, 'lat']
        lng = gages_joined_gdf.at[node, 'lng']
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
        lat = gages_joined_gdf.at[node, 'lat']
        lng = gages_joined_gdf.at[node, 'lng']
        catchment_poly = subbasins_gdf.loc[comid, 'geometry']
        node_type = subbasins_gdf.loc[node, 'type']
        isLeaf = (node_type == 'leaf')

        # Split catchment routine
        new_poly, lat_snap, lng_snap = split_catchment(node, megabasin, lat, lng, catchment_poly, isLeaf)

        # Assign the split polygon as the geometry of the new node.
        #  This is the subcatchment polygon upstream of the user-defined outlet.
        subbasins_gdf.at[node, 'geometry'] = new_poly

        # Calculate the area of the new unit catchment and insert value into GeoDataFrame
        area = calc_area(new_poly)
        subbasins_gdf.at[node, 'unitarea'] = round(area, 1)

        # Now SUBTRACT the new polygon from the target unit catchments polygon
        # and update the geometry of the unit catchment `comid`, as we have removed some of its area.
        # shapely's `difference()` method works like a cookie cutter
        target_polygon = subbasins_gdf.at[comid, 'geometry']
        revised_poly = target_polygon.difference(new_poly)
        cleaned_poly = fix_polygon(revised_poly)
        subbasins_gdf.at[comid, 'geometry'] = cleaned_poly

        # Update the area of the clipped unit catchment and insert value into GeoDataFrame
        area = calc_area(cleaned_poly)
        subbasins_gdf.at[comid, 'unitarea'] = round(area, 1)

    if PLOTS: plot_basins(subbasins_gdf, points_gdf, 'after')

    # After splitting the unit catchments, we can remove the most downstream node
    terminal_comid = new_nodes[terminal_node_id]
    subbasins_gdf = subbasins_gdf.drop(terminal_comid)

    # We should also remove it from our Graph representation
    G.remove_node(terminal_comid)

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

    # Now in both the new nodes and the unit catchment, CLIP the river using the polygon boundary.
    for node, comid in new_nodes.items():
        for id in [node, comid]:
            if id != terminal_comid:
                river_polyline = myrivers_gdf.at[id, 'geometry']
                basin_polygon = subbasins_gdf.at[id, 'geometry']
                clipped_line = river_polyline.intersection(basin_polygon)
                myrivers_gdf.at[id, 'geometry'] = clipped_line
                lengthkm = calc_length(clipped_line)
                myrivers_gdf.at[id, 'lengthkm'] = round(lengthkm, 1)

    # Now we can drop the terminal unit catchment from the rivers dataframe.
    myrivers_gdf = myrivers_gdf.drop(terminal_comid)

    # CHECK FOR Null Geometries.
    # If the user has placed one of their points very close to an existing basin outlet,
    # we can get a null geometry. This then creates problems. So remove it.

    null_nodes = subbasins_gdf.index[subbasins_gdf['unitarea'] == 0].tolist()
    subbasins_gdf.drop(null_nodes, inplace=True)

    # Remove this node from the graph
    for node in null_nodes:
        prune_node(G, node)

    if len(null_nodes) > 0:
        if VERBOSE: print(f"Pruned {len(null_nodes)} empty unit catchments from the river network.")
        if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{output_prefix}_network_pruned')

    null_rivers = myrivers_gdf.index[myrivers_gdf['lengthkm'] == 0].tolist()
    myrivers_gdf.drop(null_rivers, inplace=True)

    # SIMPLIFY the geodata?
    if SIMPLIFY:
        topo = topojson.Topology(subbasins_gdf, prequantize=False)
        subbasins_gdf = topo.toposimplify(0.0008).to_gdf()
        topo = topojson.Topology(myrivers_gdf, prequantize=False)
        myrivers_gdf = topo.toposimplify(0.0008).to_gdf()

    # SAVE NETWORK Graph.
    # Before saving, add the unitarea as an attribute in the graph.
    for node in list(G.nodes):
        area = subbasins_gdf.at[node, 'unitarea']
        G.nodes[node]['area'] = round(area, 1)

    save_network(G, output_prefix, 'json')

    # EXPORT the geodata for (1) subbasins, (2) outlets, and (3) rivers
    # (1) SUBBASINS
    subbasins_gdf.reset_index(inplace=True)
    subbasins_gdf.rename(columns={'index': 'comid'}, inplace=True)
    fname = f"{OUTPUT_DIR}/{output_prefix}_subbasins.{OUTPUT_EXT}"
    write_geodata(subbasins_gdf, fname)

    # (2) Get the OUTLETS data and write it to disk
    outlets_gdf = subbasins_gdf.copy()
    outlets_gdf['geometry'] = outlets_gdf.apply(lambda row: Point(row['lng'], row['lat']), axis=1)
    fname = f"{OUTPUT_DIR}/{output_prefix}_outlets.{OUTPUT_EXT}"
    write_geodata(outlets_gdf, fname)

    # (3) Write the RIVERS data to disk.
    myrivers_gdf.reset_index(inplace=True)
    myrivers_gdf.rename(columns={'index': 'comid'}, inplace=True)
    fname = f"{OUTPUT_DIR}/{output_prefix}_rivers.{OUTPUT_EXT}"
    write_geodata(myrivers_gdf, fname)

    if VERBOSE: print("Ran succesfully!")


def _run():
    # Create the parser for command line inputs.
    description = "Delineate subbasins using data in and input CSV file. Writes " \
                   "a set of output files beginning with the output prefix string."

    parser = argparse.ArgumentParser(description=description)

    # Add the arguments
    parser.add_argument('input_csv', help="Input CSV filename, for example 'gages.csv'")
    parser.add_argument('output_ext', help="Output prefix, a string. The output files will start with this string")

    # Parse the arguments
    args = parser.parse_args()

    # Call the main function, passing the command line arguments
    delineate(args.input_csv, args.output_ext)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Run with command-line arguments
        _run()
    else:
        # Run directly, for convenience or during development and debugging
        input_csv = 'test_inputs/outlets4.csv'
        output_prefix = 'Ice4'
        delineate(input_csv, output_prefix)
