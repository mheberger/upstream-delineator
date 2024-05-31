r"""
Delineation of watershed subbasins, using data from
MERIT-Basins and MERIT-Hydro.
Created by Matthew Heberger, May 2024.

See README for detailed instructions.

Quick guide:

First, set parameters in the file subbasins_config.py.

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
from shapely.geometry import Point
import topojson

# My stuff
from subbasins_config import *  # This file contains a bunch of variables to be set before running this script
from py.consolidate import consolidate_network, show_area_stats
from py.util import *  # Contains a bunch of functions
from py.graph_tools import *  # Functions for working with river network information as a Python NetworkX graph
from py.merit_detailed import split_catchment
from py.plot_network import draw_graph
from py.fast_dissolve import dissolve_geopandas, buffer, close_holes

# Shapely throws a bunch of FutureWarnings. Safe to ignore for now, as long as we
# are using a virtual environment, and use the library versions in requirements.txt.
warnings.simplefilter(action='ignore', category=FutureWarning)


def delineate(INPUT_CSV: str, OUTPUT_PREFIX: str):
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
    (including the watershed gid, names, and areas).

    """

    def addnode(B: list, node_id):
        """"
        Recursive function to assemble the list of upstream unit catchments
        upstream_comids is a Python List of the unit catchments that make up our watershed.
        upstream_comids is for BASIN...
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

    # Make sure the user folders from `config.py` exist (for output & plots). If not create them.
    make_folders()

    # Read the outlet points CSV file and put the data into a Pandas DataFrame
    # (I call the outlet points gages, because I usually in delineated watersheds at streamflow gages)
    gages_gdf = make_gages_gdf(INPUT_CSV)
    num_gages = len(gages_gdf)

    # Get the megabasin(s) in which the points are located (they need to all be in the same megabasin).
    megabasin = get_megabasin(gages_gdf)

    if VERBOSE:
        print('Reading data table for unit catchments in basin %s' % megabasin)
    catchments_gdf = load_gdf("catchments", megabasin, True)

    # The _network_ data is in the RIVERS file rather than the CATCHMENTS file
    # (this is just how the MERIT-Basins authors did it)
    if VERBOSE: print('Reading data table for rivers in basin %s' % megabasin)
    rivers_gdf = load_gdf("rivers", megabasin, True)
    rivers_gdf.set_index('COMID', inplace=True)
    # We wish to report the outlet point for each subbasin.
    # We can get this information from the river polylines.
    rivers_gdf['end_point'] = rivers_gdf['geometry'].apply(lambda x: x.coords[0])
    rivers_gdf['lng'] = rivers_gdf['end_point'].apply(lambda x: x[0])
    rivers_gdf['lat'] = rivers_gdf['end_point'].apply(lambda x: x[1])

    # Perform a Spatial join on gages (points) and the unit catchments (polygons)
    # to find the corresponding unit catchment for each gage
    # Adds the fields COMID and unitarea
    if VERBOSE: print(f"Performing spatial join on {num_gages} outlet points in basin #{megabasin}")
    gages_list = gages_gdf['gid'].tolist()  # Get the list before doing the join.
    gages_gdf = gpd.overlay(gages_gdf, catchments_gdf, how="intersection")
    gages_gdf.set_index('gid', inplace=True)

    # For any gages for which we could not find a unit catchment, add issue a warning
    # Basically checking which rows do not appear after doing the overlay
    gages_matched_list = gages_gdf.index.tolist()
    for gid in gages_list:
        if gid not in gages_matched_list:
            raise Warning(f"Could not assign to a unit catchment to gage with gid {gid}")

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
    for gid in gage_list:
        comid = gages_gdf.at[gid, 'COMID']
        if comid not in upstream_comids:
            gages_gdf.drop(gid, inplace=True)
            raise Warning(f"The point with gid = {gid} is not contained in the watershed of the first point.")

    # subbasins_gdf is the set of unit catchments in our watershed. This will ultimately become our output
    catchments_gdf.set_index('COMID', inplace=True)
    subbasins_gdf = catchments_gdf.loc[upstream_comids]
    # Add lat, lng, and NextDownID to subbasins_gdf.
    subbasins_gdf = subbasins_gdf.join(rivers_gdf[['lat', 'lng', 'NextDownID']])
    # Re-name the NextDownID field, and make sure it is an integer
    subbasins_gdf.rename(columns={'NextDownID': 'nextdown'}, inplace=True)
    subbasins_gdf['nextdown'] = subbasins_gdf['nextdown'].astype(int)

    if PLOTS: plot_basins(subbasins_gdf, gages_gdf, 'before')

    # Create two copies of river network data!
    # `allrivers_gdf` will contain all the available polylines in the watershed, a nice
    # visual representation of the river network.
    # `myrivers_gdf` will contain the topologically correct network where there is exactly
    # one polyline per unit catchment.
    allrivers_gdf = rivers_gdf.loc[upstream_comids]
    fname = f"{OUTPUT_DIR}/{OUTPUT_PREFIX}_allrivers.{OUTPUT_EXT}"
    write_geodata(allrivers_gdf, fname)
    del allrivers_gdf

    # With this version, we will try to create a topologically correct visual representation
    # of the river network, where there is a 1:1 mapping of river reaches to subbasins.
    myrivers_gdf = rivers_gdf.loc[upstream_comids]

    # Iterate over the gages, and run SPLIT CATCHMENT for every gage
    # For now, we are going to put the split polygon geometry into a field in gages_gdf
    gages_gdf['polygon'] = None
    gages_gdf['polygon_area'] = 0

    for gage_id in gages_gdf.index:
        comid = gages_gdf.at[gage_id, 'COMID']
        lat = gages_gdf.at[gage_id, 'lat']
        lng = gages_gdf.at[gage_id, 'lng']
        catchment_poly = subbasins_gdf.loc[comid, 'geometry']
        uparea = rivers_gdf.at[comid, 'uparea']
        # TODO: It would be better to extract this information from the network graph
        isLeaf = (uparea < 50)  # A leaf is an unit catchment with no upstream neihbor.
        # SPLIT
        new_poly, lat_snap, lng_snap = split_catchment(gage_id, megabasin, lat, lng, catchment_poly, isLeaf)
        gages_gdf.at[gage_id, 'polygon'] = new_poly
        gages_gdf.at[gage_id, 'lat_snap'] = lat_snap
        gages_gdf.at[gage_id, 'lng_snap'] = lng_snap

        # Find the area of the clipped unit catchment and insert value into GeoDataFrame
        area = calc_area(new_poly)
        gages_gdf.at[gage_id, 'polygon_area'] = round(area, 1)

    # Now, handle the case where we had more than one gage point in a single unit catchment.
    #   (1) figure out their order from upstream to downstream based on split catchment area
    #   (2) clip the polygons as appropriate
    #   (3) insert this data into subbasins_gdf, including polygon geometry and `nextdown` gage id
    comids = gages_gdf['COMID'].tolist()
    repeats = find_repeated_elements(comids)
    singles = [item for item in comids if item not in repeats]

    gages_gdf.drop(columns=['lat', 'lng', 'geometry', 'unitarea'], inplace=True)
    rnmap = {'lat_snap': 'lat',
             'lng_snap': 'lng',
             'polygon_area': 'unitarea',
             'polygon': 'geometry'
             }
    gages_gdf.rename(columns=rnmap, inplace=True)

    # First, handle all of the gages where there is only one gage in a unit catchment (standard treatment)
    # These will always be upstream of the unit catchment that we are inserting it into
    # (a) insert a new row into subbasins
    if len(singles) > 0:
        selected_rows = gages_gdf[gages_gdf['COMID'].isin(singles)]
        new_nodes = selected_rows['COMID'].to_dict()
        selected_rows.drop(columns=['COMID'], inplace=True)
        selected_rows['nextdown'] = -999

        # TODO This line seems to be creating problems with types.
        subbasins_gdf = pd.concat([subbasins_gdf, selected_rows])

        # Find and fix the 'nextdown' attribute for any nodes that were u/s of terminal catchments
        # into which we just inserted our nodes.
        for node, comid in new_nodes.items():
            rows = subbasins_gdf['nextdown'] == comid
            subbasins_gdf.loc[rows, 'nextdown'] = node
            subbasins_gdf.loc[node, 'nextdown'] = comid

            # Subtract the geometry of the new node from its parent
            if True:
                current_polygon = subbasins_gdf.loc[comid, 'geometry']
                new_poly = subbasins_gdf.loc[node, 'geometry']
                old_poly = current_polygon.difference(new_poly)
                old_poly = fix_polygon(old_poly)
                subbasins_gdf.loc[comid, 'geometry'] = old_poly

                # Split the river polyline at the new unit catchment boundary
                # First, get the piece in what is left of the old unit catchment
                river_line = rivers_gdf.at[comid, 'geometry']

                old_river = river_line.intersection(old_poly)
                old_river = multilinestring_to_linestring(old_river)
                length = calc_length(old_river)
                myrivers_gdf.at[comid, 'geometry'] = old_river
                myrivers_gdf.at[comid, 'lengthkm'] = length

                new_river = river_line.intersection(new_poly)
                new_river = multilinestring_to_linestring(new_river)
                length = calc_length(new_river)
                myrivers_gdf.at[node, 'geometry'] = new_river
                myrivers_gdf.at[node, 'lengthkm'] = length

    # Next, handle the case where there are multiple gages in a unit catchment (special treatment)
    # Handle these one target unit catchment at a time.
    for comid in repeats:
        # Find all of the gages that fall in this unit catchment.
        gages_set = gages_gdf[gages_gdf['COMID'] == comid]
        # We want to handle these in order from smallest area to largest area
        gages_set.sort_values(by='unitarea', inplace=True)
        gages_set['nextdown'] = gages_set.index.to_series().shift(-1).fillna(comid)
        subbasins_gdf = pd.concat([subbasins_gdf, gages_set])

        # Get the (whole) river polyline in the original, unsplit unit catchment
        river_line = rivers_gdf.loc[comid, 'geometry']

        # Let upstream_polygon be the cumulative upstream area that we need to subtract from our new subbasin
        count = 0
        for gage_id in gages_set.index:
            # Handle the first (smallest) upstream point in the unit catchment
            if count == 0:
                new_polygon = gages_set.at[gage_id, 'geometry']
                upstream_polygon = new_polygon
                new_river = river_line.intersection(new_polygon)

            # Handle any subsequent points in the same unit catchment
            else:
                current_polygon = gages_set.at[gage_id, 'geometry']
                old_poly = current_polygon.difference(upstream_polygon)
                upstream_polygon = current_polygon
                cleaned_poly = fix_polygon(old_poly)
                new_polygon = cleaned_poly
                new_river = river_line.intersection(new_polygon)

            # Insert the new, clipped subbasin polygon into our subbasins GeoDataFrame
            count += 1
            subbasins_gdf.at[gage_id, 'geometry'] = new_polygon
            area = calc_area(new_polygon)
            subbasins_gdf.at[gage_id, 'unitarea'] = area
            # Insert the new, clipped river reach polyline into our rivers GeoDataFrame
            myrivers_gdf.at[gage_id, 'geometry'] = new_river
            length = calc_length(new_river)
            myrivers_gdf.at[gage_id, 'lengthkm'] = length
            
        # After the last "nested gage", we need to fix the connection info for
        # any rows that previously had the unit catchment as its 'nextdown'
        rows = subbasins_gdf['nextdown'] == comid
        subbasins_gdf.loc[rows, 'nextdown'] = gages_set.index.to_list()[0]

        # Finally, we need to clip the geometry of the base unit catchment,
        # unless that unit catchment is downstream of the outlet, in which case it is not to be included
        # in the output
        if comid != terminal_comid:
            current_polygon = subbasins_gdf.at[comid, 'geometry']
            old_poly = current_polygon.difference(upstream_polygon)
            old_poly = fix_polygon(old_poly)
            subbasins_gdf.at[comid, 'geometry'] = old_poly
            area = calc_area(old_poly)
            subbasins_gdf.at[comid, 'unitarea'] = area
            new_river = gpd.overlay(river_line, new_polygon, how='overlay')
            myrivers_gdf.at[comid, 'geometry'] = new_river
            length = calc_length(new_river)
            myrivers_gdf.at[comid, 'lengthkm'] = length

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

    if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{OUTPUT_PREFIX}_network')

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
        if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{OUTPUT_PREFIX}_network_pruned')

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

    # Add attributes for subbasin area and river reach length to the
    # drainage netowork Graph
    for node in list(G.nodes):
        area = subbasins_gdf.at[node, 'unitarea']
        G.nodes[node]['area'] = round(area, 1)
        try:
            length = myrivers_gdf.at[node, 'lengthkm']
        except:
            length = 0
        G.nodes[node]['length'] = round(length, 1)

    # If the user wants larger subbasins, we can expand merge adjacent unit catchments.
    # Has to be done carefully. I call it consolidating the network.
    if CONSOLIDATE:
        G, MERGES, rivers2merge, rivers2delete = consolidate_network(G, threshold_area=MAX_AREA,
                                                                     threshold_length=THRESHOLD_LENGTH)
        if NETWORK_DIAGRAMS: draw_graph(G, f'plots/{OUTPUT_PREFIX}_network_final')
        subbasins_gdf['target'] = subbasins_gdf.index

        # Dissolve unit catchments based on information in MERGES.
        for node, target in MERGES.items():
            subbasins_gdf.at[node, 'target'] = target

        agg = {
            'unitarea': 'sum',
            'lat': 'max',
            'lng': 'max',
        }

        if VERBOSE: print("Dissolving geometries. This can take time!")
        subbasins_gdf = subbasins_gdf.dissolve(by="target", aggfunc=agg)
        subbasins_gdf.geometry = subbasins_gdf.geometry.apply(lambda p: buffer(p))
        subbasins_gdf.geometry = subbasins_gdf.geometry.apply(lambda p: close_holes(p, 0))
        subbasins_gdf.reset_index(inplace=True)
        subbasins_gdf.rename(columns={'target': 'comid'}, inplace=True)
        subbasins_gdf.set_index('comid', inplace=True)

    # Add column `nextdownid` and the orders based on data in the graph
    subbasins_gdf['nextdown'] = 0

    for idx in subbasins_gdf.index:
        try:
            nextdown = list(G.successors(idx))[0]
        except:
            nextdown = 0
        subbasins_gdf.at[idx, 'nextdown'] = nextdown

        subbasins_gdf.at[idx, 'strahler_order'] = G.nodes[idx]['strahler_order']
        subbasins_gdf.at[idx, 'shreve_order'] = G.nodes[idx]['shreve_order']

    # round all the areas to one decimal
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
            myrivers_gdf['geometry'] = myrivers_gdf['geometry'].apply(multilinestring_to_linestring)

    myrivers_gdf['end_point'] = myrivers_gdf['geometry'].apply(lambda x: x.coords[0])
    myrivers_gdf['lng'] = myrivers_gdf['end_point'].apply(lambda x: x[0])
    myrivers_gdf['lat'] = myrivers_gdf['end_point'].apply(lambda x: x[1])
    myrivers_gdf.drop(columns=['end_point'], inplace=True)

    try:
        myrivers_gdf.drop(terminal_comid, inplace=True)
    except:
        pass

    for idx in subbasins_gdf.index:
        if idx in myrivers_gdf.index:
            subbasins_gdf.at[idx, 'lat'] = myrivers_gdf.at[idx, 'lat']
            subbasins_gdf.at[idx, 'lng'] = myrivers_gdf.at[idx, 'lng']

    # Finally, write the results to disk
    write_outputs(G, myrivers_gdf, subbasins_gdf, gages_list, OUTPUT_PREFIX)

    if VERBOSE: print("Ran succesfully!")


def make_gages_gdf(input_csv: str) -> gpd.GeoDataFrame:
    """
    Reads user data from a CSV file containing information about the desired watershed outlet points.
    and returns a GeoPandas GeoDataFrame where the geometry field contains XY points in unprojected
    lat, lng (CRS 4326).
    """
    # Check that the CSV file is there
    if not os.path.isfile(input_csv):
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
    save_network(G, output_prefix, 'pkl')
    # EXPORT the geodata for (1) subbasins, (2) outlets, and (3) rivers
    # (1) SUBBASINS
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

    # Finally, return the larger watersheds for each outlet point, if desired
    if WATERSHEDS:
        if VERBOSE: print(f"Creating MERGED watersheds geodata for {len(gages_list)} outlets.")

        if FILL:
            PIXEL_AREA = 0.000000695  # Constant for the area of a single pixel in MERIT-Hydro, in decimal degrees
            area_max = FILL_THRESHOLD * PIXEL_AREA
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
                watershed_gdf.geometry = watershed_gdf.geometry.apply(lambda p: close_holes(p, area_max))

            fname = f"{OUTPUT_DIR}/wshed_{node}.{OUTPUT_EXT}"
            write_geodata(watershed_gdf, fname)


def _run():
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


def test():
    # Run directly, for convenience or during development and debugging
    input_csv = 'test_inputs/outlet_dworshak.csv'
    out_prefix = 'dworshak'
    delineate(input_csv, out_prefix)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Run with command-line arguments
        _run()
    else:
        test()
