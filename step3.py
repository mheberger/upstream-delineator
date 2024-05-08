"""
Step #3, Eliminate small "junction nodes."

The decision rule here is that the unit catchment area is small, and river reach length is short.

If that is the case, we will eliminate the node, and reconnect the graph edges as appropriate.

Afterwards, we will also make an effort to update the rivers centerline vertices, just for appearances.

The criteria for a "minor junction" that can be collapsed is:
1. Cannot be a terminal coastal catchment: nextdown !=0
2. The distance between the ends of the polyline are less than a threshold.
  (We can use the `lenghtkm` attribute -- it is close enough.)
3. It should have one or more upstream connections: numup > 0

TODO: Search for junctions that are in series, as these can cause errors. Need to look for 2, 3, 4 in a row.
  Just go recursive.

TODO: Often needs to be run multiple times because on the first pass, it is not able to handle the situation where
  there are multiple little junctions in series.

"""
import warnings

from update_network_data import update_network

warnings.filterwarnings('ignore')
import pandas as pd
from util.db import connection, cursor, engine, db_close
from util.graph_tools import make_river_network


def step3(area_threshold: float, length_threshold: float, basins_tbl: str, rivers_tbl: str):
    """
    Step #3 removes small "junction" nodes from the network.

    Inputs:
        area_threshold: Area of the unit catchment, in km²
        length_threshold: Length of the river reach, in km
        basins_tbl: name of the table to create to store catchments geodata
        rivers_tbl: name of the table to create to store rivers geodata
    """

    # First, let's read the river network into a GRAPH
    # Define your SQL query
    sql = f"""
    SELECT comid, nextdown
    FROM basins2
    """

    # Read the result of the SQL query into a GeoDataFrame
    df = pd.read_sql(sql, connection)
    df.set_index('comid', inplace=True)

    # Create a networkx object out of the rivers data.
    G = make_river_network(df)

    # Following query finds the candidates for junction nodes (small catchment area AND short river length)
    sql = f"""
    SELECT 
        R.comid
    FROM 
        basins2 as M
    INNER JOIN 
        rivers2 R ON M.comid = R.comid
    WHERE 
        R.lengthkm < {length_threshold}
        AND
        M.unitarea < {area_threshold}
        AND
        M.nextdown > 0;
    """

    cursor.execute(sql)
    rows = cursor.fetchall()
    nodes = [row[0] for row in rows]

    print(f"Found {len(nodes)} tiny junctions to collapse.")

    # Need to make sure that none of these tiny junctions are in series with one another,
    #  or this will create problems!!! The solution may be to go back to my old way of doing things
    #  and do more of a graph pruning. When I'm working with the tabular data like this, it is too hard
    #  to keep track of updated relationships. However, for now, let's opt for an easy solution.

    for node in nodes:
        # If the downstream node is also in our list, remove it from this list!
        ds_node = df.loc[node, 'nextdown']
        if ds_node in nodes:
            nodes.remove(ds_node)
            print(ds_node)

    print(f"After checking for nodes in series, we have {len(nodes)}")

    # Modify the DataFrame to add information about where merging should occur.
    df['mergeid'] = df.index

    # Now, remove these nodes from the graph.
    # As we do so, we need to store information about the new connectivity and updated areas
    new_connections = {}  # This will store any updated from: to connections, so we can fix the river centerline geometries

    for node in nodes:
        # Find nextdown
        nextdown = df.loc[node, 'nextdown']
        # We are going to merge its geometry with its downstream neighbor
        df.loc[node, 'mergeid'] = nextdown

        predecessors = list(G.predecessors(node))
        for predecessor in predecessors:
            new_connections[predecessor] = nextdown

    # We will let the SQL query take care of summing the unitareas.
    # Next, write the DataFrame to a temporary table. Then we can do a table join, and do a dissolve, to create step 5 output.
    df.reset_index(inplace=True)

    df.to_sql('a_temp', engine, if_exists='replace', index=False)

    # Next, we'll join on this table and dissolve the geometries.
    sql = f"""
    DROP TABLE IF EXISTS {basins_tbl}; 
    
    CREATE TABLE {basins_tbl} AS
    SELECT 
        mergeid as comid, 
        ST_UNION(geom) as geom, 
        SUM(M.unitarea) as unitarea
    FROM basins2 as M
    LEFT JOIN a_temp ON M.comid = a_temp.comid
    GROUP by mergeid;
    """

    cursor.execute(sql)

    # Could not figure out a good way to update the field `nextdown`, so let's update them
    # from the original table and then modify just the ones that have changed.
    sql = f"""
    ALTER TABLE {basins_tbl}
    ADD COLUMN IF NOT EXISTS nextdown INTEGER;
    
    UPDATE
        {basins_tbl}
    SET 
        nextdown = basins2.nextdown
    FROM 
        basins2
    WHERE 
        {basins_tbl}.comid = basins2.comid
    """
    cursor.execute(sql)

    # Now, go in and fix individual rows where we made new, updated connections;
    for from_node, to_node in new_connections.items():
        sql = f"""
        UPDATE {basins_tbl} 
        SET nextdown =  {to_node} 
        WHERE comid = {from_node};"""
        cursor.execute(sql)

    # NEXT, update the rivers geometries.
    # First, copy the rivers table from the previous step.
    # But omit the rivers for the nodes that we eliminated.
    nodes = [str(x) for x in nodes]
    nodes_list_str = ','.join(nodes)
    sql = f"""
    DROP TABLE IF EXISTS 
        {rivers_tbl};
    CREATE TABLE {rivers_tbl} AS
    SELECT *
    FROM rivers2
    WHERE comid NOT IN ({nodes_list_str});
    """
    cursor.execute(sql)

    # Next, update the geometries of the incoming nodes.
    # HACK: I am just moving the end point to the new outlet location!!!
    # TODO: Would probably be better to merge geometries (but then we will have overlapping rivers, confusing!)

    # For simplicity, we can just do them one at a time.

    for from_node, to_node in new_connections.items():
        # print(f"{from_node}: {to_node}")
        # First, get the coordinates of the *BEGINNING* of the *downstream* river centerline
        # This is our target
        sql = f"""
        SELECT 
            ST_X(ST_EndPoint(geom)) AS end_x,
            ST_Y(ST_EndPoint(geom)) AS end_y
        FROM {rivers_tbl}
        WHERE comid = {to_node};
        """
        cursor.execute(sql)
        row = cursor.fetchone()
        x = row[0]
        y = row[1]

        # Now modify the incoming river so it connects to this point!
        sql = f"""
        UPDATE {rivers_tbl}
        SET geom = ST_SetPoint(
                            geom,
                            0,
                            ST_MakePoint({x}, {y})
                         )
        WHERE comid = {from_node};
        """
        cursor.execute(sql)

    db_close()


if __name__ == "__main__":
    step3(2, 1, "basins_s_11", "rivers_s_11")
    update_network("basins_s_11")
