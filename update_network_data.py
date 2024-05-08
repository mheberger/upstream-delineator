"""
Update the network data: stream orders and `numup` for each node (the number of incoming edges,
or immediately upstream unit catchments.)

The Shreve Stream order gives us the TOTAL number of upstream nodes.

This is important later on for the delineation algorithms.

"""

import warnings
warnings.filterwarnings('ignore')
import pandas as pd
from util.graph_tools import calculate_strahler_stream_order, calculate_shreve_stream_order, calculate_num_incoming, make_river_network
from util.db import connection, cursor, engine, db_close


def update_network(tbl: str):
    # Define your SQL query
    sql = f"""
    SELECT comid, nextdown
    FROM {tbl};
    """

    # Read the result of the SQL query into a GeoDataFrame
    df = pd.read_sql(sql, connection)
    df.set_index('comid', inplace=True)

    # Create a networkx object out of the rivers data.
    G = make_river_network(df)
    G = calculate_strahler_stream_order(G)
    G = calculate_shreve_stream_order(G)
    G = calculate_num_incoming(G)

    # Now put the information from the graph into the DataFrame, and then use this to update the Posgres GIS layer
    df['shreve'] = 0
    df['strahler'] = 0
    df['numup'] = 0

    df = df.astype({'shreve': int})
    df = df.astype({'strahler': int})
    df = df.astype({'numup': int})

    for node in G.nodes:
        df.loc[node, 'shreve'] = G.nodes[node]['shreve_order']
        df.loc[node, 'strahler'] = G.nodes[node]['strahler_order']
        df.loc[node, 'numup'] = G.nodes[node]['num_incoming']

    # Now use this data to update the DB table
    df.reset_index(inplace=True)
    df.to_sql('a_temp', engine, if_exists='replace', index=False)

    sql = f"""
    ALTER TABLE {tbl}
    ADD COLUMN IF NOT EXISTS numup INTEGER,
    ADD COLUMN IF NOT EXISTS strahler INTEGER,
    ADD COLUMN IF NOT EXISTS shreve INTEGER;
    """

    cursor.execute(sql)

    sql = f"""
    UPDATE {tbl} AS m
    SET numup = a_temp.numup, 
     strahler = a_temp.strahler,
     shreve = a_temp.shreve
    FROM a_temp
    WHERE m.comid = a_temp.comid;
    """
    cursor.execute(sql)

    db_close()


if __name__ == "__main__":
    update_network("basins1")
