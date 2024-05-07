"""
Step #3 We just did some merging. Now we need to 'invert' the network.
For every node, we have the `nextdown` which is the to_node.
However, for each node, we also need to know the number of upstream nodes, `numup`
and have a list of the upstream nodes. Since there are a maximum of 4 in our dataset,
I found it convenient to put these in fields up1, up2, up3, and up4.

This is important later on for the delineation algorithms.

"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import pandas as pd
import networkx as nx
from sqlalchemy import create_engine
from graph_tools import calculate_strahler_stream_order, calculate_shreve_stream_order, calculate_num_incoming


# New: updated so it can do this for different files.
level = 4
id = 11

host = "localhost"
database = "basins"
user = "postgres"
password = "dbpw"
port = 5432

# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host=host,
    database=database,
    user=user,
    password=password
)

# Create a SQLAlchemy engine
engine = create_engine('postgresql://{}:{}@{}:{}/{}'.format(user, password, host, port, database))

# Define your SQL query
sql = f"""
SELECT comid, nextdown
FROM merit_basins{level}_{id};
"""

# Read the result of the SQL query into a GeoDataFrame
df = pd.read_sql(sql, conn)
df.set_index('comid', inplace=True)

# Create a networkx object out of the rivers data.
G = nx.DiGraph()

for node_id, nextdown in df['nextdown'].items():
    # Add node with comid as node ID
    G.add_node(node_id)
    # Add edge from comid to nextdown
    if nextdown > 0:
        G.add_edge(node_id, nextdown)

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
ALTER TABLE merit_basins{level}_{id}
ADD COLUMN IF NOT EXISTS numup INTEGER,
ADD COLUMN IF NOT EXISTS strahler INTEGER,
ADD COLUMN IF NOT EXISTS shreve INTEGER;
"""

cur = conn.cursor()
cur.execute(sql)

sql = f"""
UPDATE merit_basins{level}_{id} AS m
SET numup = a_temp.numup, 
 strahler = a_temp.strahler,
 shreve = a_temp.shreve
FROM a_temp
WHERE m.comid = a_temp.comid;

"""
cur.execute(sql)

conn.commit()
conn.close()

# Close the engine
engine.dispose()
