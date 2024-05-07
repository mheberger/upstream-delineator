"""
Step #4: Alternative approach.

This is meant to deal with some of the small-ish basins that were not merged in the previous step
because they were not order 1. Instead, they are internal basins that are in series.

So the rule will be

IF
    node has
    shreve > 1
    AND numup = 1
    AND area < threshold

THEN
    merge the node with its UPSTREAM neighbor.

This means
    delete the upstream node
    merge the upstream node's geometry with the downstream node's geometry
    sum the areas of the nodes that are merged

    merge the upstream node's river reach polyline with the node's polyline.


"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import pandas as pd
from sqlalchemy import create_engine
import networkx as nx

# Megabasin to process
id = 11

# Area in km² of a unit catchment below which we will eliminate it by merging it with its upstream neighbor
threshold = 100

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

# Read in data from Step #3
sql = f"""
SELECT comid as comid, unitarea, nextdown, shreve, numup
FROM merit_basins2_{id}
"""

# Read the result of the SQL query into an ordinary Pandas DataFrame
df = pd.read_sql(sql, conn)
df.set_index('comid', inplace=True)

# Create a networkx object out of the rivers data.
G = nx.DiGraph()

# Populate the graph's nodes and edges.
for node_id, nextdown in df['nextdown'].items():
    # Add node with comid as node ID
    G.add_node(node_id)
    # Add edge from comid to nextdown
    if nextdown > 0:
        G.add_edge(node_id, nextdown)

# Find the set of SERIAL NODES that are *candidates for removal*
sql = f"""
SELECT 
    comid, 
    shreve
FROM 
    merit_basins2_{id}
WHERE 
    shreve > 1
    AND numup = 1
    AND unitarea < {threshold}
"""

nodes_df = pd.read_sql(sql, conn)

nodes_list = nodes_df['comid'].values.tolist()
nodes_df.set_index('comid', inplace=True)

# Let's use the graph as an easy way to get the upstream nodes. Add this info to our nodes DataFrame
nodes_df['upnode'] = 0
nodes_df = nodes_df.astype({'upnode': int})

# Add upstream node to our candidate nodes DataFrame
for node in nodes_list:
    upnodes = list(G.predecessors(node))
    if len(upnodes) > 1:
        print("Something wrong")
    else:
        nodes_df.loc[node, 'upnode'] = upnodes[0]

# Modify the DataFrame to add information about where merging should occur.
df['mergeid'] = df.index

# Here is the algorithm for merging the serial nodes
# The rule is, merge the upstream unit catchment IFF their combined area is below our threshold.
# If we have done a merge, and there is another upstream serial node, use the same rule to consider adding it.
# Continue thusly until either the accumulated area > threshold or there are no more upstream nodes.

# First, sort nodes_df by shreve order , Z > A, so we will always be starting downstream and working up.
nodes_df.sort_values(by='shreve', ascending=False, inplace=True)

# Now, we are going to iterate over the nodes, considering what to do with them one at a time.
# Re-extract the nodes list to get it in the order of Shreve order, high to low.
nodes_list = nodes_df.index.tolist()

merges = {}
accum_area = 0
node = nodes_list.pop(0)
node_area = df.loc[node, 'unitarea']
upnode = nodes_df.loc[node, 'upnode']

while len(nodes_list) > 0:
    upnode_area = df.loc[upnode, 'unitarea']
    # If the accumulated area of the new node will be less than threshold, add it
    if accum_area + upnode_area < threshold:
        accum_area = accum_area + upnode_area  # update the accumulated area

        # Add the upstream node to the merges Dictionary
        if node not in merges:
            merges[node] = [upnode]
        else:
            merges[node].append(upnode)

        # Check to see if we can continue going upstream
        if upnode in nodes_list:
            nodes_list.remove(upnode)
            upnode = nodes_df.loc[upnode, 'upnode']
        else:
            node = nodes_list.pop(0)
            node_area = df.loc[node, 'unitarea']
            accum_area = node_area
            upnode = nodes_df.loc[node, 'upnode']
    else:
        node = nodes_list.pop(0)
        node_area = df.loc[node, 'unitarea']
        accum_area = node_area
        upnode = nodes_df.loc[node, 'upnode']

# Do a quick statistical summary of how many merges we have done and how many subbbasins are merged each time
stats = {}
for k, v in merges.items():
    num_merged = len(v)
    if num_merged not in stats:
        stats[num_merged] = 1
    else:
        stats[num_merged] = stats[num_merged] + 1

# Now, put the information about the merging into the DataFrame.
for target, sources in merges.items():
    for source in sources:
        df.loc[source, 'mergeid'] = target
        df.loc[source, 'nextdown'] = df.loc[target, 'nextdown']

# TODO: Here we are missing a critical step!
for target, sources in merges.items():
    source = sources.pop()
    rows = df[df['nextdown'] == source]
    indices = rows.index.tolist()
    for index in indices:
        df.loc[index, 'nextdown'] = target

# Next, write the DataFrame to a temporary table. Then we can do a table join, and do a dissolve, to create step 4 output.
df.reset_index(inplace=True)
df = df[['comid', 'mergeid', 'nextdown']]

df.to_sql('a_temp', engine, if_exists='replace', index=False)

sql = f"""
DROP TABLE IF EXISTS merit_basins4_{id}; 

CREATE TABLE merit_basins4_{id} AS
SELECT 
    mergeid as comid, 
    ST_UNION(geom) as geom, 
    SUM(M.unitarea) as unitarea, 
    MAX(a_temp.nextdown) as nextdown
FROM merit_basins2_{id} as M
LEFT JOIN a_temp ON M.comid = a_temp.comid
GROUP by mergeid;
"""

cur = conn.cursor()
cur.execute(sql)
conn.commit()

# Now let's make the Step #4 RIVERS

sql = f"""
DROP TABLE IF EXISTS merit_rivers4_{id};

CREATE TABLE merit_rivers4_{id} AS
SELECT 
    mergeid as comid, 
    ST_LineMerge(ST_UNION(geom_simple)) as geom, 
    SUM(R.lengthkm) as lengthkm
FROM 
    merit_rivers2_{id} as R
LEFT JOIN 
    a_temp ON R.comid = a_temp.comid
GROUP by 
    mergeid;
"""
cur.execute(sql)

# Close the database connection
conn.commit()
cur.close()
conn.close()
conn.close()
engine.dispose()

