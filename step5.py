"""
Step #5, second try: Here, we eliminate the little small junction nodes.

The decision rule here is that the unit catchment area is small, and river reach length is short.

If that is the case, we will eliminate the node, and reconnect the graph edges as appropriate.

Afterwards, we will also make an effort to update the rivers centerline vertices, just for appearances.

TODO: Try running this again, because on the first pass, it is not able to handle the situation where
  there are multiple little junctions in series.

"""
import warnings
warnings.filterwarnings('ignore')
import psycopg2
import pandas as pd
from sqlalchemy import create_engine
import networkx as nx

# Megabasin to process
level = 5
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

"""
The criteria for a "minor junction" that can be collapsed is:
1. Cannot be a terminal coastal catchment: nextdown !=0
2. The distance between the ends of the polyline are less than a threshold. 
  (We can use the `lenghtkm` attribute -- it is close enough.)
3. It should have one or more upstream connections: numup > 0 
Below we write a query to identify these, based on the results so far, from Step #4.
"""

# First, let's read the river network into a GRAPH
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

# Next, find the set of candidate nodes
area_threshold = 2  # square kilometers
length_threshold = 1  # kilometers

sql = f"""
SELECT 
    R.comid
FROM 
    merit_basins{level}_{id} as M
INNER JOIN 
    merit_rivers{level}_{id} R ON M.comid = R.comid
WHERE 
    R.lengthkm < {length_threshold}
    AND
    M.unitarea < {area_threshold}
    AND
    M.nextdown > 0;
"""
cur = conn.cursor()
cur.execute(sql)

rows = cur.fetchall()
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
df.to_csv('C:/Users/mheberger/Desktop/temp.csv')

# Next, we'll join on this table and dissolve the geometries.
sql = f"""
DROP TABLE IF EXISTS merit_basins{level + 1}_{id}; 

CREATE TABLE merit_basins{level + 1}_{id} AS
SELECT 
    mergeid as comid, 
    ST_UNION(geom) as geom, 
    SUM(M.unitarea) as unitarea
FROM merit_basins{level}_{id} as M
LEFT JOIN a_temp ON M.comid = a_temp.comid
GROUP by mergeid;
"""

cur = conn.cursor()
cur.execute(sql)

# Could not figure out a good way to update the field `nextdown`, so let's update them
# from the original table and then modify just the ones that have changed.
sql = f"""
ALTER TABLE merit_basins{level + 1}_{id}
ADD COLUMN IF NOT EXISTS nextdown INTEGER;

UPDATE
merit_basins{level + 1}_{id}
SET nextdown = merit_basins{level}_{id}.nextdown
FROM merit_basins{level}_{id}
WHERE 
merit_basins{level + 1}_{id}.comid = merit_basins{level}_{id}.comid
"""
cur.execute(sql)

# Now, go in and fix individual rows where we made new, updated connections;
for from_node, to_node in new_connections.items():
    sql = f"""
    UPDATE merit_basins{level + 1}_{id} 
    SET nextdown =  {to_node} 
    WHERE comid = {from_node};"""
    cur.execute(sql)

# NEXT, update the rivers geometries.
# First, delete the little small rivers that were associated with the nodes we just removed.
nodes = [str(x) for x in nodes]
nodes_list_str = ','.join(nodes)
sql = f"""
DROP TABLE IF EXISTS merit_rivers{level + 1}_{id};
CREATE TABLE merit_rivers{level + 1}_{id} AS
SELECT *
FROM merit_rivers{level}_{id}
WHERE comid NOT IN ({nodes_list_str});
"""
cur.execute(sql)

# Next, update the geometries of the incoming nodes.
# For simplicity, we can just do them one at a time.

for from_node, to_node in new_connections.items():
    # print(f"{from_node}: {to_node}")
    # First, get the coordinates of the beginning of the downstream river centerline
    sql = f"""
    SELECT 
        ST_X(ST_EndPoint(geom)) AS end_x,
        ST_Y(ST_EndPoint(geom)) AS end_y
    FROM merit_rivers{level + 1}_{id}
    WHERE comid = {to_node};
    """
    cur.execute(sql)
    row = cur.fetchone()
    x = row[0]
    y = row[1]

    sql = f"""
    UPDATE merit_rivers{level + 1}_{id}
    SET geom = ST_SetPoint(
                        geom,
                        0,
                        ST_MakePoint({x}, {y})
                     )
    WHERE comid = {from_node};
    """
    cur.execute(sql)

conn.commit()
conn.close()
