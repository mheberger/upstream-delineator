"""
Step #2

Testing whether I can create HUC-like basin boundaries with simple decision rules in GeoPandas

This time, I will try using the rule: For a node, merge it with its upstream nodes if its combined
area would be less than a threshold. The field `uparea` already contains the total upstream area,
including itself.

"""


import psycopg2
import pandas as pd
import networkx as nx
from sqlalchemy import create_engine


def get_all_upstream_nodes(graph, node):
    upstream_nodes = set()

    def dfs_util(current_node):
        # Mark the current node as visited
        upstream_nodes.add(current_node)

        # Traverse all predecessors (incoming edges) of the current node
        for predecessor in graph.predecessors(current_node):
            # Recursively visit all predecessors
            if predecessor not in upstream_nodes:
                dfs_util(predecessor)

    # Start DFS from the given node
    dfs_util(node)

    # Remove the starting node from the set of upstream nodes
    upstream_nodes.discard(node)

    # Return the list of all upstream nodes
    return list(upstream_nodes)

# Connect to your PostgreSQL database
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
cur = conn.cursor()

# Create a SQLAlchemy engine
engine = create_engine('postgresql://{}:{}@{}:{}/{}'.format(user, password, host, port, database))

id = 27

sql = f"""
SELECT comid, nextdown, uparea
FROM merit_basins_{id}
WHERE nextdown IS NOT Null;
"""
df = pd.read_sql(sql, conn)
df.set_index('comid', inplace=True)

df = df.sort_values(by='uparea', ascending=False)

# Get the sorted indices
node_list = df.index.tolist()

# Create a networkx object out of the rivers data.
G = nx.DiGraph()

for comid, nextdown in df['nextdown'].items():
    #print(f"{comid}: {nextdown}")
    # Add node with comid2 as node ID
    G.add_node(comid)
    G.nodes[comid]['uparea'] = df.loc[comid, 'uparea']
    # Add edge from comid2 to nextdown
    if nextdown > 0:
        G.add_edge(comid, nextdown)

merges = {}

# Do some merging of nodes where uparea < threshold
for node_id in node_list:
    if G.nodes[node_id]['uparea'] < 1000:
        # Remove its predecessors from the graph
        upstream_nodes = get_all_upstream_nodes(G, node_id)
        if len(upstream_nodes) > 1:
            for upnode in upstream_nodes:
                outgoing_edges = list(G.out_edges(upnode))
                # Remove all outgoing edges from the upstream node we will remove
                G.remove_edges_from(outgoing_edges)
                G.remove_node(upnode)
                node_list.remove(upnode)
                df.loc[upnode, 'mergeid'] = node_id

# Now I have a pruned network graph and an updated dataframe. Now do the dissolve
df.loc[df['mergeid'].isnull(), 'mergeid'] = df.index[df['mergeid'].isnull()]
df['mergeid'] = df['mergeid'].astype(int)

df.reset_index(inplace=True)

# Write DataFrame to PostgreSQL database
df.to_sql('a_temp', engine, if_exists='replace', index=False)

# Now let's do the dissolve!
sql = f"""
DROP TABLE IF EXISTS merit_basins4_{id};

CREATE TABLE merit_basins4_{id} AS 
SELECT 
  T.mergeid::INTEGER as comid,
  ST_UNION(M.geom_simple) as geom
FROM merit_basins_{id} as M
JOIN a_temp AS T ON M.comid = T.comid
GROUP BY T.mergeid;
"""

cur = conn.cursor()
cur.execute(sql)
conn.commit()

