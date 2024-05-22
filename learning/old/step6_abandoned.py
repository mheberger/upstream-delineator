"""
Step #6: Prune the network from Step 5

Exactly following the process from step 4 but with the next level up of data.

"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import pandas as pd
from sqlalchemy import create_engine
import networkx as nx
import matplotlib.pyplot as plt

# Megabasin to process
id = 77

# Area in km² of a unit catchment below which we will eliminate it by merging it with its downstream neighbor
threshold = 150

# area of a unit catchment in km² that we consider very small, and will merge even if the target basin is already quite large
miniscule = 10

# How big is too big? Beyond which we should stop merging.
max_weight = 400

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

# Read in data from Step #5
sql = f"""
SELECT comid3, unitarea, nextdown
FROM merit_basins5_{id}
"""

# This was just for testing
# WHERE comid2 in (77037439,77037445,77037448,77037453,77037464,77037485,77037505,77037552,77037569,77038085,77038086,77038087);


# Read the result of the SQL query into an ordinary Pandas DataFrame
df = pd.read_sql(sql, conn)

# For debugging
#df.to_csv('before.csv')

df.set_index('comid3', inplace=True)

# Create a networkx object out of the rivers data.
G = nx.DiGraph()

for node_id, nextdown in df['nextdown'].items():
    #print(f"{comid}: {nextdown}")
    # Add node with comid2 as node ID
    G.add_node(node_id)
    # Add edge from comid2 to nextdown
    if nextdown > 0:
        G.add_edge(node_id, nextdown)

node_list = list(G.nodes)
node_list.sort()

# Add the weights to the nodes. The weight will be the area.
for node_id, area in df['unitarea'].items():
    G.nodes[node_id]['weight'] = area
    G.nodes[node_id]['absorbed_nodes'] = [node_id]

removed_nodes = []

# Here's the plan. Give nodes a property called "absorbed nodes" which is a list.
# Any time a node gets absorbed into it's downstream neighbor, pass on this list,
# appending its own id to the list.

def remove_nodes_below_threshold(G, threshold):
    """
    Remove nodes from the graph G with weights below the threshold.
    Maintain the full connectivity of the graph.
    """

    # Copy the original graph to avoid modifying it directly
    pruned_G = G.copy()

    # Get nodes to remove based on their weights
    nodes_to_remove = [node for node, data in G.nodes(data=True) if 'weight' in data and data['weight'] < threshold]

    # Remove nodes and connect predecessors to successors
    for node in nodes_to_remove:
        predecessors = list(pruned_G.predecessors(node))
        successors = list(pruned_G.successors(node))
        if len(successors) == 0:
            # We should not ever remove terminal nodes...
            continue

        # If we remove a node, add it's weight to the downstream node.
        weight = pruned_G.nodes[node]['weight']
        successor_node = successors[0]
        successor_weight = pruned_G.nodes[successor_node]['weight'] + weight
        if successor_weight < max_weight or weight < miniscule:
            pruned_G.nodes[successor_node]['weight'] = successor_weight

            absorbed_list = pruned_G.nodes[node]['absorbed_nodes']
            pruned_G.nodes[successor_node]['absorbed_nodes'].extend(absorbed_list)

            # Remove the node from the graph
            pruned_G.remove_node(node)
            removed_nodes.append(node)

            # Connect predecessors to successors
            for predecessor in predecessors:
                for successor in successors:
                    # Check if an edge already exists to avoid duplicates
                    if not pruned_G.has_edge(predecessor, successor):
                        pruned_G.add_edge(predecessor, successor)

    return pruned_G

def calculate_strahler_stream_order(graph):
    # Step 1: Determine upstream and downstream nodes
    upstream_nodes = {node: set() for node in graph.nodes()}
    downstream_nodes = {node: set() for node in graph.nodes()}
    for u, v in graph.edges():
        downstream_nodes[u].add(v)
        upstream_nodes[v].add(u)

    # Step 2: Assign initial stream orders
    for node in graph.nodes():
        graph.nodes[node]['strahler_order'] = 1

    # Step 3: Iterate through the nodes and update stream orders
    for node in nx.topological_sort(graph):
        upstream_orders = [graph.nodes[upstream]['strahler_order'] for upstream in upstream_nodes[node]]
        max_order = max(upstream_orders) if upstream_orders else 0
        order_count = upstream_orders.count(max_order)
        if order_count == 1:
            graph.nodes[node]['strahler_order'] = max_order
        else:
            graph.nodes[node]['strahler_order'] = max_order + 1

    return graph

# Now we can prune the network!
pruned_G = remove_nodes_below_threshold(G, threshold)

# Calculate the stream order
pruned_G = calculate_strahler_stream_order(pruned_G)

# The drawing was just for debugging.
DRAW = False
if DRAW:
    # Create figure and axes for side-by-side plots
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    # Plot original graph
    axs[0].set_title('Original Graph')
    pos = nx.drawing.layout.spring_layout(G, seed=22)
    node_weights = [G.nodes[node]['weight'] for node in G.nodes]
    scaling_factor = 6
    node_sizes = [weight * scaling_factor for weight in node_weights]
    mylabels = {node: G.nodes[node]['weight'] for node in G.nodes}
    # First one for weights, second for ids
    #nx.draw(G, pos, node_color="skyblue", labels=mylabels, node_size=node_sizes, ax=axs[0], font_size=8, font_family="Arial")
    nx.draw(G, pos, node_color="skyblue", with_labels=True, node_size=node_sizes, ax=axs[0], font_size=8, font_family="Arial")

    # Plot pruned graph
    axs[1].set_title('Pruned Graph')
    node_weights = [pruned_G.nodes[node]['weight'] for node in pruned_G.nodes]
    node_sizes = [weight * scaling_factor for weight in node_weights]
    mylabels = {node: pruned_G.nodes[node]['weight'] for node in pruned_G.nodes}
    mylabels = {node: pruned_G.nodes[node]['strahler_order'] for node in pruned_G.nodes}
    nx.draw(pruned_G, pos, node_color="skyblue", node_size=node_sizes, labels=mylabels, ax=axs[1], font_size=8, font_family="Arial")
    #nx.draw(pruned_G, pos, node_color="skyblue", node_size=node_sizes, with_labels=True, ax=axs[1], font_size=8, font_family="Arial")
    plt.tight_layout()
    plt.show()

print(f"Sum of unit catchment areas BEFORE merging: {sum(df['unitarea'])}")


for node in removed_nodes:
    df.loc[node, 'unitarea'] = 0
    df.loc[node, 'nextdown'] = -999

# Update the unit areas for the revised unit catchments, and the stream order
df['sorder'] = 0
for node in pruned_G.nodes:
    weight = pruned_G.nodes[node]['weight']
    df.loc[node, 'unitarea'] = weight
    sorder = pruned_G.nodes[node]['strahler_order']
    df.loc[node, 'sorder'] = sorder

# Use information from the graph to update the field mergeid
# Information is in the node's attribute 'absorbed_nodes'
for target_node in pruned_G.nodes:
    absorbed_list = pruned_G.nodes[target_node]['absorbed_nodes']
    for absorbed_node in absorbed_list:
        df.loc[absorbed_node, 'mergeid'] = target_node


# Replace null values in 'mergeid' column with index values. This is to facilitate
# the dissolve operation. In other words, we are setting the key for the dissolve.
# And since it does not change, the geometry will not be modified.
df.loc[df['mergeid'].isnull(), 'mergeid'] = df.index[df['mergeid'].isnull()]
df['mergeid'] = df['mergeid'].astype(int)

# Also update the field `nextdownid`. This information is contained in the pruned graph's edges.
for from_node, to_node in pruned_G.edges():
    df.loc[from_node, 'nextdown'] = to_node

# For debugging
#df.to_csv('after.csv')
print(f"Sum of unit catchment areas AFTER merging: {sum(df['unitarea'])}")

df.reset_index(inplace=True)

# Write DataFrame to PostgreSQL database
df.to_sql('a_temp', engine, if_exists='replace', index=False)

# Now let's do the dissolve!
sql = f"""
DROP TABLE IF EXISTS merit_basins6_{id};

CREATE TABLE merit_basins6_{id} AS 
SELECT 
  T.mergeid::INTEGER as comid3,
  SUM(M.num_merged) as num_merged,
  SUM(T.unitarea) as unitarea,
  MAX(T.sorder) as sorder,
  MAX(T.nextdown) as nextdown,
  ST_UNION(M.geom) as geom
FROM merit_basins4_{id} as M
JOIN a_temp AS T ON M.comid2 = T.comid3
GROUP BY T.mergeid;
"""

cur = conn.cursor()
cur.execute(sql)
conn.commit()

# Close the engine
engine.dispose()

