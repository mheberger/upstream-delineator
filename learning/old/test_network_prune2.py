"""
Further experiments with using a network pruning algorithm.

Here, I have adapted the code from ChatGPT so that when I remove a node
from the network, its weight is added to its downstream neighbor. 
And I set a second threshold for the max. weight for a node. 

This ended up partially forming the basis for the solution I used. 

"""


import networkx as nx
import matplotlib.pyplot as plt


def remove_nodes_below_threshold(G, threshold):
    """
    Remove nodes from the graph G with weights below the threshold.
    Maintain the full connectivity of the graph.
    """
    max_weight = 25
    
    # Copy the original graph to avoid modifying it directly
    pruned_G = G.copy()

    # Get nodes to remove based on their weights
    nodes_to_remove = [node for node, data in G.nodes(data=True) if 'weight' in data and data['weight'] < threshold]

    # Remove nodes and connect predecessors to successors
    for node in nodes_to_remove:
        predecessors = list(pruned_G.predecessors(node))
        successors = list(pruned_G.successors(node))

        # If we remove a node, add it's weight to the downstream node.
        weight = pruned_G.nodes[node]['weight']
        successor_node = successors[0]
        successor_weight = pruned_G.nodes[successor_node]['weight'] + weight
        if successor_weight < max_weight:
            pruned_G.nodes[successor_node]['weight'] = successor_weight

            # Remove the node from the graph
            pruned_G.remove_node(node)

            # Connect predecessors to successors
            for predecessor in predecessors:
                for successor in successors:
                    # Check if an edge already exists to avoid duplicates
                    if not pruned_G.has_edge(predecessor, successor):
                        pruned_G.add_edge(predecessor, successor)

    return pruned_G


# Example usage
# Create an example acyclic directed graph with weighted edges
G = nx.DiGraph()
G.add_edges_from([
    (86, 85),
    (87, 86),
    (48, 87),
    (64, 48),
    (52, 64),
    (39, 87),
    (45, 39),
    (53, 45),
    (99, 53),
    (69, 99),
    (5, 99),
]
)

G.nodes[85]['weight'] = 10
G.nodes[86]['weight'] = 10
G.nodes[87]['weight'] = 5
G.nodes[39]['weight'] = 5
G.nodes[45]['weight'] = 5
G.nodes[53]['weight'] = 5
G.nodes[99]['weight'] = 10
G.nodes[69]['weight'] = 5
G.nodes[5]['weight']  = 10
G.nodes[48]['weight'] = 10
G.nodes[64]['weight'] = 10
G.nodes[52]['weight'] = 10


# Define the threshold
threshold = 6

# Prune the graph
pruned_G = remove_nodes_below_threshold(G, threshold)

# Print pruned graph edges and their weights
for u, v, data in pruned_G.edges(data=True):
    print(f"Edge: {u} -> {v}")

# Create figure and axes for side-by-side plots
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# Plot original graph
axs[0].set_title('Original Graph')
pos = nx.drawing.layout.spring_layout(G, seed=22)

# Extract node weights
node_weights = [G.nodes[node]['weight'] for node in G.nodes]
scaling_factor = 80
node_sizes = [weight * scaling_factor for weight in node_weights]


# Draw the river network with weights
nx.draw(G, pos, with_labels=True, node_size=node_sizes, ax=axs[0])

# Plot pruned graph
axs[1].set_title('Pruned Graph')
node_weights = [pruned_G.nodes[node]['weight'] for node in pruned_G.nodes]
node_sizes = [weight * scaling_factor for weight in node_weights]
nx.draw(pruned_G, pos, with_labels=True, node_size=node_sizes, ax=axs[1])

plt.tight_layout()
plt.show()
