"""
Learning: Testing using NetworkX, and an automatic network pruning function

Idea is to automatically remove nodes from a network where an attribute has 
a certain value. For example, the weight attribute falls below a threshold. 

Code was largely from ChatGPT. 

Was useful for learning, and some code I later adapted. 

"""

import networkx as nx
import matplotlib.pyplot as plt


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
        predecessors = list(G.predecessors(node))
        successors = list(G.successors(node))

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
    (2, 1),
    (3, 1),
    (4, 2),
    (5, 2),
    (6, 3),
    (7, 4),
    (8, 6),
    (9, 7),
    (10, 6),
    (1, 11),
    (12, 11),
    (11, 13)
]
)

G.nodes[1]['weight'] = 9
G.nodes[2]['weight'] = 3
G.nodes[3]['weight'] = 8
G.nodes[4]['weight'] = 9
G.nodes[5]['weight'] = 7
G.nodes[6]['weight'] = 6
G.nodes[7]['weight'] = 8
G.nodes[8]['weight'] = 3
G.nodes[9]['weight'] = 9
G.nodes[10]['weight'] = 8
G.nodes[11]['weight'] = 1
G.nodes[12]['weight'] = 8
G.nodes[13]['weight'] = 8


# Define the threshold
threshold = 4

# Prune the graph
pruned_G = remove_nodes_below_threshold(G, threshold)

# Print pruned graph edges and their weights
for u, v, data in pruned_G.edges(data=True):
    print(f"Edge: {u} -> {v}")

# Create figure and axes for side-by-side plots
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# Plot original graph
axs[0].set_title('Original Graph')
pos = nx.drawing.layout.spring_layout(G, seed=42)

# Extract node weights
node_weights = [G.nodes[node]['weight'] for node in G.nodes]
scaling_factor = 40
node_sizes = [weight * scaling_factor for weight in node_weights]


# Draw the river network with stream orders
nx.draw(G, pos, with_labels=True, node_size=node_sizes, ax=axs[0])

# Plot pruned graph
axs[1].set_title('Pruned Graph')
nx.draw(pruned_G, pos, with_labels=True, ax=axs[1])

plt.tight_layout()
plt.show()
