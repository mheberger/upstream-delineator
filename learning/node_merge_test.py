# Try this algorithm, from CoPilot.
# Just wondering whether it would give a better answer than what I got before. 
# The answer is no, it does not. This is basically crap.

import networkx as nx
from plot_network import draw_graph
import pickle

# Assuming you have your original graph 'G' with nodes and weights
# Define your weight threshold (e.g., 100)

def merge_adjacent_nodes(graph, threshold):
    simplified_graph = nx.DiGraph()

    for node in graph.nodes():
        # Check if the combined weight of adjacent nodes is below the threshold
        outgoing_edges = graph.out_edges(node)
        incoming_edges = graph.in_edges(node)
        total_weight = sum(graph[u][v]['weight'] for u, v in outgoing_edges)
        total_weight += sum(graph[u][v]['weight'] for u, v in incoming_edges)

        if total_weight < threshold:
            # Add the merged node to the simplified graph
            simplified_graph.add_node(node, weight=total_weight)
            # Update outgoing and incoming edges
            for u, v in outgoing_edges:
                simplified_graph.add_edge(node, v)
            for u, v in incoming_edges:
                simplified_graph.add_edge(u, node)
        else:
            # If weight exceeds threshold, add the node as-is
            simplified_graph.add_node(node, weight=graph.nodes[node]['weight'])

    return simplified_graph


def main():
    G = pickle.load(open("graph.pkl", "rb"))
    draw_graph(G, 'plots/before.png')
    

if __name__ == "__main__":
    main()


# Now 'simplified_graph' contains the simplified river network
# Verify that it preserves all information and has the same total weight

# Note: Adjust the threshold value based on your specific requirements.
