# Set of routines for simplifying and consolidating river networks

import pickle
import networkx as nx
from py.plot_network import draw_graph
from py.graph_tools import *

# Set to True to draw a bunch of network graphs. Mostly for debugging.
DRAW = False


def find_keys_by_value(dictionary, value):
    keys = []
    for key, val in dictionary.items():
        if val == value:
            keys.append(key)
    return keys


def update_merges(MERGES, node, target):
    """
    Update the merges dictionary as necessary whenever we add a new node: target pair

    :param MERGES: a dictionary mapping node: target
    :param node: the node that is being pruned, or removed from the network
    :param target: the target node, to which the area will be merged
    :return: MERGES
    """
    if node in MERGES.values():
        keys = find_keys_by_value(MERGES, node)
        for key in keys:
            MERGES[key] = target

    return MERGES


def step1(G: nx.DiGraph, threshold_area: int or float, MERGES: dict):
    # Step #1, merge "leaves", or unit catchments with order = 1.
    leaves = [n for n, attr in G.nodes(data=True) if attr.get('shreve_order') == 1 and 'type' not in attr]
    for leaf in leaves:
        if G.has_node(leaf):
            successor = list(G.successors(leaf))[0]
            neighbors = list(G.predecessors(successor))
            neighbors.remove(leaf)
            for neighbor in neighbors:
                if 'type' not in G.nodes[neighbor] and G.nodes[neighbor]['shreve_order'] == 1:
                    merged_area = G.nodes[leaf]['area'] + G.nodes[neighbor]['area']
                    if merged_area < threshold_area:
                        G.nodes[leaf]['area'] = merged_area
                        G.remove_node(neighbor)
                        # This step is essential; we are tracking merges, but
                        # what happens when we delete a node that was previously a target?
                        MERGES = update_merges(MERGES, neighbor, leaf)
                        MERGES[neighbor] = leaf

    G = calculate_shreve_stream_order(G)
    return G, MERGES


def step2(G: nx.DiGraph, threshold_area: int or float, MERGES: dict):
    """
    Step #2, merge "stems", these are nodes that have exactly 1 incoming edges
    If merging it with its downstream neighbor means the merged node has a size
    below the threshold, go ahead and merge it!

    """
    stems = [node for node in G.nodes if G.in_degree(node) == 1 and 'type' not in G.nodes[node]]
    for stem in stems:
        if G.has_node(stem):
            successor = list(G.successors(stem))[0]
            merged_area = G.nodes[stem]['area'] + G.nodes[successor]['area']
            if merged_area < threshold_area:
                predecessor = list(G.predecessors(stem))
                G.nodes[successor]['area'] = merged_area
                if predecessor:
                    if successor:
                        G.add_edge(predecessor[0], successor)
                G.remove_node(stem)
                MERGES = update_merges(MERGES, stem, successor)
                MERGES[stem] = successor

    G = calculate_shreve_stream_order(G)
    return G, MERGES


def step3(G: nx.DiGraph, threshold_length: int or float, MERGES: dict):
    """
    Step #3, where we eliminate the small "junction" nodes that occur where there are two or
    more confluences close to one another. If the river length in the node is very small, we
    can consider the confluences to occur in essentially the same location and eliminate the
    unit catchment. Requires special handling for the river reaches.
    """
    junctions = [node for node in G.nodes if G.nodes[node]['length'] < threshold_length and 'type' not in G.nodes[node]]
    for junction in junctions:
        if G.has_node(junction):
            successor = list(G.successors(junction))[0]
            merged_area = G.nodes[junction]['area'] + G.nodes[successor]['area']
            G.nodes[successor]['area'] = merged_area
            G = prune_node(G, junction)
            MERGES = update_merges(MERGES, junction, successor)
            MERGES[junction] = successor

    G = calculate_shreve_stream_order(G)
    return G, MERGES


def step4(G: nx.DiGraph, threshold_area: int or float, MERGES: dict):
    # Step #4, prune small solo "leaves", or unit catchments with order = 1 and area < threshold

    # This step will get us all the candidate leaves.
    leaves = [n for n, attr in G.nodes(data=True) if attr.get('shreve_order') == 1
              and attr.get('area') < threshold_area and 'type' not in attr]

    for leaf in leaves:
        if G.has_node(leaf):
            successor = list(G.successors(leaf))[0]
            merged_area = G.nodes[leaf]['area'] + G.nodes[successor]['area']
            if merged_area < threshold_area:
                G.nodes[successor]['area'] = merged_area
                G = prune_node(G, leaf)
                MERGES = update_merges(MERGES, leaf, successor)
                MERGES[leaf] = successor

    G = calculate_shreve_stream_order(G)
    return G, MERGES


def consolidate_network(G: nx.DiGraph, MERGES: dict, threshold_area, threshold_length):
    """
    Consolidates the nodes in a river network graph, merging nodes to make them larger, while
    preserving the overall shape and connectivity of the graph.
    Saves a record of nodes that have been merged in the dictionary MERGES that we will later
    use to dissolve their geometries (unit catchment polygons) in GeoPandas.

    Inputs: G: graph of the river network. Nodes should have the attributes 'area' and 'length'
      MERGES: a dictionary, can be empty, or may already exist if we call this funtion more than once.

    Outputs: returns the inputs. G should be smaller (fewer nodes) and MERGES will be bigger (more entries)
    """

    # Make sure that the network has the Shreve stream order attribute for every node
    G = calculate_shreve_stream_order(G)
    if DRAW: draw_graph(G, filename='plots/test_net', title="Original Network")

    # Step 1, merge leaves!
    G, MERGES = step1(G, threshold_area, MERGES)
    if DRAW: draw_graph(G, filename='plots/test_pruned', title="Pruned Network after Step 1")

    # Step 2, merge any stems
    G, MERGES = step2(G, threshold_area, MERGES)
    if DRAW: draw_graph(G, filename='plots/test_step2', title="Stems merged, after Step 2")

    # Step 3, collapse small junctions
    G, MERGES = step3(G, threshold_length, MERGES)
    if DRAW:  draw_graph(G, filename='plots/test_step3', title="Tiny junctions removed, after Step 3")

    # Step #4, prune small solo leaves
    G, MERGES = step4(G, threshold_area, MERGES)
    if DRAW:  draw_graph(G, filename='plots/test_step4', title="Small solo leaves pruned, after Step 4")

    # Iterate through the process until there are no more changes to the network?

    while True:
        previous_num_nodes = G.number_of_nodes()
        G, MERGES = step1(G, threshold_area, MERGES)
        G, MERGES = step2(G, threshold_area, MERGES)
        G, MERGES = step3(G, threshold_length, MERGES)
        G, MERGES = step4(G, threshold_area, MERGES)
        num_nodes = G.number_of_nodes()
        if num_nodes == previous_num_nodes:
            break

    if DRAW:  draw_graph(G, filename='plots/test_step5', title="After iteration, Step 5")

    return G, MERGES


if __name__ == "__main__":
    fname = 'output/ice2_graph.pkl'
    G = pickle.load(open(fname, "rb"))
    MERGES = {}
    # Set a threshold for merging (e.g., 100 for the sum of areas)
    threshold_area = 300
    threshold_length = 2
    G, MERGES = consolidate_network(G, MERGES, threshold_area=100, threshold_length=2)
