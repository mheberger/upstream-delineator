# Set of routines for simplifying and consolidating river networks

import pickle
import networkx as nx
from py.plot_network import draw_graph
from py.graph_tools import *
import numpy as np
from subbasins_config import VERBOSE
from scipy.stats import skew
import matplotlib.pyplot as plt

# Set to True to draw a bunch of network graphs. Mostly for debugging.
DRAW = False


def show_area_stats(G:nx.Graph):
    """
    Print some stats summarizing the area of the subbasins in our graph
    :param G:
    :return:
    """
    areas = [data['area'] for node, data in G.nodes(data=True)]

    # Step 3: Calculate statistics using NumPy
    n = len(areas)
    mean_area = np.mean(areas)
    median_area = np.median(areas)
    std_area = np.std(areas)
    cv = std_area / mean_area
    skewness = skew(areas)

    # Print the results
    print(f"n  | Media | Mean | Std. Dev. | CV  | Skew" )
    print(f"{n},   {median_area:.3g},   {mean_area:.3g},   {std_area:.3g},   {cv:.3g},   {skewness:.3g}")

    plt.hist(areas, bins='auto', edgecolor='black')
    plt.title('Histogram of Areas')
    plt.xlabel('Area')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()


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


def remove_entries_with_value(d, value_to_remove):
    """
    Remove all entries from the dictionary where the specified value is present.

    Parameters:
    d (dict): The dictionary to modify.
    value_to_remove: The value to remove from the dictionary.

    Returns:
    dict: A new dictionary with the specified entries removed.
    """
    return {k: v for k, v in d.items() if v != value_to_remove}


def step1(G: nx.DiGraph, threshold_area: int or float, MERGES: dict, rivers2merge, rivers2delete):
    # Step #1, merge "leaves", or unit catchments with order = 1.
    leaves = [n for n, attr in G.nodes(data=True) if attr.get('shreve_order') == 1 and 'custom' not in attr]
    for leaf in leaves:
        if G.has_node(leaf):
            successors = list(G.successors(leaf))
            if len(successors) > 0:
                successor = list(G.successors(leaf))[0]
                neighbors = list(G.predecessors(successor))
                neighbors.remove(leaf)
                for neighbor in neighbors:
                    if 'custom' not in G.nodes[neighbor] and G.nodes[neighbor]['shreve_order'] == 1:
                        merged_area = G.nodes[leaf]['area'] + G.nodes[neighbor]['area']
                        if merged_area < threshold_area:
                            G.nodes[leaf]['area'] = merged_area
                            G.remove_node(neighbor)
                            # This step is essential; we are tracking merges, but
                            # what happens when we delete a node that was previously a target?
                            MERGES = update_merges(MERGES, neighbor, leaf)
                            MERGES[neighbor] = leaf
                            rivers2delete.append(neighbor)
                            if neighbor in rivers2merge:
                                rivers2delete.extend(rivers2merge[neighbor])
                                del rivers2merge[neighbor]

    G = calculate_shreve_stream_order(G)
    return G, MERGES, rivers2merge, rivers2delete


def step2(G: nx.DiGraph, threshold_area: int or float, MERGES: dict, rivers2merge):
    """
    Step #2, merge "stem" nodes with their upstream neighbor.
    A stem node is one that has exactly 1 incoming edge, in other words it is not a junction
    If merging it with its upstream neighbor means the merged node has a size
    below the threshold, go ahead and merge it!

    """
    stems = [node for node in G.nodes if G.in_degree(node) == 1 and 'custom' not in G.nodes[node]]
    for stem in stems:
        if G.has_node(stem):
            predecessor = list(G.predecessors(stem))[0]
            if 'custom' in G.nodes[predecessor]:
                continue
            else:
                merged_area = G.nodes[stem]['area'] + G.nodes[predecessor]['area']
                if merged_area < threshold_area:
                    G.nodes[predecessor]['area'] = merged_area
                    successors = list(G.successors(stem))
                    if len(successors) > 0:
                        successor = successors[0]
                        G.add_edge(predecessor, successor)
                    G.remove_node(stem)
                    MERGES = update_merges(MERGES, stem, predecessor)
                    MERGES[stem] = predecessor

                    if stem in rivers2merge:
                        node_list = rivers2merge[stem]
                        del rivers2merge[stem]
                    else:
                        node_list = []

                    if predecessor in rivers2merge:
                        rivers2merge[predecessor].append(stem)
                    else:
                        rivers2merge[predecessor] = [stem]

                    rivers2merge[predecessor].extend(node_list)

    G = calculate_shreve_stream_order(G)
    return G, MERGES, rivers2merge


def step3(G: nx.DiGraph, threshold_length: int or float, MERGES: dict):
    """
    Step #3, where we eliminate the small "junction" nodes that occur where there are two or
    more confluences close to one another. If the river length in the node is very small, we
    can consider the confluences to occur in essentially the same location and eliminate the
    unit catchment. Requires special handling for the river reaches.
    """
    junctions = [node for node in G.nodes if G.nodes[node]['length'] < threshold_length and 'custom' not in G.nodes[node]]
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


def step4(G: nx.DiGraph, threshold_area: int or float, MERGES: dict, rivers2merge, rivers2delete):
    # Step #4, prune small solo "leaves", or unit catchments with order = 1 and area < threshold
    # TODO: Instead of merging the leaf with the downstream node, merge it with its neighbor,
    #   but only if there is just one (otherwise may try to merge with a polygon it does not touch)
    # This step will get us all the *candidate* leaves.
    # Those that have a shreve order = 1, area < threshold, and are NOT a custom node.
    leaves = [n for n, attr in G.nodes(data=True) if attr.get('shreve_order') == 1
              and attr.get('area') < threshold_area and 'custom' not in attr]

    # Actually, all of the candidate leaves are OK to merge?
    # Do not merge them with a downstream node if they have a neighbor!
    for leaf in leaves:
        if G.has_node(leaf):
            successor = list(G.successors(leaf))[0]
            neighbors = list(G.predecessors(successor))
            if len(neighbors) == 2:
                neighbors.remove(leaf)
                neighbor = neighbors[0]
                merged_area = G.nodes[leaf]['area'] + G.nodes[neighbor]['area']
                if merged_area < threshold_area:
                    G.nodes[neighbor]['area'] = merged_area
                    G = prune_node(G, leaf)
                    MERGES = update_merges(MERGES, leaf, neighbor)
                    MERGES[leaf] = neighbor
                    rivers2delete.append(leaf)
                    if leaf in rivers2merge:
                        rivers2delete.extend(rivers2merge[leaf])
                        del rivers2merge[leaf]

    G = calculate_shreve_stream_order(G)
    return G, MERGES, rivers2merge, rivers2delete


def consolidate_network(G: nx.DiGraph, threshold_area, threshold_length):
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
    MERGES = {}
    rivers2merge = {}
    rivers2delete = []

    if DRAW: draw_graph(G, filename='plots/test_net', title="Original Network")
    if VERBOSE:
        print('Iteration #1')

    # Step 1, merge leaves!
    G, MERGES, rivers2merge, rivers2delete = step1(G, threshold_area, MERGES, rivers2merge, rivers2delete)
    if DRAW: draw_graph(G, filename='plots/test_pruned', title="Pruned Network after Step 1")

    # Step 2, merge any stems
    G, MERGES, rivers2merge = step2(G, threshold_area, MERGES, rivers2merge)
    if DRAW: draw_graph(G, filename='plots/test_step2', title="Stems merged, after Step 2")

    # Step 3, collapse small junctions
    #G, MERGES = step3(G, threshold_length, MERGES)
    #if DRAW:  draw_graph(G, filename='plots/test_step3', title="Tiny junctions removed, after Step 3")

    # Iterate through the process until there are no more changes to the network?

    # Step #4, prune small solo leaves
    G, MERGES, rivers2merge, rivers2delete = step4(G, threshold_area, MERGES, rivers2merge, rivers2delete)
    if DRAW:  draw_graph(G, filename='plots/test_step4', title="Small solo leaves pruned, after Step 4")

    i = 1
    while True:
        previous_num_nodes = G.number_of_nodes()
        G, MERGES, rivers2merge, rivers2delete = step1(G, threshold_area, MERGES, rivers2merge, rivers2delete)
        G, MERGES, rivers2merge = step2(G, threshold_area, MERGES, rivers2merge)
        #G, MERGES = step3(G, threshold_length, MERGES)
        G, MERGES, rivers2merge, rivers2delete = step4(G, threshold_area, MERGES, rivers2merge, rivers2delete)
        num_nodes = G.number_of_nodes()
        i += 1
        if VERBOSE:
            print(f"Iteration #{i}")

        if num_nodes == previous_num_nodes:
            break

    if DRAW:  draw_graph(G, filename='plots/test_step5', title="After iteration, Step 5")
    if VERBOSE: show_area_stats(G)

    return G, MERGES, rivers2merge, rivers2delete


def main():
    fname = 'output/ice2_graph.pkl'
    G = pickle.load(open(fname, "rb"))
    # Set a threshold for merging (e.g., 100 for the sum of areas)
    threshold_area = 300
    threshold_length = 2
    G, MERGES, rivers2merge, rivers2delete = consolidate_network(G, threshold_area=threshold_area,
                                                                 threshold_length=threshold_length)

if __name__ == "__main__":
