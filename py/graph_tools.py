# A set of functions for manipulating river network data as Python NetworkX Graphs.
# Functions to create a graph, calculate the stream order (Strahler and Shreve)
# These are useful for visualizing and processing the river network data.

import networkx as nx
from pandas import DataFrame


def make_river_network(df: DataFrame, terminal_node=None) -> nx.Graph:
    """
    Creates a network graph of our river network data.
    Input is a Pandas DataFrame. The index should be a unique id for the node in the network (river reach, unit catchment)
    and there must be a field `nextdown`, which is id the downstream neighbor, or direction of flow.

    TODO: Do not add an outgoing edge to the river network's terminal node?
    """
    G = nx.DiGraph()

    # Populate the graph's nodes and edges.
    for node_id, nextdown in df['nextdown'].items():
        # Add node with comid as node ID
        G.add_node(node_id)
        # Add edge from comid to nextdown
        if nextdown > 0 and node_id != terminal_node:
            G.add_edge(node_id, nextdown)

    # Simple way to do make sure the terminal node is not added?

    return G


def calculate_strahler_stream_order(graph):
    """
        Calculates the Strahler stream order for a river network that is
        represented as a NetworkX acyclic directed graph.
        Adds the attribute 'strahler_order' for nodes in the network.

        Arguments:
            a Graph
        Returns:
            a Graph

        """
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


def calculate_shreve_stream_order(graph: nx.Graph) -> nx.Graph:
    """
    Calculates the Shreve stream order for a river network that is
    represented as a NetworkX acyclic directed graph.
    Adds the attribute 'shreve_order' for nodes in the network.

    Arguments:
        a Graph
    Returns:
        a Graph

    """
    upstream_nodes = {node: set() for node in graph.nodes()}
    downstream_nodes = {node: set() for node in graph.nodes()}
    for u, v in graph.edges():
        downstream_nodes[u].add(v)
        upstream_nodes[v].add(u)

    for node in graph.nodes():
        graph.nodes[node]['shreve_order'] = 1

    for node in nx.topological_sort(graph):
        if upstream_nodes[node]:
            graph.nodes[node]['shreve_order'] = max([graph.nodes[upstream]['shreve_order'] for upstream in upstream_nodes[node]]) + 1

    return graph


def calculate_num_incoming(G):
    """calculates the number of upstream or incoming edges to each node
    and puts that in an attribute"""

    for node in G.nodes():
        num_incoming = G.in_degree(node)
        G.nodes[node]['num_incoming'] = num_incoming

    return G


def add_distance_from_outlet(G):
    """
    Adds an attribute `distance` to each node, indicating the number of edges
    between the node and the outlet (terminal node)

    I just made up this attribute.
    It is kind of the inverse of the shreve order.
    So perhaps it is not really necessary.
    I initially thought it could be useful for plotting the network, until I discovered GraphViz.
    Code from ChatGPT
    TODO: Not thoroughly tested!!!
    """

    # Ensure the graph is a Directed Acyclic Graph (DAG)
    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("The graph must be a Directed Acyclic Graph (DAG)")

    # Perform topological sort
    topo_order = list(nx.topological_sort(G))

    # Initialize distance attribute for all nodes
    nx.set_node_attributes(G, 0, 'distance')

    # Compute distance for each node
    for node in topo_order:
        for pred in G.predecessors(node):
            G.nodes[node]['distance'] = max(G.nodes[node]['distance'], G.nodes[pred]['distance'] + 1)

    return G


def insert_node(G: nx.Graph, node, comid) -> nx.Graph:
    """
    Custom function to insert a new node in my flow network graph at a given location.

    To insert new river reaches and unit catchments into the existing flow network, it
    was helpful to create a Python NetworkX Graph.

    Adding an outlet point means inserting a new node to the graph.
    We selectively add and remove edges to maintain the connectivity of the graph.

    The logic is different if we are inserting a node
    into a unit catchment that is a "leaf" node, or one with a Strahler order 1
    or into a "stem" node, one with Strahler order > 1.
    """

    # Get the Strahler order of the node we're inserting into.
    order = G.nodes[comid]['strahler_order']

    if order == 1:
        # LEAF NODE
        # Add the new node to the network
        G.add_node(node)
        # Set some attributes for the node that will be useful later
        G.nodes[node]['new'] = True
        G.nodes[node]['type'] = 'leaf'
        # Add the network connection (edge) from the new node to the target unit catchment with comid = comid.
        G.add_edge(node, comid)

    else:
        # STEM NODE
        # Step 1: Add the new node
        G.add_node(node)
        G.nodes[node]['new'] = True
        G.nodes[node]['type'] = 'stem'

        # Step 2: Find incoming edges and remove them.
        predecessors = list(G.predecessors(comid))
        incoming_edges = list(G.in_edges(comid))
        for u, v in incoming_edges:
            G.remove_edge(u, v)

        # Step 3: Add an edge from the new node to the comid (node it is being inserted upstream of)
        G.add_edge(node, comid)

        # Step #4: Add new incoming edges to the new node
        for predecessor in predecessors:
            G.add_edge(predecessor, node)

    return G

