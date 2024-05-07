# A couple of functions to calculate the stream order (Strahler and Shreve)
# These are useful for visualizing and processing the river network data.

import networkx as nx


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


def calculate_shreve_stream_order(graph):
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
