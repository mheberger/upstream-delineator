import networkx as nx
import matplotlib.pyplot as plt


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

# Example usage
# Create a river network graph (replace with your graph creation code)
river_network = nx.DiGraph()
river_network.add_edges_from([
    (86, 85),
    (87, 86),
    (48, 87),
    (64, 48),
    (52, 64),
    (65, 64),
    (39, 87),
    (45, 39),
    (53, 45),
    (99, 53),
    (69, 99),
    (5, 99),
]
)

# Calculate Strahler stream order
river_network_with_order = calculate_shreve_stream_order(river_network)

# Print the stream orders for each node
for node in river_network_with_order.nodes():
    #print(f"Node {node}: Strahler order {river_network_with_order.nodes[node]['strahler_order']}")
    print(f"Node {node}: Shreve order {river_network_with_order.nodes[node]['shreve_order']}")

# Plot the river network with stream order labels
pos = nx.spring_layout(river_network_with_order)  # positions for all nodes

plt.figure(figsize=(10, 6))
labels = nx.get_node_attributes(river_network_with_order, 'shreve_order')
nx.draw(river_network_with_order, pos, labels=labels, node_size=2000, node_color='skyblue', font_size=10)


plt.title("River Network with Shreve Stream Order")
#plt.title("River Network with Strahler Stream Order")
plt.show()
