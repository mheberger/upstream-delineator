import networkx as nx
import graphviz
from PIL import Image


def draw_graph(G: nx.DiGraph, filename: str, vertical=False):
    """
    Plots a NetworkX directed acyclic graph (dag)
    using the GraphViz library. Note that you need to install this software and it needs to be
    on the system path.

    Args:
        G: the river network directed graph.
        filename: the filename (bare, no extension) to create
        vertical: if True, will arrange plot from top to bottom,
          Default is left to right.
    """
    # Initialize a new directed graph in Graphviz
    if vertical:
        dot = graphviz.Digraph()
    else:
        dot = graphviz.Digraph(graph_attr={'rankdir': 'LR'})

    # Add nodes with distance attributes as labels
    for node in G.nodes():
        if str(node)[0:2] == '27':
            label = str(node)[-4:]
        else:
            label = f"{node}"

        if G.nodes[node].get('new'):
            dot.node(str(node), label=label, style='filled', fillcolor='lightblue', fontname='Arial')
        else:
            dot.node(str(node), label=label, fontname='Arial')

    # Add edges
    for u, v in G.edges():
        dot.edge(str(u), str(v))

    # Save the graph as a file and render it
    dot.render(filename, format='png', cleanup=True)

    # Display the graph using PIL
    #image = Image.open(filename + ".png")
    #image.show()


if __name__ == "__main__":
    # Example usage:
    # Create a directed acyclic graph (DAG)
    G = nx.DiGraph()
    G.add_edges_from([(1, 2), (1, 3), (3, 4), (2, 4), (4, 5)])

    # Plot the DAG
    draw_graph(G, "my_dag")
