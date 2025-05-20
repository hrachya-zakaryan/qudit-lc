import numpy as np
import itertools
import networkx as nx
import matplotlib.pyplot as plt

from networkx.drawing.nx_agraph import graphviz_layout
import pygraphviz as pgv
import networkx as nx
import matplotlib.pyplot as plt

def draw_graphs(encoded_graphs, n, d, cols=3):
    """
    Draws multiple graphs given their bit encoding.

    Parameters:
        encoded_graphs (list of int): A list of bit-encoded graphs.
        n (int): Number of vertices in the graphs.
        d (int): Modulo value for weights (used for decoding).
        cols (int): Number of columns in the subplot grid (default: 3).
    """
    num_graphs = len(encoded_graphs)
    rows = (num_graphs + cols - 1) // cols  # Compute number of rows needed

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))  # Create subplots
    axes = axes.flatten() if num_graphs > 1 else [axes]  # Flatten axes for easy iteration

    for i, encoded_graph in enumerate(encoded_graphs):
        ax = axes[i]
        adjacency_matrix = bitpack_decode(encoded_graph, n, d)
        
        G = nx.Graph()
        for v in range(n):
            G.add_node(v, label=f'Node {v}')
            for u in range(v + 1, n):
                if adjacency_matrix[v, u] > 0:
                    G.add_edge(v, u, weight=adjacency_matrix[v, u])

        pos = nx.spring_layout(G)  
        edge_labels = {(u, v): f"{w}" for u, v, w in G.edges.data("weight")}
        
        nx.draw(G, pos, ax=ax, with_labels=True, node_color="lightblue", node_size=500, font_weight="bold")
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax)
        
        ax.set_title(f"Graph {i+1}")
        ax.axis("off")  # Hide axis

    # Hide any unused subplot axes
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()


layout_engines = ['dot', 'neato', 'circo', 'fdp', 'sfdp', 'twopi']

def edge_crossings(pos, edges):
    crossings = 0
    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            (u1, v1), (u2, v2) = edges[i], edges[j]
            if len({u1, v1, u2, v2}) < 4:
                continue  # skip if edges share a node
            a, b = pos[u1], pos[v1]
            c, d = pos[u2], pos[v2]
            if lines_intersect(a, b, c, d):
                crossings += 1
    return crossings

def lines_intersect(p1, p2, q1, q2):
    def ccw(a, b, c):
        return (c[1]-a[1]) * (b[0]-a[0]) > (b[1]-a[1]) * (c[0]-a[0])
    return ccw(p1, q1, q2) != ccw(p2, q1, q2) and ccw(p1, p2, q1) != ccw(p1, p2, q2)

def layout_score(pos, graph):
    crossings = edge_crossings(pos, list(graph.edges()))
    spread = np.var([p[0] for p in pos.values()] + [p[1] for p in pos.values()])
    return crossings * 10 - spread  # fewer crossings and more spread = better

def best_layout(graph):
    best_pos = None
    best_score = float('inf')
    for engine in layout_engines:
        try:
            pos = graphviz_layout(graph, prog=engine)
            score = layout_score(pos, graph)
            if score < best_score:
                best_score = score
                best_pos = pos
        except:
            continue
    return best_pos if best_pos else nx.spring_layout(graph)


def circular_subgraph_layout(n, center, radius):
    """ Arrange n points in a circular layout inside a given center and radius. """
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return {i: center + radius * np.array([np.cos(angle), np.sin(angle)]) for i, angle in enumerate(angles)}

def draw_representatives(encoded_graphs, d, cols=3):
    """
    Draws multiple graphs given their bit encoding.

    Parameters:
        encoded_graphs (list of int): A list of bit-encoded graphs.
        n (int): Number of vertices in the graphs.
        d (int): Modulo value for weights (used for decoding).
        cols (int): Number of columns in the subplot grid (default: 3).
    """
    num_graphs = len(encoded_graphs)
    rows = (num_graphs + cols - 1) // cols  # Compute number of rows needed

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 1, rows * 1))  # Create subplots
    axes = axes.flatten() if num_graphs > 1 else [axes]  # Flatten axes for easy iteration

    for i, encoded_graph in enumerate(encoded_graphs):
        ax = axes[i]
        n=encoded_graph[1][0]
        adjacency_matrix = bitpack_decode(encoded_graph[0], n, d)
        
        subG = nx.from_numpy_array(adjacency_matrix)  # Convert adjacency matrix to NetworkX graph
        #pos = circular_subgraph_layout(len(subG.nodes), np.array(center), radius * 0.7)  # Keep subgraph size, shrink nodes

        # Get weighted edges
        edges = [(u, v) for u, v in subG.edges()]

        # Set edge colors based on the weight
        edge_colors = ['black' if adjacency_matrix[u, v] == 1 else 'red' for u, v in edges]
        
        #pos = nx.nx_agraph.graphviz_layout(subG, prog='circo')
        #pos = nx.circular_layout(subG)
        node_order, min_cross = best_permutation(subG)
        pos = circular_positions(node_order)
        # Draw subgraph edges with the appropriate color
        nx.draw_networkx_edges(subG, pos, ax=ax, edge_color=edge_colors, alpha=1, width=2)
        
        # Draw subgraph nodes (light blue fill, black outline)
        nx.draw_networkx_nodes(subG, pos, ax=ax, node_size=40, node_color='blue', edgecolors='black', linewidths=1)
        
        ax.set_title(f"Graph {i+1}")
        ax.axis("off")  # Hide axis
        ax.set_aspect('equal')

    # Hide any unused subplot axes
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig("Representatives.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def circular_positions(nodes):
    """Assign fixed circular positions to nodes."""
    n = len(nodes)
    angle_step = 2 * np.pi / n
    return {node: (np.cos(i * angle_step), np.sin(i * angle_step)) for i, node in enumerate(nodes)}

def count_crossings(G, pos):
    """Count the number of edge crossings."""
    crossings = 0
    edges = list(G.edges())

    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])

    def intersect(e1, e2):
        a, b = pos[e1[0]], pos[e1[1]]
        c, d = pos[e2[0]], pos[e2[1]]
        return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)

    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            if len(set(edges[i]) & set(edges[j])) == 0:  # disjoint edges only
                if intersect(edges[i], edges[j]):
                    crossings += 1
    return crossings

def best_permutation(G):
    """Find permutation with minimal crossings."""
    nodes = list(G.nodes())
    min_crossings = float('inf')
    best_order = nodes

    for perm in itertools.permutations(nodes):
        pos = circular_positions(perm)
        crossings = count_crossings(G, pos)
        if crossings < min_crossings:
            min_crossings = crossings
            best_order = perm

    return best_order, min_crossings

def draw_graph(encoded_graph, n, d):
    """
    Draws a graph given its bit encoding.
    
    Parameters:
        encoded_graph (int): The bit-encoded graph.
        n (int): Number of vertices in the graph.
        d (int): Modulo value for weights (used for decoding).
    """
    # Decode the graph to an adjacency matrix
    adjacency_matrix = bitpack_decode(encoded_graph, n, d)
    
    # Create a NetworkX graph
    G = nx.Graph()
    
    # Add nodes and edges to the graph
    for i in range(n):
        G.add_node(i, label=f'Node {i}')
        for j in range(i + 1, n):
            if adjacency_matrix[i, j] > 0:
                G.add_edge(i, j, weight=adjacency_matrix[i, j])
    
    # Draw the graph with labels
    pos = nx.spring_layout(G)  # Use a spring layout for visualization
    edge_labels = {(u, v): f"{d['weight']}" for u, v, d in G.edges(data=True)}
    nx.draw(G, pos, with_labels=True, node_color="lightblue", node_size=500, font_weight="bold")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    
    # Show the plot
    plt.title("Graph Visualization")
    plt.show()



def bitpack_encode(matrix, d):
    n = matrix.shape[0]
    bit_length = (d-1).bit_length()  # Number of bits per weight
    packed = 0
    shift = 0
    
    # Iterate through the upper triangular part (excluding diagonal)
    for i in range(n):
        for j in range(i + 1, n):
            weight = matrix[i, j]
            packed |= (weight << shift)  # Shift and add weight to packed
            shift += bit_length  # Update shift for the next weight
    
    return int(packed)

def bitpack_decode(packed, n, d):
    bit_length = (d-1).bit_length()
    matrix = np.zeros((n, n), dtype=int)
    shift = 0
    
    # Decode weights from packed integer
    for i in range(n):
        for j in range(i + 1, n):
            weight = (packed >> shift) & ((1 << bit_length) - 1)  # Extract bits
            matrix[i, j] = weight
            matrix[j, i] = weight  # Reflect symmetry
            shift += bit_length  # Move to the next weight
    
    return matrix


# Local scaling

def local_scaling(G, v: int, k: int, d: int):
    n = G.shape[0]
    for u in range(n):
        if u != v:
            G[v, u] = (k * G[v, u]) % d
            G[u, v] = G[v, u]  # Symmetric

# Local complementation


def local_complementation(G, v: int, k: int, d: int):
    n = G.shape[0]
    for u in range(n):
        for w in range(n):
            if u!=w and u != v and w != v and G[v, u] != 0 and G[v, w] != 0:
                G[u, w] = (G[u, w] + k * G[v, u] * G[v, w]) % d
