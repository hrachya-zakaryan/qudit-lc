import numpy as np
import itertools
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from networkx.drawing.nx_agraph import graphviz_layout
import pygraphviz as pgv
import networkx as nx
import matplotlib.pyplot as plt


def circular_subgraph_layout(n, center, radius):
    """ Arrange n points in a circular layout inside a given center and radius. """
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return {i: center + radius * np.array([np.cos(angle), np.sin(angle)]) for i, angle in enumerate(angles)}

def draw_subgraph_inside_circle(ax, center, radius, adjacency_matrix):
    """ Draws a small subgraph inside a circular boundary. """
    subG = nx.from_numpy_array(adjacency_matrix) 
    pos = circular_subgraph_layout(len(subG.nodes), np.array(center), radius * 0.7)
    edges = [(u, v) for u, v in subG.edges()]
    edge_colors = ['black' if adjacency_matrix[u, v] == 1 else 'red' for u, v in edges]
    nx.draw_networkx_edges(subG, pos, ax=ax, edge_color=edge_colors, alpha=1, width=1.3)
    nx.draw_networkx_nodes(subG, pos, ax=ax, node_size=23, node_color='blue', edgecolors='black', linewidths=0.6)


def plot_graph(G, n, d, seed, orbit):
    """
    Plot a graph G using NetworkX with each node as a circle containing a graph.

    Args:
        G(nx.Graph): Graph in networkx form
        n(int): Number of vertices
        d(int): Local dimension
        seed(int): Seed for the spring_layout
        orbit(int): Orbit number 
    """
    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(G, seed=seed, k=1/np.sqrt(len(G.nodes)))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='blue', alpha=0.7, width=0.8)
    for node in G.nodes():
        x, y = pos[node]
        circle = plt.Circle((x, y), 0.1, color='white', ec='black', lw=1.2)
        ax.add_patch(circle)
        adjacency_matrix = bitpack_decode(int(node), n, d)
        draw_subgraph_inside_circle(ax, (x, y), 0.08, adjacency_matrix)

    ax.set_xlim(min(x for x, y in pos.values()) - 0.2, max(x for x, y in pos.values()) + 0.2)
    ax.set_ylim(min(y for x, y in pos.values()) - 0.2, max(y for x, y in pos.values()) + 0.2)
    ax.set_aspect('equal')
    ax.axis('off')

    legend_elements = [
        mlines.Line2D([], [], color='black', lw=1.3, label='Edge Weight = 1'),
        mlines.Line2D([], [], color='red', lw=1.3, label='Edge Weight = 2')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=7, frameon=True)
    plt.tight_layout()
    plt.savefig(f"n{n}d{d}_o{orbit}_s{seed}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_graphs(encoded_graphs, n, d, cols=3):
    """
    Draws multiple graphs given their bit encoding.

    Parameters:
        encoded_graphs (list[int]): A list of bit-encoded graphs.
        n (int): Number of vertices in the graphs.
        d (int): Modulo value for weights (used for decoding).
        cols (int): Number of columns in the subplot grid (default: 3).
    """

    num_graphs = len(encoded_graphs)
    rows = (num_graphs + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4)) 
    axes = axes.flatten() if num_graphs > 1 else [axes]

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
        ax.axis("off") 

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()

def draw_representatives(encoded_graphs, d, cols=3):
    """
    Draws the representatives of the orbits.

    Parameters:
        encoded_graphs (list[int]): A list of bit-encoded graphs.
        d (int): Local dimension.
        cols (int): Number of columns in the subplot grid (default: 3).
    """

    num_graphs = len(encoded_graphs)
    rows = (num_graphs + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 1, rows * 1)) 
    axes = axes.flatten() if num_graphs > 1 else [axes]

    for i, encoded_graph in enumerate(encoded_graphs):
        ax = axes[i]
        n=encoded_graph[1][0]
        adjacency_matrix = bitpack_decode(encoded_graph[0], n, d)
        
        subG = nx.from_numpy_array(adjacency_matrix)
        
        edges = [(u, v) for u, v in subG.edges()]

        edge_colors = ['black' if adjacency_matrix[u, v] == 1 else 'red' for u, v in edges]
        
        node_order, min_cross = best_permutation(subG)
        pos = circular_positions(node_order) #Circular layout with the optimal node ordering
        nx.draw_networkx_edges(subG, pos, ax=ax, edge_color=edge_colors, alpha=1, width=2)
        nx.draw_networkx_nodes(subG, pos, ax=ax, node_size=40, node_color='blue', edgecolors='black', linewidths=1)
        
        ax.set_title(f"No. {i+1}")
        ax.axis("off")
        ax.set_aspect('equal')
    
    legend_elements = [
        mlines.Line2D([], [], color='black', lw=1.3, label='Edge Weight = 1'),
        mlines.Line2D([], [], color='red', lw=1.3, label='Edge Weight = 2')
    ]
    fig.legend(handles=legend_elements, loc='lower right', fontsize=13, frameon=True,bbox_to_anchor=(0.98, 0.015))
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
            if len(set(edges[i]) & set(edges[j])) == 0:
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
    adjacency_matrix = bitpack_decode(encoded_graph, n, d)
    G = nx.Graph()

    for i in range(n):
        G.add_node(i, label=f'Node {i}')
        for j in range(i + 1, n):
            if adjacency_matrix[i, j] > 0:
                G.add_edge(i, j, weight=adjacency_matrix[i, j])
    
    pos = nx.spring_layout(G)
    edge_labels = {(u, v): f"{d['weight']}" for u, v, d in G.edges(data=True)}
    nx.draw(G, pos, with_labels=True, node_color="lightblue", node_size=500, font_weight="bold")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    
    plt.title("Graph")
    plt.show()


def total_weight(graph):
    """Calculate the total weight of a graph. Here only for d=3."""
    bit_length=int(graph).bit_length()
    pair_sum = []
    for i in range(0, bit_length, 2):
        pair = (graph >> i) & 0b11 
        pair_sum.append(pair)  
    return sum(pair_sum)

def bitpack_encode(matrix, d):
    "Create bitpacked encoding of the adjacency matrix."
    n = matrix.shape[0]
    bit_length = (d-1).bit_length()  # Number of bits per weight
    packed = 0
    shift = 0
    
    # Iterate through the upper triangular part (excluding diagonal)
    for i in range(n):
        for j in range(i + 1, n):
            weight = matrix[i, j]
            packed |= (weight << shift)  # Shift and add weight to packed
            shift += bit_length
    return int(packed)

def bitpack_decode(packed, n, d):
    """Decode the bitpacked integer into an adjacency matrix."""
    bit_length = (d-1).bit_length()
    matrix = np.zeros((n, n), dtype=int)
    shift = 0
    
    # Decode weights from packed integer
    for i in range(n):
        for j in range(i + 1, n):
            weight = (packed >> shift) & ((1 << bit_length) - 1)  # Extract bits
            matrix[i, j] = weight
            matrix[j, i] = weight
            shift += bit_length
    
    return matrix


def local_scaling(G, v, k, d):
    """
    Perform local scaling.

    Args:
        G(np.array): Adjaceny matrix of the graph
        v(int): Vertex for the scaling
        k(int): Scaling factor
        d(int): Local dimension
    """

    n = G.shape[0]
    for u in range(n):
        if u != v:
            G[v, u] = (k * G[v, u]) % d
            G[u, v] = G[v, u]


def local_complementation(G, v, k, d):
    """
    Perform local complementation.

    Args:
        G(np.array): Adjaceny matrix of the graph
        v(int): Vertex for the complementation
        k(int): Scaling factor
        d(int): Local dimension
    """

    n = G.shape[0]
    for u in range(n):
        for w in range(n):
            if u!=w and u != v and w != v and G[v, u] != 0 and G[v, w] != 0:
                G[u, w] = (G[u, w] + k * G[v, u] * G[v, w]) % d
