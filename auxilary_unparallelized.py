import numpy as np
from itertools import product, combinations
import networkx as nx
import matplotlib.pyplot as plt


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
    bit_length = d.bit_length()  # Number of bits per weight
    packed = 0
    shift = 0
    
    # Iterate through the upper triangular part (excluding diagonal)
    for i in range(n):
        for j in range(i + 1, n):
            weight = matrix[i, j]
            packed |= (weight << shift)  # Shift and add weight to packed
            shift += bit_length  # Update shift for the next weight
    
    return packed

def bitpack_decode(packed, n, d):
    bit_length = d.bit_length()
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

def create_graphs(n: int,d: int):
    num_entries = n * (n - 1) // 2
    counter=0
    graphs=[]
    for combination in product(range(d), repeat=num_entries):
            # Reconstruct the symmetric matrix
        matrix = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                matrix[i, j] = combination[idx]
                matrix[j, i] = combination[idx]
                idx += 1
        if nx.is_connected(nx.from_numpy_array(matrix)):
            graphs.append(bitpack_encode(matrix,d))
            counter+=1

    return graphs

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



def save_to_file(graphs, filename="encoded_graphs.txt"):
    """
    Save the encoded graphs to a plain text file, one per line.

    Parameters:
        graphs (list[int]): List of bit-encoded graphs.
        filename (str): Name of the file to save the encoded graphs.
    """
    with open(filename, "w") as f:
        for encoded_graph in graphs:
            f.write(f"{encoded_graph}\n")
    
    print(f"Saved {len(graphs)} encoded graphs to {filename}")
    
def load_from_file(filename="encoded_graphs.txt"):
    """
    Load encoded graphs from a plain text file.

    Parameters:
        filename (str): Name of the file to load the encoded graphs from.

    Returns:
        list[int]: List of bit-encoded graphs.
    """
    with open(filename, "r") as f:
        return [int(line.strip()) for line in f.readlines()]

def save_orbits_to_file(orbits, filename="orbits.txt"):
    """
    Save the orbits (equivalence classes) to a text file.

    Parameters:
        orbits (list[set[int]]): List of orbits, each as a set of encoded graphs.
        filename (str): Name of the file to save the orbits.
    """
    with open(filename, "w") as f:
        for orbit in orbits:
            line = ",".join(map(str, orbit))  # Convert each orbit to a comma-separated line
            f.write(f"{line}\n")
    
    print(f"Saved {len(orbits)} orbits to {filename}")

def save_representatives_to_file(orbits, filename="representatives.txt"):
    """
    Save the representatives (first element of each orbit) to a text file.

    Parameters:
        orbits (list[set[int]]): List of orbits, each as a set of encoded graphs.
        filename (str): Name of the file to save the representatives.
    """
    with open(filename, "w") as f:
        for orbit in orbits:
            representative = next(iter(orbit))  # Get the first element of the set
            f.write(f"{representative}\n")
    
    print(f"Saved {len(orbits)} representatives to {filename}")

def load_orbits_from_file(filename="orbits.txt"):
    """
    Load the orbits (equivalence classes) from a text file.

    Parameters:
        filename (str): Name of the file containing the orbits.

    Returns:
        list[set[int]]: A list of orbits, where each orbit is represented as a set of encoded graphs.
    """
    orbits = []
    with open(filename, "r") as f:
        for line in f:
            orbit = set(map(int, line.strip().split(',')))  # Parse each line into a set of integers
            orbits.append(orbit)
    print(f"Loaded {len(orbits)} orbits from {filename}")
    return orbits



def generate_connected_graphs(n):
    """Generate all connected graphs with n vertices, removing isomorphic duplicates."""
    # Generate all possible graphs
    nodes = range(n)
    all_edges = list(combinations(nodes, 2))  # All possible edges for n vertices
    unique_graphs = []
    seen_canonicals = set()

    # Iterate over all subsets of edges to form graphs
    for k in range(len(all_edges) + 1):
        for edges in combinations(all_edges, k):
            G = nx.Graph()
            G.add_nodes_from(nodes)
            G.add_edges_from(edges)

            # Check if the graph is connected
            if nx.is_connected(G):
                # Use a canonical form to detect isomorphism
                canonical = nx.canonical_graph(G)
                canonical_string = nx.to_graph6_bytes(canonical).decode()

                if canonical_string not in seen_canonicals:
                    seen_canonicals.add(canonical_string)
                    unique_graphs.append(G)

    return unique_graphs