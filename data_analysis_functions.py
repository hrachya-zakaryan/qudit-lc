"""
Script to perform various graph analysis and operations using NetworkX.

This script includes a collection of functions for graph-related tasks such as computing the maximum degree,
chromatic number, determining if a graph is a tree, calculating distance matrices, and other graph properties.
It also includes functions for decoding bit-encoded graph representations, plotting distance matrices, and checking 
planarity or self-loops. One should keep in mind that some of the functions are compuatationally heavy, for this reason
approximation techniques were applied. In any case, they are useful to be implemented.

Dependencies:
- networkx
- numpy
- matplotlib
- auxilary_unparallelized (for custom operations)

Functions:
- max_degree: Compute the maximum degree in a graph, removing self-loops.
- chromatic_number: Compute the chromatic number of a graph using a greedy coloring algorithm.
- is_tree: Check if a graph is a tree.
- distance_matrix: Compute the shortest path distance matrix of a graph.
- bfs_distance_matrix: Compute the shortest path distance matrix using BFS for unweighted graphs.
- avg_of_distance_matrix: Compute the average shortest path length of a graph.
- plot_distance_matrix: Plot the shortest path distance matrix of a graph.
- average_of_distance_matrix: Compute the average shortest path distance of a graph.
- max_in_distance_matrix: Compute the maximum shortest path distance in a graph.
- is_planar_graph: Check if a graph is planar.
- count_self_loops: Count the number of self-loops in a graph.
- decode_to_net_G: Decode a bit-encoded graph into a NetworkX MultiGraph.
- find_red_nodes: Find the minimum number of red or blue nodes in a bipartite graph.
- chromatic_number_and_color_counts: Compute the chromatic number and color counts for a graph.
- min_chromatic_number_in_OG: Compute the minimum chromatic number in a graph orbit.

Notes:
- The auxiliary functions such as `local_complementation`, `local_scaling`, `bitpack_decode`, and `bitpack_encode` are imported 
  from an external module, which should be accessible for full functionality.
"""

import networkx as nx
import numpy as np
from collections import deque, Counter
import matplotlib.pyplot as plt
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph

# Calculate the maximum degree in a graph
def max_degree(g):
    """
    Compute the maximum degree of a graph after removing self-loops.

    Parameters:
    g (networkx.Graph): A NetworkX graph object. The graph can be directed or undirected
                        and may contain self-loops.

    Returns:
    int: The maximum degree among all nodes in the graph after self-loops have been removed.
    """
    g.remove_edges_from(nx.selfloop_edges(g))
    # Compute the maximum degree 
    delta = max(dict(g.degree()).values())
    return delta



# Chromatix number of graph
def chromatic_number(G):
    """
    Compute the chromatic number of a graph G.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        int: The chromatic number of the graph.
    """
    coloring = nx.coloring.greedy_color(G, strategy="largest_first")
    return max(coloring.values()) + 1  # The highest color index + 1


# Determine if a graph is a tree or not
def is_tree(G):
    """
    Check if a given graph G is a tree.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        bool: True if the graph is a tree, False otherwise.
    """
    return nx.is_connected(G) and G.number_of_edges() == G.number_of_nodes() - 1


# Distance matrix of graph. A disclamer is that for very large graphs like OG for n>6 it is possible that this calculation take very long time.
def distance_matrix(G):
    """
    Compute the shortest path distance matrix of a graph G.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        np.ndarray: Distance matrix where entry (i, j) is the shortest path distance between nodes i and j.
    """
    return nx.floyd_warshall_numpy(G)

# We did not used this function at the end, but it could be helpful to keep for the interested user 
def bfs_distance_matrix(G):
    """
    Compute the shortest path distance matrix of an unweighted graph G using BFS.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        np.ndarray: Distance matrix where entry (i, j) is the shortest path distance between nodes i and j.
    """
    nodes = list(G.nodes())
    n = len(nodes)
    index_map = {node: i for i, node in enumerate(nodes)}  # Map nodes to indices
    dist_matrix = np.full((n, n), np.inf)  # Unreachable nodes set to infinity
    np.fill_diagonal(dist_matrix, 0)       # Distance to self is always zero

    
    # BFS function
    def bfs(start):
        dist = {node: np.inf for node in nodes}
        dist[start] = 0
        queue = deque([start])

        while queue:
            u = queue.popleft()
            for v in G.neighbors(u):
                if dist[v] == np.inf:  # If not visited
                    dist[v] = dist[u] + 1
                    queue.append(v)
        
        return [dist[node] for node in nodes]

    # Compute BFS for each node
    for i, node in enumerate(nodes):
        dist_matrix[i] = bfs(node)

    return dist_matrix

def avg_of_distance_matrix(graph):
    return nx.average_shortest_path_length(graph)



# Plot the distance matrix. One should keep in mind that this is even slower than then calculation of the distance matrix.
def plot_distance_matrix(G):
    """
    Compute and plot the shortest path distance matrix of a graph G using matplotlib.
    The matrix contains only integer distances.

    Parameters:
        G (networkx.Graph): Input graph.
    """
    # Get the list of nodes
    nodes = list(G.nodes())
    n = len(nodes)
    
    # Initialize an empty distance matrix
    dist_matrix = np.zeros((n, n), dtype=int)
    
    # Compute the shortest path lengths
    for i, node_i in enumerate(nodes):
        lengths = nx.shortest_path_length(G, source=node_i)
        for j, node_j in enumerate(nodes):
            if node_j in lengths:
                dist_matrix[i, j] = int(lengths[node_j])  # Ensure integers
            else:
                dist_matrix[i, j] = np.inf  # Use infinity for unreachable nodes

    # Plot the matrix using imshow
    plt.figure(figsize=(8, 6))
    plt.imshow(dist_matrix, cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Shortest Path Distance")
    plt.title("Distance Matrix")
    plt.xlabel("Nodes")
    plt.ylabel("Nodes")
    plt.xticks(ticks=range(n), labels=nodes, rotation=90)
    plt.yticks(ticks=range(n), labels=nodes)

    # Show the plot
    plt.show()

#If one has determined the distance matrix, one can use the following two functions. In our work, we resorted to heuristic methods in the corresponding files. 
#However, it is still useful to keep them because they were used to produce some figures in our work.
def average_of_distance_matrix(G):
    """
    Compute the average shortest path distance of a graph G, including the diagonal elements.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        float: The average shortest path distance.
    """
    dist_matrix = distance_matrix(G)  
    return np.mean(dist_matrix)  # Include all elements, including diagonal

def max_in_distance_matrix(G):
    """
    Compute the maximum shortest path distance in a graph G, including the diagonal elements.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        int: The maximum shortest path distance.
    """
    dist_matrix = distance_matrix(G)  
    return int(np.max(dist_matrix))  # Ensure the result is an integer

# Check if a graph is planar
def is_planar_graph(G):
    """
    Check if a graph is planar.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        bool: True if the graph is planar, False otherwise.
    """
    is_planar, _ = nx.check_planarity(G) # the second output is a counter example, which by default is False
    return is_planar

#Count self-loops. In 1910.03969 (Adcock et al.) they are answering if the OG graph has loops or not.
#We can extend the idea and find how many loops we have. 
def count_self_loops(G):
    """
    Count the number of self-loops in a graph.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        int: Number of self-loops in the graph.
    """
    return nx.number_of_selfloops(G)

# From encoded graphs  to networkx.MultiGraphs
def decode_to_net_G(bit_encoded, n, d):
    """
    Decodes a bit-encoded representation into a NetworkX graph.

    This function decodes a bit-packed representation of a graph into a 
    NumPy array, and then constructs a multi-graph using NetworkX. The graph 
    is represented as a MultiGraph, which allows multiple edges between nodes.

    Parameters:
    bit_encoded (str or list): The bit-encoded representation of the graph, 
                                which is to be decoded.
    n (int): The number of nodes in the graph.
    d (int): The dimension of the graph or the encoding depth.

    Returns:
    networkx.MultiGraph: A NetworkX MultiGraph object constructed from the 
                          decoded bit-encoded data, with multiple edges allowed.

    Example:
    g = decode_to_net_G(encoded_bits, 5, 3)
    """
    temp = bitpack_decode(bit_encoded, n, d)
    g = nx.from_numpy_array(temp, parallel_edges=True, create_using=nx.MultiGraph())  # Explicitly use MultiGraph
    return g


def find_red_nodes(bit_encoded, n, d):
    """
    Finds the minimum number of red or blue nodes in a bipartite graph 
    represented by a bit-encoded string.

    This function decodes a bit-encoded representation into a NetworkX graph 
    and then identifies the bipartite sets (red and blue nodes). It returns 
    the smaller of the two sets' sizes, which represents the minimum number 
    of nodes in either set.

    Parameters:
    bit_encoded (str or list): The bit-encoded representation of the graph 
                                to be decoded.
    n (int): The number of nodes in the graph.
    d (int): The dimension of the graph or the encoding depth.

    Returns:
    int: The smaller of the two sets' sizes (the red and blue sets in the 
         bipartite graph).

    Example:
    min_red_nodes = find_red_nodes(encoded_bits, 5, 3)
    """
    temp = decode_to_net_G(bit_encoded, n, d)
    red_nodes, blue_nodes = nx.algorithms.bipartite.sets(temp)
    n_red_nodes = len(red_nodes)
    n_blue_nodes = len(blue_nodes)
    return min(n_red_nodes, n_blue_nodes)


def chromatic_number_and_color_counts(graph):
    """
    Determine the chromatic number of a graph and the count of each color.

    Parameters:
    -----------
    graph : networkx.Graph
        The input graph.

    Returns:
    --------
    tuple
        A tuple containing:
        - chromatic_number (int): The minimum number of colors required to color the graph.
        - color_counts (dict): A dictionary where keys are color indices and values are the counts of nodes with that color.
    """
    # Color the graph using a greedy algorithm
    color_mapping = nx.coloring.greedy_color(graph, strategy="largest_first")

    # Count the occurrences of each color
    color_counts = Counter(color_mapping.values())

    # The chromatic number is the highest color index + 1 (because colors are 0-indexed)
    chromatic_number = max(color_mapping.values()) + 1

    return chromatic_number, dict(color_counts)



def min_chromatic_number_in_OG(orbit_graph, n, d):
    """
    Computes the minimum chromatic number for a set of graphs in an orbit.

    This function decodes a list of encoded graphs in the orbit into NetworkX 
    graphs, computes the chromatic number for each graph, and returns the 
    smallest chromatic number among them.

    Parameters:
    orbit_graph (networkx.Graph): A NetworkX graph where nodes represent 
                                  encoded graphs in the orbit.
    n (int): The number of nodes in each graph.
    d (int): The dimension of the graphs or the encoding depth.

    Returns:
    int: The minimum chromatic number among the decoded graphs in the orbit.

    Example:
    min_chromatic = min_chromatic_number_in_OG(orbit, 5, 3)
    """
    encoded_graphs_in_orbit = list(orbit_graph.nodes())
    graphs_of_orbit_in_G_form = []
    for i in range(len(encoded_graphs_in_orbit)):
        graphs_of_orbit_in_G_form.append(decode_to_net_G(encoded_graphs_in_orbit[i], n, d))
    
    list_of_chromatic_number = []
    for i in range(len(graphs_of_orbit_in_G_form)):
        list_of_chromatic_number.append(chromatic_number(graphs_of_orbit_in_G_form[i]))
    
    return min(list_of_chromatic_number)

    