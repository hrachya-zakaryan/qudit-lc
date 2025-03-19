import networkx as nx
import numpy as np
from collections import deque, Counter
import matplotlib.pyplot as plt
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph

# def plot_graph(G):
#     """
#     Plot a given graph G using NetworkX and Matplotlib with improved aesthetics.
    
#     Parameters:
#         G (networkx.Graph): The graph to be plotted.
#     """
#     # Choose a layout for better structure
#     pos = nx.spring_layout(G, seed=42, k=1/np.sqrt(len(G.nodes)))  # k controls spacing between nodes
    
#     # Get node sizes based on degree
#     degrees = dict(G.degree())
#     node_sizes = [100 + degrees[node] * 50 for node in G.nodes()]
    
#     # Draw nodes
#     nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=node_sizes, edgecolors='black')
    
#     # Draw edges with transparency
#     nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.5)
    
#     # Draw labels
#     nx.draw_networkx_labels(G, pos, font_size=10, font_family='sans-serif')
    
#     # Add a title with graph info
#     plt.title(f"Graph with {len(G.nodes)} Nodes and {len(G.edges)} Edges", fontsize=14)
    
#     # Increase figure size and resolution for better visibility
#     plt.gcf().set_size_inches(12, 10)
#     plt.show()

# chromatix number of graph
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


#determine if a graph is a tree or not
def is_tree(G):
    """
    Check if a given graph G is a tree.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        bool: True if the graph is a tree, False otherwise.
    """
    return nx.is_connected(G) and G.number_of_edges() == G.number_of_nodes() - 1


#distance matrix of graph
def distance_matrix(G):
    """
    Compute the shortest path distance matrix of a graph G.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        np.ndarray: Distance matrix where entry (i, j) is the shortest path distance between nodes i and j.
    """
    return nx.floyd_warshall_numpy(G)


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



#plot it
def plot_distance_matrix(G):
    """
    Compute and plot the shortest path distance matrix of a graph G using matplotlib.
    The matrix contains only integer distances.

    Parameters:
        G (networkx.Graph): Input graph.
    """
    #TODO: Use the above function
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

#check if a graph is planar
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

#count self loops. In Adcocks paper they are just answering if the OG graph has loops or not. We can extend the idea and find how many loops we have. We can get the same information as before is this function returns 0
def count_self_loops(G):
    """
    Count the number of self-loops in a graph.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        int: Number of self-loops in the graph.
    """
    return nx.number_of_selfloops(G)


#count the number of circles in the graph, can be deternmined in polynomial time O(V+E), see: https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.cycles.cycle_basis.html
def count_cycles(G):
    """
    Count the number of independent cycles (cycle basis) in an undirected graph.

    Parameters:
        G (networkx.Graph): Input undirected graph.

    Returns:
        int: The number of independent cycles in the graph.
    """
    return len(nx.cycle_basis(G))

#has Eulerian circle, see https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.euler.is_eulerian.html
def has_eulerian_cycle(G):
    """
    Check if a graph has an Eulerian cycle.

    Parameters:
        G (networkx.Graph): Input graph.

    Returns:
        bool: True if the graph has an Eulerian cycle, False otherwise.
    """
    return nx.is_eulerian(G)


def find_orbit_graph(start_graph, n, d):
    """
    Constructs an undirected graph where each node represents a graph state,
    and an edge exists if one graph can be reached from another using
    local scaling or local complementation.
    
    Parameters:
    - start_graph: The starting graph state (bitpack encoded)
    - n: Number of vertices in the graph
    - d: Local dimension 
    
    Returns:
    - G: NetworkX undirected graph representing the orbit
    - min_graph: The minimum graph state in the orbit
    """
    G = nx.Graph()  
    queue = deque([start_graph])  
    visited = set()  # Set to track visited graphs
    visited.add(start_graph)
    G.add_node(start_graph)
    min_graph = start_graph
    
    scaling_factors = range(2, d)  
    complementing_factors = range(1, d)  

    while queue:
        current = queue.popleft()
        current_matrix = bitpack_decode(current, n, d)
        
        # Update minimum graph
        if current.bit_count() <= min_graph.bit_count() and current < min_graph:
            min_graph = current
        
        # Try all local transformations
        for v in range(n):
            # Local scaling
            for k in scaling_factors:
                scaled_matrix = current_matrix.copy()
                local_scaling(scaled_matrix, v, k, d)
                encoded_scaled = bitpack_encode(scaled_matrix, d)
                
                # Add the scaled state as a node and create an edge
                if encoded_scaled not in visited:
                    visited.add(encoded_scaled)
                    queue.append(encoded_scaled)
                    G.add_node(encoded_scaled)
                
                G.add_edge(current, encoded_scaled)  # Undirected edge
                
                # Add a loop if it leads back to itself
                if encoded_scaled == current:
                    G.add_edge(current, current)
            
            # Local complementation
            for k in complementing_factors:
                complemented_matrix = current_matrix.copy()
                local_complementation(complemented_matrix, v, k, d)
                encoded_complemented = bitpack_encode(complemented_matrix, d)
                
                # Add the complemented state as a node and create an edge
                if encoded_complemented not in visited:
                    visited.add(encoded_complemented)
                    queue.append(encoded_complemented)
                    G.add_node(encoded_complemented)
                
                G.add_edge(current, encoded_complemented)  # Undirected edge
                
                # Add a loop if it leads back to itself
                if encoded_complemented == current:
                    G.add_edge(current, current)
    
    return G, min_graph


#has Hamiltonian circle this is an NP-complete problem, we are not going to implement it, even the paper they are not including it into the correlation analysis.
#Automorphism group is NP-hard to compute, we are not going to implement it, even the paper they are not including it into the correlation analysis.





#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Here we implement the functions to calculate the properties of the graphs that they are element of the OG
#Schmidt measure. It is interesting because in LC-orbits paper, they are not calculating this by themselves, they just use the previous results from the optimal prep paper. 
#The optimal paper itself seems just to provide some bounds on it and not really calculating the final result. TODO: Check the web for other implentations 


#encoded to networkx.G
def decode_to_net_G(bit_encoded,n,d):
    temp = bitpack_decode(bit_encoded,n,d)
    g = nx.from_numpy_array(temp)
    return g
    
def find_red_nodes(bit_encoded,n,d):
    temp = decode_to_net_G(bit_encoded,n,d)
    red_nodes, blue_nodes = nx.algorithms.bipartite.sets(temp)
    n_red_nodes = len(red_nodes)
    n_blue_nodes = len(blue_nodes)
    return min(n_red_nodes, n_blue_nodes)


def schmidt_measure(orbit, n, d, print_nodes=False):
    """
    Calculate the Schmidt measure for a given orbit of encoded graphs.

    This function identifies all two-colorable graphs in the given orbit,
    determines the minimum number of red nodes for each graph, and returns
    the minimum value among them as the Schmidt measure.

    Parameters:
    -----------
    orbit : list of int
        A list of integers representing encoded graphs.
    n : int
        The number of nodes in each graph.
    d : int
        The local dimension of the graph state.
    print_nodes : bool, optional
        If True, prints the encoding and the number of red nodes for each 
        two-colorable graph. Default is False.

    Returns:
    --------
    int
        The minimum number of red nodes among the two-colorable graphs in the orbit,
        representing the normalized Schmidt measure.

    Notes:
    ------
    This function assumes that `decode_to_net_G` and `find_red_nodes` are defined elsewhere
    and that the graphs are encoded such that their chromatic number can be determined.

    Example:
    --------
    >>> orbit = [21, 42, 1365, 2709]
    >>> schmidt_measure(orbit, 4, 3, print_nodes=True)
    encoding: 21, has 2 red nodes
    encoding: 42, has 1 red node
    1
    """

    encoded_two_colorable = []
    for i in range(len(orbit)):
        # g = decode_to_net_G(orbit[i],n,d)
        g = decode_to_net_G(orbit[i],n,d)
        color = chromatic_number(g)
        if color == 2:
            encoded_two_colorable.append(orbit[i])

    array_of_minimum_number_of_red = []        
    # print_nodes = False
    for i  in range(len(encoded_two_colorable)):
        temp = find_red_nodes(encoded_two_colorable[i],n,d)
        array_of_minimum_number_of_red.append(temp)
        if print_nodes == True:    
            if temp == 1:
                print(f"encoding: {encoded_two_colorable[i]}, has {temp} red node")
            else:
                print(f"encoding: {encoded_two_colorable[i]}, has {temp} red nodes")

    schmidt_measure_normalized = min(array_of_minimum_number_of_red)

    return schmidt_measure_normalized




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




def min_chromatic_number_in_OG(orbit_graph,n,d):
    encoded_graphs_in_orbit = list(orbit_graph.nodes())
    graphs_of_orbit_in_G_form = []
    for i in range(len(encoded_graphs_in_orbit)):
        graphs_of_orbit_in_G_form.append(decode_to_net_G(encoded_graphs_in_orbit[i],n,d))
    list_of_chromatic_number = []
    for i in range(len(graphs_of_orbit_in_G_form)):
        list_of_chromatic_number.append(chromatic_number(graphs_of_orbit_in_G_form[i]))
    
    return min(list_of_chromatic_number)
    