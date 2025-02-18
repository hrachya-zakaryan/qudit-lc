from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph, draw_graphs
from isomorphism_method import generate_permuted, orbit_search_isomorphic_from_file, find_orbit_isomorphic, orbit_search
import networkx as nx
from collections import deque
from matplotlib import pyplot as plt
import numpy as np
import os

def create_orbital_graph(start_graph, n, d, G):   
    
    
    # Define all possible scaling and complementation factors
    scaling_factors = range(2, d)  # Local scaling factors (mod d)
    complementing_factors = range(1, d)
    

    current = start_graph
    current_matrix = bitpack_decode(current, n, d)
    
    for v in range(n):
        # Local scaling
        for k in scaling_factors:
            scaled_matrix = current_matrix.copy()
            local_scaling(scaled_matrix, v, k, d)
            encoded_scaled = bitpack_encode(scaled_matrix, d)
            perms=set(generate_permuted(bitpack_decode(encoded_scaled,n,d),d))
            inter_nodes=perms.intersection(set(G.nodes))
            encoded_scaled=inter_nodes.pop()
            G.add_edge(encoded_scaled,current)
        
        # Local complementation
        for k in complementing_factors:
            complemented_matrix = current_matrix.copy()
            local_complementation(complemented_matrix, v, k, d)
            encoded_complemented = bitpack_encode(complemented_matrix, d)
            perms=set(generate_permuted(bitpack_decode(encoded_complemented,n,d),d))
            inter_nodes=perms.intersection(set(G.nodes))
            encoded_complemented=inter_nodes.pop()
            G.add_edge(encoded_complemented,current)
            



def orbit_search_isomorphic_from_file(n,d):
    # Read the graph6 file
    
    graphs = set()

    # List all files in the directory
    directory=os.path.join(os.getcwd(), str(n))
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Ensure we only read files (skip directories)
        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                for line in file:
                    graph = int(line.strip())  # Convert graph from string to integer
                    graphs.add(graph)
    orbits_reps=[]
    orbits=[]
    while graphs:
        current_graph = graphs.pop()
        temp_orbit, current_graph = find_orbit_isomorphic(graphs,current_graph,n,d)
        orbits_reps.append(current_graph)
        orbits.append(temp_orbit)
        
    return orbits, orbits_reps



def circular_subgraph_layout(n, center, radius):
    """ Arrange n points in a circular layout inside a given center and radius. """
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return {i: center + radius * np.array([np.cos(angle), np.sin(angle)]) for i, angle in enumerate(angles)}

def draw_subgraph_inside_circle(ax, center, radius, adjacency_matrix):
    """ Draws a small subgraph inside a circular boundary. """
    subG = nx.from_numpy_array(adjacency_matrix)  # Convert adjacency matrix to NetworkX graph
    pos = circular_subgraph_layout(len(subG.nodes), np.array(center), radius * 0.7)  # Keep subgraph size, shrink nodes

    # Get weighted edges
    edges = [(u, v) for u, v in subG.edges()]
    edge_weights = {(u, v): int(adjacency_matrix[u, v]) for u, v in edges if adjacency_matrix[u, v] > 0}

    # Draw subgraph edges (light blue)
    nx.draw_networkx_edges(subG, pos, ax=ax, edge_color='blue', alpha=0.8, width=0.5)

    # Draw subgraph nodes (light blue fill, black outline)
    nx.draw_networkx_nodes(subG, pos, ax=ax, node_size=20, node_color='blue', edgecolors='black', linewidths=0.6)

    # Draw edge weights (small font)
    nx.draw_networkx_edge_labels(subG, pos, edge_labels=edge_weights, ax=ax, font_size=4, font_color='black')

def plot_graph(G, n, d):
    """
    Plot a graph G using NetworkX with each node as a circle containing a small graph.
    Parameters:
        G (networkx.Graph): The main graph to be plotted.
        bitpack_decode (function): Function to decode node labels into adjacency matrices.
        n (int): Number of vertices in each subgraph.
        d (int): Parameter for decoding function.
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    pos = nx.spring_layout(G, seed=42, k=1/np.sqrt(len(G.nodes)))  # Layout for main graph

    # Draw main graph edges (light blue)
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='blue', alpha=0.7, width=0.8)

    # Draw main graph nodes as circles
    for node in G.nodes():
        x, y = pos[node]  # Get node position
        
        # Draw a large circle for the node (white fill, black border)
        circle = plt.Circle((x, y), 0.1, color='white', ec='black', lw=1.2)
        ax.add_patch(circle)
        
        # Get adjacency matrix for the subgraph
        adjacency_matrix = bitpack_decode(node, n, d)
        
        # Draw the subgraph inside the node circle
        draw_subgraph_inside_circle(ax, (x, y), 0.08, adjacency_matrix)

    # Set axis limits and aspect ratio
    ax.set_xlim(min(x for x, y in pos.values()) - 0.2, max(x for x, y in pos.values()) + 0.2)
    ax.set_ylim(min(y for x, y in pos.values()) - 0.2, max(y for x, y in pos.values()) + 0.2)
    ax.set_aspect('equal')
    ax.axis('off')

    # Add a title
    plt.title(f"Graph with {len(G.nodes)} Nodes and {len(G.edges)} Edges", fontsize=14)
    plt.show()
n=3
d=3

orbits=orbit_search_isomorphic_from_file(n,d)[0]

for orbit in orbits:
    G=nx.Graph()
    for g in orbit:
        G.add_node(g)
    for g in orbit:
        create_orbital_graph(g,n,d,G)
        
    plot_graph(G,n,d)
#print(G.edges)
#draw_graphs([41,42],n,d)
