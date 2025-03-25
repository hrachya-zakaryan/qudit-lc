import networkx as nx
from collections import deque
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph, draw_graphs
import numpy as np
import time
import multiprocessing
import os
from matplotlib import pyplot as plt
import subprocess
import matplotlib.lines as mlines

def total_weight(graph):
    bit_length=int(graph).bit_length()
    pair_sum = []
    for i in range(0, bit_length, 2):
        pair = (graph >> i) & 0b11 
        pair_sum.append(pair)  

    return sum(pair_sum)


def call_generate_graphs(n,d, encoded_value,path,index,ts):
    #print("start:",index,":",time.time()-ts)
    result = subprocess.run(
        [path, str(n), str(d), str(encoded_value)], 
        capture_output=True, text=True
    )
    print("end", index,":",time.time()-ts)
    output_dir = os.path.join(os.getcwd(),f"d{d}_c",str(n))
    graphs = {int(line) for line in result.stdout.splitlines()}
    output_file = os.path.join(output_dir, f"{index}.txt")
    with open(output_file, "w") as f:
        for graph in graphs:
            f.write(str(graph) + "\n")

def generate_graphs_c(n,d,path):
    filename=f"d3n{n}.txt"
    with open(filename, "r") as file:
        graph6_lines = [line.strip() for line in file if line.strip()]
    iso_graphs = [bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines]
    ts=time.time()

    output_dir = os.path.join(os.getcwd(),f"d{d}_c",str(n))
    os.makedirs(output_dir, exist_ok=True)
    num_workers = multiprocessing.cpu_count()
    with multiprocessing.Pool(num_workers) as pool:
        index=0
        while iso_graphs:
            batch_size = min(len(iso_graphs), num_workers)
            batch = [iso_graphs.pop(0) for _ in range(batch_size)]
            print(f"{len(iso_graphs)}:{time.time()-ts}")
            #weighted_graphs = pool.starmap(call_generate_graphs, [(n, d, g, path,index*num_workers+batch.index(g),ts) for g in batch])
            pool.starmap(call_generate_graphs, [(n, d, g, path,index*num_workers+batch.index(g),ts) for g in batch])
            # for i in range(batch_size):
            #     graphs=weighted_graphs[i]
            #     output_file = os.path.join(output_dir, f"{index*num_workers+i}.txt")
            #     with open(output_file, "w") as f:
            #         for graph in graphs:
            #             f.write(str(graph) + "\n")  # Write each graph on a new line
            #     #print(f"{index*num_workers+i} end:{time.time()-ts}")
            index+=1
    


def call_complementation_layer(n, d, encoded_value,path):
    """ Calls the C program and returns results as a set. """
    result = subprocess.run(
        [path, str(n), str(d), str(encoded_value)], 
        capture_output=True, text=True
    )
    return {int(line) for line in result.stdout.splitlines()}


def orbit_atlas_c(path,n,d):
    graphs=set()
    ts=time.time()
    directory=os.path.join(os.getcwd(), f"d{d}_c", str(n))
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        
        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                for line in file:
                    graph = int(line.strip()) 
                    graphs.add(graph)
    G=nx.Graph()
    G.add_nodes_from(graphs)

    num_workers = multiprocessing.cpu_count()
    with multiprocessing.Pool(num_workers) as pool:
        while graphs:
            batch_size = min(len(graphs), num_workers)
            batch = [graphs.pop() for _ in range(batch_size)]
            print(f"{len(graphs)}:{time.time()-ts}")
            complements = pool.starmap(call_complementation_layer, [(n, d, g, path) for g in batch])
            for i in range(batch_size):
                for complement in complements[i]:
                    G.add_edge(batch[i],complement)
    return G


def separate_orbits(G):
    components = list(nx.connected_components(G))
    subgraphs = [G.subgraph(nodes).copy() for nodes in components]
    return subgraphs

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

    # Set edge colors based on the weight
    edge_colors = ['green' if adjacency_matrix[u, v] == 1 else 'red' if adjacency_matrix[u, v] == 2 else 'blue' if adjacency_matrix[u, v] == 3 else 'black' for u, v in edges]

    # Draw subgraph edges with the appropriate color
    nx.draw_networkx_edges(subG, pos, ax=ax, edge_color=edge_colors, alpha=1, width=1.3)

    # Draw subgraph nodes (light blue fill, black outline)
    nx.draw_networkx_nodes(subG, pos, ax=ax, node_size=23, node_color='blue', edgecolors='black', linewidths=0.6)

    # Draw edge weights (small font)
    #nx.draw_networkx_edge_labels(subG, pos, edge_labels=edge_weights, ax=ax, font_size=4, font_color='black')

def plot_graph(G, n, d):
    """
    Plot a graph G using NetworkX with each node as a circle containing a small graph.
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    pos = nx.spring_layout(G, seed=37, k=1/np.sqrt(len(G.nodes)))  # Layout for main graph

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

    # Create a legend for subgraph edge colors
    legend_elements = [
        mlines.Line2D([], [], color='green', lw=1.3, label='Edge Weight = 1'),
        mlines.Line2D([], [], color='red', lw=1.3, label='Edge Weight = 2'),
        mlines.Line2D([], [], color='blue', lw=1.3, label='Edge Weight = 3'),
        mlines.Line2D([], [], color='black', lw=1.3, label='Edge Weight = 4')
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10, frameon=True)
    plt.show()
n=5
d=7

#generate_graphs_c(n,d,"c/generate_graphs")
g=orbit_atlas_c("c/complement",n,d)
# if __name__ == "__main__":
#     multiprocessing.freeze_support()
#     g=orbit_atlas_c(n,d)
# g=nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_6", data=False)
# print("Vertices:",g.number_of_nodes())
# print("Edges:", g.number_of_edges())
# encoded_value=bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes("EU~w".encode()),dtype=int),d)
# call_generate_graphs(n,d, encoded_value,"c/generate_graphs",96,time.time())
orbits=separate_orbits(g)

directory = f"orbits_d{d}_n{n}_separated"
os.makedirs(directory, exist_ok=True)
for i in range(len(orbits)):
    nx.write_edgelist(orbits[i],f"{directory}/orbit_{i}", data=False)

# # print(time.time()-ts)
# for o in orbits:
#   plot_graph(o,n,d)
#print(len(orbits))

#orbit_atlas_c(n,d)
