import networkx as nx
from collections import deque
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph, draw_graphs
from itertools import permutations
import numpy as np
from itertools import product
import time
from networkx.algorithms.isomorphism import GraphMatcher
import asyncio
import concurrent.futures
import multiprocessing
import os
from matplotlib import pyplot as plt
import subprocess

def generate_permuted(adj_matrix,d):
    """Generate all possible vertex permutations of an adjacency matrix."""
    G=nx.from_numpy_array(adj_matrix)
    node_permutations = list(permutations(G.nodes()))
    
    # Generate relabeled graphs
    relabeled_graphs = set()
    for perm in node_permutations:
        mapping = {old: new for old, new in zip(G.nodes(), perm)}
        new_G = nx.relabel_nodes(G, mapping)
        sorted_nodes = sorted(new_G.nodes())
        g=bitpack_encode(nx.to_numpy_array(new_G,dtype=int,nodelist=sorted_nodes,weight="weight"),d)
        
        relabeled_graphs.add(g)

    return relabeled_graphs


def generate_weighted(adj_matrix,d):
    """Generate all possible weighted adjacency matrices from a binary adjacency matrix."""
    
    n = len(adj_matrix)
    upper_indices = [(i, j) for i in range(n) for j in range(i + 1, n) if adj_matrix[i, j] == 1]

    # Generate all possible weight combinations (1 or 2) for these positions
    num_edges = len(upper_indices)
    weight_combinations = product(range(1,d), repeat=num_edges)

    weighted_matrices = []
    for weights in weight_combinations:
        # Create a copy of the original matrix
        new_matrix = adj_matrix.copy()
        # Assign weights symmetrically
        for (pos, weight) in zip(upper_indices, weights):
            i, j = pos
            new_matrix[i, j] = weight
            new_matrix[j, i] = weight
        weighted_matrices.append(bitpack_encode(new_matrix,d))

    return weighted_matrices

executor = concurrent.futures.ProcessPoolExecutor()  # Use processes instead of threads

# Wrap CPU-bound functions to run in threads
async def async_generate_weighted(adj_matrix, d):
    """Run generate_weighted in a thread."""
    loop = asyncio.get_running_loop()
    return await loop.run_in_executor(executor, generate_weighted, adj_matrix, d)

async def async_generate_permuted(adj_matrix, d):
    """Run generate_permuted in a thread."""
    loop = asyncio.get_running_loop()
    return await loop.run_in_executor(executor, generate_permuted, adj_matrix, d)


async def process_graph(ig,n,d,output_dir, index):
        temp_ts=time.time()
        results=[] 
        weighted_graphs = set(await async_generate_weighted(bitpack_decode(ig,n,d),d))
    
        while weighted_graphs:
            g=weighted_graphs.pop()
            
            perms=await async_generate_permuted(bitpack_decode(g,n,d),d)
            results.append(min(perms))
            weighted_graphs.difference_update(perms)
         # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Write to a file named by index
        output_file = os.path.join(output_dir, f"{index}.txt")
        with open(output_file, "w") as f:
            for graph in results:
                f.write(str(graph) + "\n")  # Write each graph on a new line
        print(f"{index} end:{time.time()-temp_ts}")

async def generate_graphs(filename, n, d):
    
    with open(filename, "r") as file:
        graph6_lines = [line.strip() for line in file if line.strip()]

    iso_graphs = [bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines]
    ts=time.time()

    output_dir = os.path.join(os.getcwd(), str(n))
    os.makedirs(output_dir, exist_ok=True)
    await asyncio.gather(*[process_graph(ig,n,d,output_dir,i) for i,ig in enumerate(iso_graphs)])
    print(f"end: {time.time()-ts}")

def total_weight(graph):
    bit_length=int(graph).bit_length()
    pair_sum = []
    for i in range(0, bit_length, 2):
        pair = (graph >> i) & 0b11 
        pair_sum.append(pair)  

    return sum(pair_sum)

def find_orbit(start_graph, n, d):
    queue = deque([start_graph])  # BFS queue
    visited = set()  # Set of visited graphs
    visited.add(start_graph)
    
    # Define all possible scaling and complementation factors
    scaling_factors = range(2, d)  # Local scaling factors (mod d)
    complementing_factors = range(1, d)
    

    while queue:
        current = queue.popleft()
        current_matrix = bitpack_decode(current, n, d)
       
      
        for v in range(n):
            # Local scaling
            for k in scaling_factors:
                scaled_matrix = current_matrix.copy()
                local_scaling(scaled_matrix, v, k, d)
                encoded_scaled = bitpack_encode(scaled_matrix, d)
                if encoded_scaled not in visited:
                    visited.add(encoded_scaled)
                    queue.append(encoded_scaled)
            
            # Local complementation
            for k in complementing_factors:
                complemented_matrix = current_matrix.copy()
                local_complementation(complemented_matrix, v, k, d)
                encoded_complemented = bitpack_encode(complemented_matrix, d)
                if encoded_complemented not in visited:
                    visited.add(encoded_complemented)
                    queue.append(encoded_complemented)
    #print(len(visited))
    return visited


def all_complementations(current,current_level, n, d, visited):
    """
    Process a single graph: apply local scaling and complementation.
    Returns a set of new unique graphs found in this step.
    """
    current_matrix = bitpack_decode(current, n, d)
    scaling_factors = range(2, d)
    complementing_factors = range(1, d)
    new_graphs = set()

    for v in range(n):
        # Local scaling
        for k in scaling_factors:
            scaled_matrix = current_matrix.copy()
            local_scaling(scaled_matrix, v, k, d)
            encoded_scaled = bitpack_encode(scaled_matrix, d)
            if encoded_scaled not in visited and encoded_scaled not in current_level:
                new_graphs.add(encoded_scaled)

        # Local complementation
        for k in complementing_factors:
            complemented_matrix = current_matrix.copy()
            local_complementation(complemented_matrix, v, k, d)
            encoded_complemented = bitpack_encode(complemented_matrix, d)
            if encoded_complemented not in visited and encoded_complemented not in current_level:
                new_graphs.add(encoded_complemented)

    return new_graphs



def find_orbit_isomorphic(start_graph, graphs, n, d):
    queue = deque([start_graph])  # BFS queue
    visited = set()  # Set of visited graphs
    visited.add(start_graph)
    ts=time.time()
    min_graph=start_graph
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    while queue:
        current_level = list(queue)
        queue.clear()
        print(f"{len(current_level)}:{time.time()-ts}")
        results = pool.starmap(all_complementations, [(g, current_level, n, d, visited) for g in current_level])
        temp_level=set().union(*results)
        
        while temp_level:
            g=temp_level.pop()
            visited.add(g)
            perms=generate_permuted(bitpack_decode(g,n,d),d)
            temp_level.difference_update(perms)
            graphs.difference_update(perms)
            if g not in temp_level:
                queue.append(g)
    
    return visited, min_graph

# for current in current_level:
            
#             current_matrix = bitpack_decode(current, n, d)
#             # if current.bit_count()<=min_graph.bit_count():
#             #     if total_weight(current)<=total_weight(min_graph):
#             #         if current<min_graph:
#             #             min_graph=current
            
#             for v in range(n):
#                 # Local scaling
#                 for k in scaling_factors:
#                     scaled_matrix = current_matrix.copy()
#                     local_scaling(scaled_matrix, v, k, d)
#                     encoded_scaled = bitpack_encode(scaled_matrix, d)
#                     if (encoded_scaled not in visited) or (encoded_scaled not in current_level):
#                         perms=generate_permuted(scaled_matrix,d)
#                         encoded_scaled = min(perms)
#                         if encoded_scaled not in visited:
#                             temp_level.add(encoded_scaled)
#                             visited.add(encoded_scaled)
                
#                 # Local complementation
#                 for k in complementing_factors:
#                     complemented_matrix = current_matrix.copy()
#                     local_complementation(complemented_matrix, v, k, d)
#                     encoded_complemented = bitpack_encode(complemented_matrix, d)
#                     if (encoded_complemented not in visited) or (encoded_complemented not in current_level):
#                         perms=generate_permuted(complemented_matrix,d)
#                         encoded_complemented = min(perms)
#                         if encoded_complemented not in visited:
#                             temp_level.add(encoded_complemented)
#                             visited.add(encoded_complemented)

#             while temp_level:
#                 g=temp_level.pop()
#                 queue.append(g)



def orbit_search(n,d):
    # Read the graph6 file
    ts=time.time()
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
        print(f"Before orbit: {time.time()-ts}")
        temp_orbit, current_graph = find_orbit(current_graph,n,d)
        orbits.append(current_graph)
        print(f"{current_graph}:{len(temp_orbit)}")
        print(f"After orbit: {time.time()-ts}")
        for g in temp_orbit:
           graphs.difference_update(generate_permuted(bitpack_decode(g, n, d), d))
        temp_orbit=set(temp_orbit)
        orbits_reps.append(current_graph)
        orbits.append(temp_orbit)
        
        

    print(f"End: {time.time()-ts}")
    return orbits, orbits_reps

async def orbit_search_isomorphic(filename,n,d):
    # Read the graph6 file
    ts=time.time()
    if d==2:
        with open(filename, "r") as file:
            graph6_lines = [line.strip() for line in file if line.strip()]

        graphs = set([bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines])
    else:
        graphs=set(await generate_graphs(filename,n,d))
    orbits=[]
    while graphs:
        current_graph = graphs.pop()
        print(f"Before orbit: {time.time()-ts}")
        temp_orbit, current_graph = find_orbit_isomorphic(graphs,current_graph,n,d)
        orbits.append(current_graph)
        print(f"{current_graph}:{len(temp_orbit)}")
        print(f"After orbit: {time.time()-ts}")        

    print(f"End: {time.time()-ts}")
    return orbits

# def orbit_search_isomorphic_from_file(n,d):
#     # Read the graph6 file
#     graphs = set()
#     ts=time.time()
#     # List all files in the directory
#     directory=os.path.join(os.getcwd(), str(n))
#     for filename in os.listdir(directory):
#         file_path = os.path.join(directory, filename)

#         # Ensure we only read files (skip directories)
#         if os.path.isfile(file_path):
#             with open(file_path, "r") as file:
#                 for line in file:
#                     graph = int(line.strip())  # Convert graph from string to integer
#                     graphs.add(graph)
#     orbit_reps=[]
#     orbits=[]
#     while graphs:
#         current_graph = graphs.pop()
#         print(f"Before orbit: {time.time()-ts}")
#         orbit, min_graph = find_orbit_isomorphic(current_graph,graphs,n,d)
#         print(f"{min_graph}:{len(orbit)}:{orbit}")
#         orbit_reps.append(min_graph)
#         orbits.append(orbit)
#         print(f"After orbit: {time.time()-ts}")
#         graphs.difference_update(orbit)
#     print(len(orbit_reps))
#     print(f"End: {time.time()-ts}")
#     return orbits, orbit_reps
def orbit_search_isomorphic_from_file(n,d):
    # Read the graph6 file
    graphs = set()
    ts=time.time()
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
    orbit_reps=[]
    orbits=[]
    while graphs:
        current_graph = graphs.pop()
        print(f"Before orbit: {time.time()-ts}")
        orbit= find_orbit(current_graph,n,d)
        
        
       
        print(f"After orbit: {time.time()-ts}")
        orbit_filtered=set()
        while orbit:
            g=orbit.pop()
            perms=generate_permuted(bitpack_decode(g,n,d),d)
            orbit_filtered.add(min(perms))
            orbit.difference_update(perms)
            graphs.difference_update(perms)
        orbits.append(orbit_filtered)
        min_graph=min(orbit_filtered)
        orbit_reps.append(min_graph)
        
        
    print(len(orbit_reps))
    print(f"End: {time.time()-ts}")
    return orbits, orbit_reps


def create_complementation_layer(current,n,d):
    scaling_factors = range(2, d)
    complementing_factors = range(1, d)
    sub_G=nx.Graph()
    sub_G.add_node(current)
    current_matrix=bitpack_decode(current,n,d)
    for v in range(n):
        # Local scaling
        for k in scaling_factors:
            scaled_matrix = current_matrix.copy()
            local_scaling(scaled_matrix, v, k, d)
            encoded_scaled = bitpack_encode(scaled_matrix, d)
            if encoded_scaled not in sub_G.nodes:
                perms=generate_permuted(scaled_matrix,d)
                encoded_scaled=min(perms)
                sub_G.add_node(encoded_scaled)
                sub_G.add_edge(current,encoded_scaled)

        # Local complementation
        for k in complementing_factors:
            complemented_matrix = current_matrix.copy()
            local_complementation(complemented_matrix, v, k, d)
            encoded_complemented = bitpack_encode(complemented_matrix, d)
            if encoded_complemented not in sub_G.nodes:
                perms=generate_permuted(complemented_matrix,d)
                encoded_complemented=min(perms)
                sub_G.add_node(encoded_complemented)
                sub_G.add_edge(current,encoded_complemented)
    return sub_G

def orbit_atlas(n,d):
    graphs=set()
    ts=time.time()
    directory=os.path.join(os.getcwd(), str(n))
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Ensure we only read files (skip directories)
        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                for line in file:
                    graph = int(line.strip())  # Convert graph from string to integer
                    graphs.add(graph)
    G=nx.Graph()
    G.add_nodes_from(graphs)

    num_workers = multiprocessing.cpu_count()
    with multiprocessing.Pool(num_workers) as pool:
        while graphs:
            batch_size = min(len(graphs), num_workers)
            batch = [graphs.pop() for _ in range(batch_size)]
            print(f"{len(graphs)}:{time.time()-ts}")
            subgraphs = pool.starmap(create_complementation_layer, [(g, n, d) for g in batch])
            for sub_G in subgraphs:
                G.update(sub_G)
    return G


def call_complementation_layer(n, d, encoded_value):
    """ Calls the C program and returns results as a set. """
    result = subprocess.run(
        ["./c/a.exe", str(n), str(d), str(encoded_value)],  # Convert args to strings
        capture_output=True, text=True
    )

    # Convert output lines to a set of integers
    return {int(line) for line in result.stdout.splitlines()}


def orbit_atlas_c(n,d):
    graphs=set()
    ts=time.time()
    directory=os.path.join(os.getcwd(), str(n))
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Ensure we only read files (skip directories)
        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                for line in file:
                    graph = int(line.strip())  # Convert graph from string to integer
                    graphs.add(graph)
    G=nx.Graph()
    G.add_nodes_from(graphs)

    num_workers = multiprocessing.cpu_count()
    with multiprocessing.Pool(num_workers) as pool:
        while graphs:
            batch_size = min(len(graphs), num_workers)
            batch = [graphs.pop() for _ in range(batch_size)]
            print(f"{len(graphs)}:{time.time()-ts}")
            complements = pool.starmap(call_complementation_layer, [(n, d, g) for g in batch])
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
n=7
d=3
# if __name__ == "__main__":
#     multiprocessing.freeze_support()
#     g=orbit_atlas_c(n,d)
#     nx.write_edgelist(g,f"orbits_d{d}_n{n}", data=False)


#orbit_search_isomorphic_from_file(n,d)
ts=time.time() 


g=nx.read_edgelist(f"orbits_d{d}_n{n}", nodetype=int)
print(time.time()-ts)
orbits=separate_orbits(g)
print(time.time()-ts)
directory = f"orbits_d{d}_n{n}_separated"
os.makedirs(directory, exist_ok=True)
for i in range(len(orbits)):
    nx.write_edgelist(orbits[i],f"{directory}/orbit_{i}", data=False)

# print(time.time()-ts)
#for o in orbits:
#   plot_graph(o,n,d)

#orbit_atlas_c(n,d)
