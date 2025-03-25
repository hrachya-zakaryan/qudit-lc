import networkx as nx
import numpy as np
from collections import deque, Counter
import matplotlib.pyplot as plt
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph
from itertools import permutations, combinations
from data_analysis_functions import *
import math
import random
from isomorphism_method import *
import os
import multiprocessing
import subprocess #get  C file
import time
# from orbital_graphs import create_orbital_graph, orbit_search_isomorphic_from_file, circular_subgraph_layout, draw_subgraph_inside_circle,plot_graph
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from concurrent.futures import ProcessPoolExecutor, as_completed





#TODO: Create a new file with the pure data used and push it to github for safety.
#TODO: Pipeline for plotting the correlation plots for each combination
#TODO: Correlators for each of the data set.
#TODO: Triangle plot for meaningful data


#TODO In mathematica what is the difference between graph plot and the plotting of calculate graph
#TODO: Check if the weights are really implemented
#TODO: Check that I am calculcating with the correct order the elements OR put everything together
#TODO:  Implement the mathematica files for the calculation of the 2 type of chromatic index we need.
#TODO: Is the center of the graph physical?





n=7
d=3

# def approximate_aspl(graph, sample_size=10000):
#     """
#     Approximates the Average Shortest Path Length (ASPL) by sampling node pairs.

#     Parameters:
#     graph (networkx.Graph): Input graph.
#     sample_size (int): Number of node pairs to sample.

#     Returns:
#     float: Approximate ASPL of the graph.
#     """
#     nodes = list(graph.nodes())
#     total_dist = 0
#     valid_pairs = 0

#     for _ in range(sample_size):
#         u, v = random.sample(nodes, 2)  # Randomly select two nodes
#         try:
#             total_dist += nx.shortest_path_length(graph, source=u, target=v)
#             valid_pairs += 1
#         except nx.NetworkXNoPath:
#             continue  # Skip if no path exists

#     return total_dist / valid_pairs if valid_pairs > 0 else float('inf')

# def approximate_diameter(graph, sample_size=10000):
#     """
#     Approximates the graph diameter using Monte Carlo sampling.

#     Parameters:
#     graph (networkx.Graph): Input graph.
#     sample_size (int): Number of node pairs to sample.

#     Returns:
#     int: Approximate diameter of the graph.
#     """
#     nodes = list(graph.nodes())
#     max_distance = 0

#     for _ in range(sample_size):
#         u, v = random.sample(nodes, 2)  # Randomly select two nodes
#         try:
#             dist = nx.shortest_path_length(graph, source=u, target=v)
#             max_distance = max(max_distance, dist)
#         except nx.NetworkXNoPath:
#             continue  # Skip if no path exists

#     return max_distance

# ts = time.time()
# # Usage example
# g = nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_2", nodetype=int)

# approx_diameter_value = approximate_diameter(g,sample_size=100000)
# print(f"Approx ASPL : {approx_diameter_value}")
# print("Time of calc:", time.time()-ts)

# ts = time.time()
# g = nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_0", nodetype=int)
# ex_diam = nx.diameter(g)
# print("Exact diameter",ex_diam)
# print("Time of calc:", time.time()-ts)






#_______________
# def compute_distance(u, v, graph):
#     """
#     Compute the shortest path length between two nodes u and v in the graph.
#     Returns -1 if no path exists.
#     """
#     try:
#         return nx.shortest_path_length(graph, source=u, target=v)
#     except nx.NetworkXNoPath:
#         return -1  # Return -1 if no path exists

# def approximate_diameter(graph, cpus_used, sample_size=10000):
#     """
#     Approximates the graph diameter using Monte Carlo sampling with parallel computation.

#     Parameters:
#     graph (networkx.Graph): Input graph.
#     sample_size (int): Number of node pairs to sample.

#     Returns:
#     int: Approximate diameter of the graph.
#     """
#     nodes = list(graph.nodes())
#     max_distance = 0
#     with multiprocessing.Pool(cpus_used) as pool:
#         # Prepare all pairs to sample
#         pairs = [(random.choice(nodes), random.choice(nodes)) for _ in range(sample_size)]

#         # Use starmap to parallelize the computation of distances
#         results = pool.starmap(compute_distance, [(u, v, graph) for u, v in pairs])

#         # Find the maximum distance from the results
#         max_distance = max(results)

#     return max_distance

# # Usage example
# if __name__ == "__main__":
#     ts = time.time()
#     g = nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_0", nodetype=int)
#     cpus_to_use = 7
#     approx_diameter_value = approximate_diameter(g, cpus_to_use, sample_size=3000000)
#     print(f"Approx Diameter : {approx_diameter_value}")
#     print("Time of calc:", time.time()-ts)
#_______________













def bfs_max_distance(graph, start):
    """
    Perform BFS from a starting node and return the maximum distance found.
    """
    visited = set()
    queue = deque([(start, 0)])
    max_distance = 0

    while queue:
        node, dist = queue.popleft()

        if node in visited:
            continue

        visited.add(node)
        max_distance = max(max_distance, dist)

        for neighbor in graph.neighbors(node):
            if neighbor not in visited:
                queue.append((neighbor, dist + 1))

    return max_distance

def parallel_bfs_max_distance(args):
    """
    Wrapper function for parallel execution of bfs_max_distance.
    """
    graph, start = args
    return bfs_max_distance(graph, start)

def approximate_diameter(graph, cpus_used, sample_size=None):
    """
    Approximates the graph diameter using landmark-based BFS with parallel computation.

    Parameters:
    graph (networkx.Graph): Input graph.
    cpus_used (int): Number of CPU cores to use for parallelization.
    sample_size (int): Number of landmark nodes to sample (default: sqrt(n)).

    Returns:
    int: Approximate diameter of the graph.
    """
    nodes = list(graph.nodes())

    # Determine the number of landmarks (default: sqrt(n) if not specified)
    if sample_size is None:
        sample_size = int(len(nodes) ** 0.5)

    # Randomly select landmark nodes
    landmarks = random.sample(nodes, min(sample_size, len(nodes)))

    # Use multiprocessing for parallel BFS computation
    with multiprocessing.Pool(cpus_used) as pool:
        results = pool.map(parallel_bfs_max_distance, [(graph, node) for node in landmarks])

    # Return the maximum distance as the approximate diameter
    return max(results)

# Usage example
if __name__ == "__main__":
    array_diameter_mc = []
    cpus_to_use = 7
    size_sample = 1000
    for cnt in range(3, n+1):
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            print(f"Starting orbits_d{d}_n{cnt}_separated/{file}")
            g = nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            ts = time.time()
            approx_diameter_value = approximate_diameter(g, cpus_to_use, size_sample)
            array_diameter_mc.append(approx_diameter_value)
            print(f"Approx Diameter of Starting orbits_d{d}_n{cnt}_separated/{file}: {approx_diameter_value}")
            print(f"Time of calc of orbits_d{d}_n{cnt}_separated/{file}: {time.time() - ts}")
    
    with open(f"og_data/mc_{size_sample}_og_diameters_n{n}.txt", "w") as f:
        f.write(str(array_diameter_mc))
    print(array_diameter_mc)












# approx_diameter_value = approximate_diameter(g, sample_size=100000, n_jobs=2)
# print(f"Approx ASPL : {approx_diameter_value}")
# print("Time of calc:", time.time()-ts)





calculate_og_chromatic_numbers = False #Done
calculate_og_mean_distance_matrix = False #Done (Statistically Proven the Approximation is pretty close to reality)
calculate_og_is_planar = False #Done
calculate_og_number_of_loops = False #Done
calculate_og_number_of_circles = False
calculate_og_minimum_chromatic_number_in_OG = False #Done
less_than_10_mb = False #Done
all_data = False #Done statistically
calculate_og_density = False #Done

if calculate_og_chromatic_numbers and calculate_og_is_planar and calculate_og_number_of_loops and calculate_og_minimum_chromatic_number_in_OG and calculate_og_density:
    array_og_chromatic_number = []
    array_og_is_planar = []
    array_og_number_of_loops = []
    array_og_minimum_chromatic_number_in_OG = []
    array_og_density = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            og_chromatic_number = chromatic_number(g)
            array_og_chromatic_number.append(og_chromatic_number)
            print(f"OG chromatic number for {file} found")
            og_is_planar = is_planar_graph(g)
            array_og_is_planar.append(og_is_planar)
            print(f"OG is planar for {file} found")
            og_number_of_loops = count_self_loops(g)
            array_og_number_of_loops.append(og_number_of_loops)
            print(f"OG number of loops for {file} found")
            og_minimum_chromatic_number_in_OG = min_chromatic_number_in_OG(g,n,d)
            array_og_minimum_chromatic_number_in_OG.append(og_minimum_chromatic_number_in_OG)
            print(f"Minimum chromatic number in OG for {file} found")
            den = nx.density(g)
            array_og_density.append(den)
            print(f"Density of OG for {file} found")
                

      
            
                
                
            
    with open("og_data/join_og_chromatic_numbers.txt", "w") as f:
        f.write(str(array_og_chromatic_number))
    with open("og_data/join_og_is_planar.txt", "w") as f:
        f.write(str(array_og_is_planar)) 
    with open("og_data/join_og_number_of_loops.txt", "w") as f:
        f.write(str(array_og_number_of_loops)) 
    with open("og_data/join_og_minimum_chromatic_number_in_OG.txt", "w") as f:
        f.write(str(array_og_minimum_chromatic_number_in_OG))
    with open("og_data/join_og_density.txt", "w") as f:
        f.write(str(array_og_density))
        
        
        
       
                    
                    
                     
if calculate_og_mean_distance_matrix and less_than_10_mb:
    array_og_mean_distance_matrix = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            file_size_mb = os.path.getsize(f"orbits_d{d}_n{cnt}_separated/{file}") / (1024 * 1024)
            if file_size_mb <= 10:
                print(f"Star Calculating Mean Path Distance of OG for {file}")
                g = nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
                og_mean_distance_matrix = avg_of_distance_matrix(g)
                array_og_mean_distance_matrix.append(og_mean_distance_matrix)
            else:
                print(f"Skip Calculating Mean Path Distance of OG for {file} due to file size. Zero assinged")
                array_og_mean_distance_matrix.append(0)
    with open("og_data/join_og_mean_path_length.txt", "w") as f:
        f.write(str(array_og_mean_distance_matrix))
        
        
  

                

