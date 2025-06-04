"""
Script to compute the observables for each entanglement class.

This script processes graph orbits stored as edgelist files generated from Mathematica-based enumeration.
It computes approximate diameters using a Monte Carlo BFS sampling method and optionally evaluates additional
graph properties such as max degree, planarity, tree structure, and more.

Features:
- Efficient approximation of graph diameter via parallel landmark-based BFS
- Optional generation of Mathematica-readable graph input files
- Extraction of orbit-wide properties including minimum of maximum degrees
- Modular toggle for computing individual observables
- Scalable to large datasets using multiprocessing

Dependencies:
- networkx
- numpy
- matplotlib
- scipy
- multiprocessing
- Custom modules: auxilary_unparallelized, data_analysis_functions, isomorphism_method

Inputs:
- Directory structure: "orbits_d{d}_n{cnt}_separated/"
  Each file is an edgelist corresponding to one orbit of graphs.

Outputs:
- Monte Carlo estimates: "og_data_final/mc_{sample_size}_og_diameters_n{n}.txt"
- Additional observables in "og_data_final/", e.g., minimum max degrees per orbit

Usage:
- Adjust `n` and `d` to select the local dimension and number of particles.
- Set relevant computation flags (e.g., `calculate_og_is_tree`) to True.
- Make sure the input directories are structured correctly and accessible.

Notes:
- Some observables require decoding packed integer encodings of graphs.
- Mathematica files are only generated if `produce_mathematica_files` is set.
"""



import networkx as nx
import numpy as np
from collections import deque, Counter
import matplotlib.pyplot as plt
from auxilary_functions import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph
from itertools import permutations, combinations
from data_analysis_functions import *
import math
import random
import os
import multiprocessing
import subprocess #get  C file
import time
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from concurrent.futures import ProcessPoolExecutor, as_completed

n=7
d=3

# #START the Monte Carlo simulation for the diameter

# def bfs_max_distance(graph, start):
#     """
#     Perform BFS from a starting node and return the maximum distance found.
#     """
#     visited = set()
#     queue = deque([(start, 0)])
#     max_distance = 0

#     while queue:
#         node, dist = queue.popleft()

#         if node in visited:
#             continue

#         visited.add(node)
#         max_distance = max(max_distance, dist)

#         for neighbor in graph.neighbors(node):
#             if neighbor not in visited:
#                 queue.append((neighbor, dist + 1))

#     return max_distance

# def parallel_bfs_max_distance(args):
#     """
#     Wrapper function for parallel execution of bfs_max_distance.
#     """
#     graph, start = args
#     return bfs_max_distance(graph, start)

# def approximate_diameter(graph, cpus_used, sample_size=None):
#     """
#     Approximates the graph diameter using landmark-based BFS with parallel computation.

#     Parameters:
#     graph (networkx.Graph): Input graph.
#     cpus_used (int): Number of CPU cores to use for parallelization.
#     sample_size (int): Number of landmark nodes to sample (default: sqrt(n)).

#     Returns:
#     int: Approximate diameter of the graph.
#     """
#     nodes = list(graph.nodes())

#     # Determine the number of landmarks (default: sqrt(n) if not specified)
#     if sample_size is None:
#         sample_size = int(len(nodes) ** 0.5)

#     # Randomly select landmark nodes
#     landmarks = random.sample(nodes, min(sample_size, len(nodes)))

#     # Use multiprocessing for parallel BFS computation
#     with multiprocessing.Pool(cpus_used) as pool:
#         results = pool.map(parallel_bfs_max_distance, [(graph, node) for node in landmarks])

#     # Return the maximum distance as the approximate diameter
#     return max(results)

# # Usage 
# if __name__ == "__main__":
#     array_diameter_mc = []
#     cpus_to_use = 7
#     size_sample = 1000
#     for cnt in range(3, n+1):
#         files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
#         files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
#         for file in files_for_given_orbit:
#             print(f"Starting orbits_d{d}_n{cnt}_separated/{file}")
#             g = nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
#             ts = time.time()
#             approx_diameter_value = approximate_diameter(g, cpus_to_use, size_sample)
#             array_diameter_mc.append(approx_diameter_value)
#             print(f"Approx Diameter of Starting orbits_d{d}_n{cnt}_separated/{file}: {approx_diameter_value}")
#             print(f"Time of calc of orbits_d{d}_n{cnt}_separated/{file}: {time.time() - ts}")
    
#     with open(f"og_data/mc_{size_sample}_og_diameters_n{n}.txt", "w") as f:
#         f.write(str(array_diameter_mc))
#     print(array_diameter_mc)


# approx_diameter_value = approximate_diameter(g, sample_size=100000, cpus_used=2)
# print(f"Approx ASPL : {approx_diameter_value}")
# print("Time of calc:", time.time()-ts)


#END the Monte Carlo simulation for the diameter

#START: Produce mathematica files
# We have to adapt the input of the orbit files so that we can just use it in mathematica 
produce_mathematica_files = False

if produce_mathematica_files:
    files_for_given_orbit = os.listdir(f"orbits_d{d}_n{n}_separated")
    for file in files_for_given_orbit:
        #insert the encoded orbits
        g=nx.read_edgelist(f"orbits_d{d}_n{n}_separated/{file}", nodetype=int)
        encoded_graphs_in_orbit = list(g.nodes())
        graphs_of_orbit_in_G_form = [bitpack_decode(encoded_graphs_in_orbit[i], n, d) for i in range(len(encoded_graphs_in_orbit))]
        file_name =f"files_for_mathematica/graphs_in_orbits_d{d}_n{n}_separated_{file}.txt" 
        with open(file_name, "w") as f:
            # Loop through the list of graphs (or matrices) and write them in the desired format
            for graph in graphs_of_orbit_in_G_form:
                # Convert each matrix (graph) to the correct format
                matrix_str = "{ " + ", ".join([f"{{{', '.join(map(str, row))}}}" for row in graph]) + " }"
                f.write(matrix_str + "\n\n")

#END: Produce mathematica files

#Start the calculation of the observables and produce the files to enter in the data analysis notebook
calculate_og_chromatic_numbers = False #Final Done
calculate_og_mean_distance_matrix = False #Final Done (Statistically Proven the Approximation is pretty close to reality)
calculate_og_is_planar = False #Final Done
calculate_og_number_of_loops = False #Final Done
calculate_og_minimum_chromatic_number_in_OG = False #Final Done
less_than_10_mb = False #Final Done
all_data = False #Done statistically
calculate_og_density = False #Final Done
calculate_og_is_tree = False #Final Done
calculate_max_degree = True #Final Done
calcualte_min_of_max_degree_in_OG = False #Final Done


if calcualte_min_of_max_degree_in_OG:
    array_min_of_max_degree_in_OG = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated_sorted")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            # Load the original graph (assumed to be a standard graph with self-loops)
            # file_path = f"orbits_d{d}_n{n}_separated/orbit_0"
            og = nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated_sorted/{file}", nodetype=int)
            # Get a list of all node indices
            node_indices = list(og.nodes)
            # print("Number of nodes:", len(node_indices))
            # print("Node indices:", node_indices)

            list_of_max_degrees = []
            for total_nodes in range(len(node_indices)):
                # Select a specific node 
                specific_graph = node_indices[total_nodes]

                # Convert it to a MultiGraph
                multi_graph_ex = decode_to_net_G(specific_graph, n, d)

                # Compute and print the degree of nodes (accounting for multiple edges)
                max_degree = max(dict(multi_graph_ex.degree()).values())
                # print("max degree:",max_degree)
                list_of_max_degrees.append(max_degree)
                # draw_graph(specific_graph,n,d)

            min_of_max_degrees = min(list_of_max_degrees)
            array_min_of_max_degree_in_OG.append(min_of_max_degrees)

            print(f"min of max degree for {file}: {min_of_max_degrees}")
    with open("og_data_final/join_min_of_max_degree_in_OG.txt", "w") as f:
        f.write(str(array_min_of_max_degree_in_OG))




if calculate_max_degree:
    array_og_max_degree = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated_sorted")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated_sorted/{file}", nodetype=int)
            temp = max_degree(g)
            array_og_max_degree.append(temp)
            print(f"OG max for {file}: {temp}")
    with open("og_data_final/join_og_max_degree.txt", "w") as f:
        f.write(str(array_og_max_degree))

if calculate_og_is_tree:
    array_og_is_tree = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated_sorted")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated_sorted/{file}", nodetype=int)
            temp = is_tree(g)
            array_og_is_tree.append(temp)
            print(f"OG is tree for {file} found")
    with open("og_data_final/join_og_is_tree.txt", "w") as f:
        f.write(str(array_og_is_tree))
    
    




if calculate_og_chromatic_numbers and calculate_og_is_planar and calculate_og_number_of_loops and calculate_og_minimum_chromatic_number_in_OG and calculate_og_density:
    array_og_chromatic_number = []
    array_og_is_planar = []
    array_og_number_of_loops = []
    array_og_minimum_chromatic_number_in_OG = []
    array_og_density = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated_sorted")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated_sorted/{file}", nodetype=int)
            og_chromatic_number = chromatic_number(g)
            array_og_chromatic_number.append(og_chromatic_number)
            print(f"OG chromatic number for {file} found: {og_chromatic_number}")
            og_is_planar = is_planar_graph(g)
            array_og_is_planar.append(og_is_planar)
            print(f"OG is planar for {file} found: {og_is_planar}")
            og_number_of_loops = count_self_loops(g)
            array_og_number_of_loops.append(og_number_of_loops)
            print(f"OG number of loops for {file} found: {og_number_of_loops}")
            og_minimum_chromatic_number_in_OG = min_chromatic_number_in_OG(g,n,d)
            array_og_minimum_chromatic_number_in_OG.append(og_minimum_chromatic_number_in_OG)
            print(f"Minimum chromatic number in OG for {file} found: {og_minimum_chromatic_number_in_OG}")
            den = nx.density(g)
            array_og_density.append(den)
            print(f"Density of OG for {file} found: {den}")
                
            
    with open("og_data_final/join_og_chromatic_numbers.txt", "w") as f:
        f.write(str(array_og_chromatic_number))
    with open("og_data_final/join_og_is_planar.txt", "w") as f:
        f.write(str(array_og_is_planar)) 
    with open("og_data_final/join_og_number_of_loops.txt", "w") as f:
        f.write(str(array_og_number_of_loops)) 
    with open("og_data_final/join_og_minimum_chromatic_number_in_OG.txt", "w") as f:
        f.write(str(array_og_minimum_chromatic_number_in_OG))
    with open("og_data_final/join_og_density.txt", "w") as f:
        f.write(str(array_og_density))
        
        
        
       
                    
                    
                     
if calculate_og_mean_distance_matrix and less_than_10_mb:
    array_og_mean_distance_matrix = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated_sorted")
        files_for_given_orbit.sort(key=lambda x: int(x.split('_')[-1]))  # Sort the list based on the orbit number
        for file in files_for_given_orbit:
            file_size_mb = os.path.getsize(f"orbits_d{d}_n{cnt}_separated_sorted/{file}") / (1024 * 1024)
            if file_size_mb <= 10:
                print(f"Star Calculating Mean Path Distance of OG for {file}")
                g = nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated_sorted/{file}", nodetype=int)
                og_mean_distance_matrix = avg_of_distance_matrix(g)
                array_og_mean_distance_matrix.append(og_mean_distance_matrix)
                print(f"Result: {og_mean_distance_matrix}")
            else:
                print(f"Skip Calculating Mean Path Distance of OG for {file} due to file size. Zero assinged")
                array_og_mean_distance_matrix.append(0)
    with open("og_data_final/join_og_mean_path_length.txt", "w") as f:
        f.write(str(array_og_mean_distance_matrix))
        
if calculate_og_mean_distance_matrix and all_data:
    array_og_mean_distance_matrix = []
    array_not_cal_og_mean_distance_matrix = []
    for cnt in range(3,n+1):
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated_sorted")
        for file in files_for_given_orbit:
            print(f"n = {cnt} and {file}")
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated_sorted/{file}", nodetype=int)
            og_mean_distance_matrix = avg_of_distance_matrix(g)
            print(f"avg of distance matrix: {og_mean_distance_matrix}")
            array_og_mean_distance_matrix.append(og_mean_distance_matrix)
    with open(f"og_data_final/new_og_mean_distance_matrix_d{d}_n{n}.txt", "w") as f:
        f.write(str(array_og_mean_distance_matrix))
  

                

