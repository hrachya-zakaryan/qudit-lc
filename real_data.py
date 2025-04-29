import networkx as nx
import numpy as np
from collections import deque, Counter
import matplotlib.pyplot as plt
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph
from itertools import permutations, combinations
from data_analysis_functions import *
import math
from isomorphism_method import *
import os
import multiprocessing
import subprocess #get  C file
import concurrent.futures 
import time
# from orbital_graphs import create_orbital_graph, orbit_search_isomorphic_from_file, circular_subgraph_layout, draw_subgraph_inside_circle,plot_graph



#TODO In mathematica what is the difference between graph plot and the plotting of calculate graph
n=3
d=3

# # This will be used for the general 
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


# g=nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_0", nodetype=int)
# print(len(list(g.nodes())))
# print(len(list(g.edges())))
# adj = nx.to_numpy_array(g)
# print(adj)

def bfs_og_mean_path_length(v,g):
    level={v}
    dist=1
    temp_g=g.copy()
    temp_sum=0
    while level:
        current_level = set(level)
        level = set()
        for vertex in current_level:
            if vertex in temp_g.nodes:
                neighbors = set([i for i in nx.neighbors(temp_g,vertex)])
                temp_g.remove_node(vertex)
                if vertex in neighbors:
                    neighbors.remove(vertex)
                neighbors.difference_update(current_level)
                neighbors.difference_update(level)
                temp_sum+=len(neighbors)*dist
                level.update(neighbors)
                
        dist+=1
    return temp_sum


parallel_bfs_og_mean_dist = False
if parallel_bfs_og_mean_dist:
    # Read the graph
    g = nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_0", nodetype=int)
    # ts_new = time.time()
    # print("netwrokx avg:",avg_of_distance_matrix(g))    
    # print("networkx:", time.time()-ts_new)    
    if __name__== "__main__":
        multiprocessing.freeze_support()
        sorted_nodes = sorted(g.nodes)
        relabel_dict = {node: i for i, node in enumerate(sorted_nodes)}
        g = nx.relabel_nodes(g,relabel_dict)
        total_sum=0
        ts = time.time()
        print("number of nodes:",len(g.nodes))
        num_workers = multiprocessing.cpu_count()
        nodes=list(g.nodes)
        print(num_workers)
        with multiprocessing.Pool(num_workers) as pool:
            while nodes:
                print(time.time()-ts)
                batch_size = min(len(nodes),num_workers)
                batch = [nodes.pop(0) for _ in range(batch_size)]
                res_bfs = pool.starmap(bfs_og_mean_path_length,[(v,g) for v in batch])
                total_sum+=sum(res_bfs)
        N=len(g.nodes)
        mean_distance=total_sum/(N*(N-1))
        print("mean distance:",mean_distance)
        print("bfs time", time.time()-ts)       
    






# # Get the adjacency matrix
# adj = nx.to_numpy_array(g)
# np.savetxt("adjacency_matrix.txt", adj, fmt='%.0f')




# dist_mat = distance_matrix(g)
# plot_distance_matrix(g)
# test = avg_of_distance_matrix(g)
# print(test)
     

# #TODO: Check if the weights are really implemented
# #TODO Check if the functions bellow really work for weighted graphs. Check the results with mathematica also.



calculate_og_chromatic_numbers = False #Done
calculate_og_mean_distance_matrix = False
calculate_og_is_planar = False #Done
calculate_og_number_of_loops = False #Done
calculate_og_number_of_circles = False
calculate_og_minimum_chromatic_number_in_OG = False #Done
less_than_10_mb = False
all_data = False



if calculate_og_chromatic_numbers:
    array_og_chromatic_number = []
    for cnt in range(3,n+1):
        print(f"Start calculating for array_og_chromatic_number n={cnt}")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            og_chromatic_number = chromatic_number(g)
            array_og_chromatic_number.append(og_chromatic_number)
    with open("og_data/og_chromatic_numbers.txt", "w") as f:
        f.write(str(array_og_chromatic_number))
        
def calculate_og_mean_path_length(filename,count):
    file_path = f"orbits_d{d}_n{count}_separated/{filename}"
    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
    if file_size_mb <= 10 :
        print(f"n = {count} and {filename}")
        g=nx.read_edgelist(file_path, nodetype=int)
        og_mean_distance_matrix = avg_of_distance_matrix(g)
        return og_mean_distance_matrix
    else:
        print(f"For n = {count} and {filename} mean distance not calculated")
        return str(filename)
        
        
if calculate_og_mean_distance_matrix and all_data:
    array_og_mean_distance_matrix = []
    array_not_cal_og_mean_distance_matrix = []
    for cnt in range(3,n+1):
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        for file in files_for_given_orbit:
            print(f"n = {cnt} and {file}")
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            og_mean_distance_matrix = avg_of_distance_matrix(g)
            print(f"avg of distance matrix: {og_mean_distance_matrix}")
            array_og_mean_distance_matrix.append(og_mean_distance_matrix)
    with open(f"og_data/new_og_mean_distance_matrix_d{d}_n{n}.txt", "w") as f:
        f.write(str(array_og_mean_distance_matrix))




if calculate_og_is_planar:
    array_og_is_planar = []
    for cnt in range(3,n+1):
        print(f"Start calculating for array_og_is_planar n={cnt}")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            og_is_planar = is_planar_graph(g)
            array_og_is_planar.append(og_is_planar)
        
            
    with open("og_data/og_is_planar.txt", "w") as f:
        f.write(str(array_og_is_planar))        
        


if calculate_og_number_of_loops:
    array_og_number_of_loops = []
    for cnt in range(3,n+1):
        print(f"Start calculating for array_og_number_of_loops n={cnt}")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            og_number_of_loops = count_self_loops(g)
            array_og_number_of_loops.append(og_number_of_loops)
    
    
    with open("og_data/og_number_of_loops.txt", "w") as f:
        f.write(str(array_og_number_of_loops))        
        

    

if calculate_og_number_of_circles:
    array_og_number_of_circles = []
    for cnt in range(3,n+1):
        print(f"Start calculating for array_og_number_of_circles n={cnt}")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            og_number_of_circles = count_cycles(g)
            array_og_number_of_circles.append(og_number_of_circles)
    with open("og_data/og_number_of_circles.txt", "w") as f:
        f.write(str(array_og_number_of_circles))        



    

if calculate_og_minimum_chromatic_number_in_OG:
    array_og_minimum_chromatic_number_in_OG = []
    for cnt in range(3,n+1):
        print(f"Start calculating for array_og_minimum_chromatic_number_in_OG n={cnt}")
        files_for_given_orbit = os.listdir(f"orbits_d{d}_n{cnt}_separated")
        for file in files_for_given_orbit:
            g=nx.read_edgelist(f"orbits_d{d}_n{cnt}_separated/{file}", nodetype=int)
            og_minimum_chromatic_number_in_OG = min_chromatic_number_in_OG(g,n,d)
            array_og_minimum_chromatic_number_in_OG.append(og_minimum_chromatic_number_in_OG)
    with open("og_data/og_minimum_chromatic_number_in_OG.txt", "w") as f:
        f.write(str(array_og_minimum_chromatic_number_in_OG))        
        
        
        
        
        
        



# # plot_graph(g,n,d)
# print("------------------------------- Checked Observables -------------------------------")
# print(f"------------------------------- d={d}, n={3} -------------------------------")

# #1
# og_chromatic_number = chromatic_number(g)
# print('chromatic number of G:',og_chromatic_number)


# # og_distance_matrix = distance_matrix(g)

# # plot_distance_matrix(g)
# #2
# og_avg_distance_matrix = average_of_distance_matrix(g)
# print("avg distance matrix:", og_avg_distance_matrix)

# #3
# og_max_distance_matrix = max_in_distance_matrix(g)
# print("max in OG distance matrix:", og_max_distance_matrix)

# og_planar = is_planar_graph(g)
# print("OG is planar:", og_planar)

# #4
# og_number_of_loops = count_self_loops(g)
# print("og number of loops:", og_number_of_loops)

# #5
# number_of_circles = count_cycles(g)
# print("number of circles:", number_of_circles)

# print("has Eulerian circle:", has_eulerian_cycle(g))

# #6 OG chromatic index comes from mathematica

# #7 Minimum number of colorability among the graphs belonging the SAME OG
# minimum_chromatic_number_in_OG = min_chromatic_number_in_OG(g,n,d)
# print("min_chromatic_number_in_OG:",minimum_chromatic_number_in_OG)
