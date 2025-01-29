from auxilary_unparallelized import bitpack_decode, bitpack_encode, create_graphs, local_scaling, local_complementation, draw_graph, save_to_file, load_from_file, save_orbits_to_file, save_representatives_to_file, load_orbits_from_file
from collections import deque
import numpy as np
from itertools import product, combinations


def find_orbit(start_graph, n, d):
    queue = deque([start_graph])  # BFS queue
    visited = []  # Set of visited graphs
    visited.append(start_graph)
    
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
                    visited.append(encoded_scaled)
                    queue.append(encoded_scaled)
            
            # Local complementation
            for k in complementing_factors:
                complemented_matrix = current_matrix.copy()
                local_complementation(complemented_matrix, v, k, d)
                encoded_complemented = bitpack_encode(complemented_matrix, d)
                if encoded_complemented not in visited:
                    visited.append(encoded_complemented)
                    queue.append(encoded_complemented)
    #print(len(visited))
    return visited  # Return all visited graphs as the result


def full_classification(n, d):
    # Generate the initial list of graphs
    graphs = create_graphs(n, d)
    classifications = []  # To store equivalence classes

    while graphs:
        start_graph = graphs[0]  # Take the first graph in the list
        visited = find_orbit(start_graph, n, d)  # Perform BFS from this graph
        classifications.append(visited)  # Add the equivalence class to the results
        
        # Remove all visited graphs from the original list
        graphs = [graph for graph in graphs if graph not in visited]

    return classifications


def full_classification_loaded(n, d, filename):
    # Generate the initial list of graphs
    graphs = load_from_file(filename)
    classifications = []  # To store equivalence classes

    while graphs:
        start_graph = graphs[0]  # Take the first graph in the list
        visited = find_orbit(start_graph, n, d)  # Perform BFS from this graph
        classifications.append(visited)  # Add the equivalence class to the results
        
        # Remove all visited graphs from the original list
        graphs = [graph for graph in graphs if graph not in visited]

    return classifications


def check_orbit_completeness(orbits, n, d):
    """
    Check if all graphs with the same structure (weight-1 edges only) and all weight combinations
    are present in each orbit.

    Parameters:
        orbits (list[set[int]]): List of orbits, where each orbit is a set of encoded graphs.
        n (int): Number of vertices in the graph.
        d (int): Modulo value for weights.

    Returns:
        list[tuple[int, bool]]: List of tuples where each entry contains the orbit index
                                and a boolean indicating if the orbit is complete.
    """
    incomplete_orbits = []

    for orbit_index, orbit in enumerate(orbits):
        structure_dict = {}

        for encoded_graph in orbit:
            # Decode the graph and determine its structure (weight-1 edges)
            adjacency_matrix = bitpack_decode(encoded_graph, n, d)
            structure_matrix = (adjacency_matrix > 0).astype(int)  # Binary structure matrix
            structure_key = bitpack_encode(structure_matrix, 2)  # Encode structure as a key

            # Store all graphs for this structure in a set
            if structure_key not in structure_dict:
                structure_dict[structure_key] = set()
            structure_dict[structure_key].add(encoded_graph)

        # Check if all weight combinations for each structure are present
        for structure_key, encoded_graphs in structure_dict.items():
            # Generate all possible weight combinations for this structure
            structure_matrix = bitpack_decode(structure_key, n, 3)  # Decode binary structure
            num_edges = np.sum(np.triu(structure_matrix))  # Count edges in the structure

            # Generate all combinations of weights
            for weight_combination in product(range(1,d), repeat=num_edges):
                idx = 0
                weighted_matrix = structure_matrix.copy()
                for i in range(n):
                    for j in range(i + 1, n):
                        if structure_matrix[i, j] == 1:  # If edge exists
                            weighted_matrix[i, j] = weight_combination[idx]
                            weighted_matrix[j, i] = weight_combination[idx]
                            idx += 1

                # Encode and check if this graph is in the orbit
                encoded_weighted_graph = bitpack_encode(weighted_matrix, d)
                if encoded_weighted_graph not in encoded_graphs:
                    incomplete_orbits.append((orbit_index, False))
                    break
            else:
                # All combinations are present for this structure
                continue
            break
        else:
            # All structures and weights are complete
            incomplete_orbits.append((orbit_index, True))

    return incomplete_orbits

n=7
d=2

#res=full_classification(n,d)
#print(len(res))
#for graphs in res:
    #draw_graph(graphs[0],n,d)
#res=full_classification_loaded(n,d,"n6d3")
#save_orbits_to_file(res,filename="n6d2_orbits")
#save_representatives_to_file(res,filename="n_6,d_3_rep")
#res=load_orbits_from_file("n_4,d_3_orbits")
#print(check_orbit_completeness(res,n,d))
#print(bitpack_decode(2725,n,d))
#print(bitpack_decode(1445,n,d))