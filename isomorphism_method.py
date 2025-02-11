import networkx as nx
from collections import deque
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph
from itertools import permutations
import numpy as np
from itertools import product
import time
from networkx.algorithms.isomorphism import GraphMatcher



def generate_non_isomorphic_graphs(base_graph):
    edges = list(base_graph.edges)
    unique_graphs = set()
    unique_graphs_bitpacked = []
    for i in range(2 ** len(edges)):
        # Generate a weight combination
        weight_config = [(edges[j], 1 + ((i >> j) & 1)) for j in range(len(edges))]

        # Create a new graph with this weight configuration
        G = nx.Graph(base_graph)
        for edge, weight in weight_config:
            G[edge[0]][edge[1]]['weight'] = weight

        # Get the canonical form using isomorphism check
        can_form = nx.convert_node_labels_to_integers(G)
        if not any(GraphMatcher(G, other,edge_match=lambda x,y: x['weight']==y['weight']).is_isomorphic() for other in unique_graphs):
            unique_graphs.add(G)
            unique_graphs_bitpacked.append(bitpack_encode(nx.to_numpy_array(G,dtype=int),d))
    return unique_graphs_bitpacked


def create_graphs(filename, n, d):
    
    with open(filename, "r") as file:
        graph6_lines = [line.strip() for line in file if line.strip()]

    iso_graphs = [bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines]

    weighted_iso_graphs=[]
    for ig in iso_graphs:
        weighted_iso_graphs.extend(generate_non_isomorphic_graphs(nx.from_numpy_array(bitpack_decode(ig,n,d))))
    return weighted_iso_graphs



def generate_permuted(adj_matrix,d):
    """Generate all possible vertex permutations of an adjacency matrix."""
    n = len(adj_matrix)
    permuted_matrices = set()
    
    for perm in permutations(range(n)):
        permuted_matrix = adj_matrix[np.ix_(perm, perm)]
        permuted_matrices.add(bitpack_encode(permuted_matrix,d))
        
    return list(permuted_matrices)


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
        if all(weight == 1 for weight in weights):
            continue
        new_matrix = adj_matrix.copy()
        # Assign weights symmetrically
        for (pos, weight) in zip(upper_indices, weights):
            i, j = pos
            new_matrix[i, j] = weight
            new_matrix[j, i] = weight
        weighted_matrices.append(bitpack_encode(new_matrix,d))

    return weighted_matrices


def generate_graphs(filename, n, d):
    
    with open(filename, "r") as file:
        graph6_lines = [line.strip() for line in file if line.strip()]

    iso_graphs = [bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines]


    weighted_iso_graphs=[]
    for ig in iso_graphs:
        weighted_iso_graphs.append(ig) 
        weighted_graphs = generate_weighted(bitpack_decode(ig,n,d),d)
        i=0
        while i<len(weighted_graphs):
            weighted_iso_graphs.append(weighted_graphs[i])
            perms=generate_permuted(bitpack_decode(weighted_graphs[i],n,d),d)
            weighted_graphs=[x for x in weighted_graphs if x not in perms]
            i+=1
        weighted_iso_graphs.append(weighted_graphs[-1])
    return weighted_iso_graphs


def total_weight(graph):
    bit_length=graph.bit_length()
    pair_sum = []
    for i in range(0, bit_length, 2):
        # Extract the current pair of bits
        pair = (graph >> i) & 0b11  # Mask the last two bits
        pair_sum.append(bin(pair).count('1'))  # Count the number of set bits in the pair

    return pair_sum[::-1]  # Reverse to maintain left-to-right order

def find_orbit(start_graph, n, d):
    queue = deque([start_graph])  # BFS queue
    visited = []  # Set of visited graphs
    visited.append(start_graph)
    
    # Define all possible scaling and complementation factors
    scaling_factors = range(2, d)  # Local scaling factors (mod d)
    complementing_factors = range(1, d)
    min_graph=visited[0]

    while queue:
        current = queue.popleft()
        current_matrix = bitpack_decode(current, n, d)
        if current.bit_count()<=min_graph.bit_count() and current<min_graph:
            min_graph=current
      
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
    return visited, min_graph


def orbit_search(filename,n,d):
    # Read the graph6 file
    ts=time.time()
    if d==2:
        with open(filename, "r") as file:
            graph6_lines = [line.strip() for line in file if line.strip()]

        graphs = set([bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines])
    else:
        graphs=set(generate_graphs(filename,n,d))
    orbits=[]
    full_orbits_filtered=[]
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
        temp_orbit_filtered=set()
        while temp_orbit:
            g = temp_orbit.pop()
            temp_orbit_filtered.add(g)
            temp_orbit.difference_update(generate_permuted(bitpack_decode(g,n,d),d))
        full_orbits_filtered.append(temp_orbit_filtered)

    print(f"End: {time.time()-ts}")
    return orbits, full_orbits_filtered

# Convert each line from graph6 to a NetworkX graph
n=6
d=3
#o,f_o=orbit_search("d3n4.txt",n,d)

ts=time.time()
print(time.time()-ts)
gg=generate_graphs("d3n6.txt",n,d)
print(len(gg))
print(time.time()-ts)
ts=time.time()
print(time.time()-ts)
cg=create_graphs("d3n6.txt",n,d)
print(len(cg))
print(time.time()-ts)

#for i in range(len(cg)):
#    draw_graph(gg[i],n,d)
#    draw_graph(cg[i],n,d)