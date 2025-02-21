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
    G=nx.from_numpy_array(adj_matrix)
    node_permutations = list(permutations(G.nodes()))
    
    # Generate relabeled graphs
    relabeled_graphs = []
    for perm in node_permutations:
        mapping = {old: new for old, new in zip(G.nodes(), perm)}
        new_G = nx.relabel_nodes(G, mapping)
        sorted_nodes = sorted(new_G.nodes())
        g=bitpack_encode(nx.to_numpy_array(new_G,dtype=int,nodelist=sorted_nodes,weight="weight"),d)
        
        relabeled_graphs.append(g)

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

# def generate_graphs(filename, n, d):
    
#     with open(filename, "r") as file:
#         graph6_lines = [line.strip() for line in file if line.strip()]

#     iso_graphs = [bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines]
#     ts=time.time()
    
#     #draw_graphs(iso_graphs,n,d)
#     weighted_iso_graphs=[]
#     i=1
#     for ig in iso_graphs:
#         temp_ts=time.time()
#         weighted_iso_graphs.append(ig) 
#         weighted_graphs = set(generate_weighted(bitpack_decode(ig,n,d),d))
    
#         while len(weighted_graphs)>0:
#             g=weighted_graphs.pop()
#             weighted_iso_graphs.append(g)
#             perms=generate_permuted(bitpack_decode(g,n,d),d)
#             weighted_graphs.difference_update(perms)

#         print(f"{i}: {time.time()-temp_ts}")
#         i+=1
   
#     #print(f"end: {time.time()-ts}")
#     return weighted_iso_graphs

async def process_graph(ig,n,d,output_dir, index):
        temp_ts=time.time()
        results=[]
        results.append(ig) 
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
    #draw_graphs(iso_graphs,n,d)
    #i=1
    
        #print(f"{i}: {time.time()-temp_ts}")
        #i+=1
    await asyncio.gather(*[process_graph(ig,n,d,output_dir,i) for i,ig in enumerate(iso_graphs)])
    print(f"end: {time.time()-ts}")

def total_weight(graph):
    bit_length=int(graph).bit_length()
    pair_sum = []
    for i in range(0, bit_length, 2):
        # Extract the current pair of bits
        pair = (graph >> i) & 0b11  # Mask the last two bits
        pair_sum.append(pair)  # Count the number of set bits in the pair

    return sum(pair_sum)  # Reverse to maintain left-to-right order

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
        if current.bit_count()<=min_graph.bit_count():
                if total_weight(current)<=total_weight(min_graph):
                    if current<min_graph:
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


def find_orbit_isomorphic(graphs,start_graph, n, d):
    queue = deque([start_graph])  # BFS queue
    visited = set()  # Set of visited graphs
    visited.add(start_graph)
    
    # Define all possible scaling and complementation factors
    scaling_factors = range(2, d)  # Local scaling factors (mod d)
    complementing_factors = range(1, d)
    min_graph=start_graph

    while queue:
        current_level = list(queue)
        queue.clear()
        temp_level=set()
        for current in current_level:
            #print("current")
            #draw_graph(current,n,d)
            
            current_matrix = bitpack_decode(current, n, d)
            if current.bit_count()<=min_graph.bit_count():
                if total_weight(current)<=total_weight(min_graph):
                    if current<min_graph:
                        min_graph=current
            
            for v in range(n):
                # Local scaling
                for k in scaling_factors:
                    scaled_matrix = current_matrix.copy()
                    local_scaling(scaled_matrix, v, k, d)
                    encoded_scaled = bitpack_encode(scaled_matrix, d)
                    if encoded_scaled not in visited:
                        temp_level.add(encoded_scaled)
                
                # Local complementation
                for k in complementing_factors:
                    complemented_matrix = current_matrix.copy()
                    local_complementation(complemented_matrix, v, k, d)
                    encoded_complemented = bitpack_encode(complemented_matrix, d)
                    if encoded_complemented not in visited:
                        temp_level.add(encoded_complemented)
        while temp_level:
            g=temp_level.pop()
            #print(nx.from_numpy_array(bitpack_decode(g,n,d)).edges.data())
            #print("g")
            #draw_graph(g,n,d)
            perms=set(generate_permuted(bitpack_decode(g,n,d),d))
            inter_visited=bool(perms.intersection(visited))
            if not inter_visited:
                visited.add(g)
                queue.append(g)
                graphs.difference_update(perms)
            temp_level.difference_update(perms)
                    
    return visited, min_graph


def find_orbit_isomorphic_no_global_filter(start_graph, n, d):
    queue = deque([start_graph])  # BFS queue
    visited = set()  # Set of visited graphs
    visited.add(start_graph)
    
    # Define all possible scaling and complementation factors
    scaling_factors = range(2, d)  # Local scaling factors (mod d)
    complementing_factors = range(1, d)
    min_graph=start_graph
    ts=time.time()
    while queue:
        print(len(queue))
        print(f"{time.time()-ts}")
        current_level = list(queue)
        queue.clear()
        temp_level=set()
        for current in current_level:
            #print("current")
            #draw_graph(current,n,d)
            
            current_matrix = bitpack_decode(current, n, d)
            if current.bit_count()<=min_graph.bit_count():
                if total_weight(current)<=total_weight(min_graph):
                    if current<min_graph:
                        min_graph=current
            
            for v in range(n):
                # Local scaling
                for k in scaling_factors:
                    scaled_matrix = current_matrix.copy()
                    local_scaling(scaled_matrix, v, k, d)
                    encoded_scaled = bitpack_encode(scaled_matrix, d)
                    if encoded_scaled not in visited:
                        temp_level.add(encoded_scaled)
                
                # Local complementation
                for k in complementing_factors:
                    complemented_matrix = current_matrix.copy()
                    local_complementation(complemented_matrix, v, k, d)
                    encoded_complemented = bitpack_encode(complemented_matrix, d)
                    if encoded_complemented not in visited:
                        temp_level.add(encoded_complemented)
        while temp_level:
            g=temp_level.pop()
            #print(nx.from_numpy_array(bitpack_decode(g,n,d)).edges.data())
            #print("g")
            #draw_graph(g,n,d)
            #perms=set(generate_permuted(bitpack_decode(g,n,d),d))
            #inter_visited=bool(perms.intersection(visited))
            if g not in temp_level:
                visited.add(g)
                queue.append(g)
            #temp_level.difference_update(perms)
                    
    return visited, min_graph



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

def orbit_search_isomorphic_from_file(n,d):
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
    print(len(graphs))
    print(f"Start:{time.time()-ts}")
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


# Convert each line from graph6 to a NetworkX graph
n=7
d=3
# o,f_o=orbit_search(f"d3n{n}.txt",n,d)
# print(len(o))
# print("---------------------")
#o=(orbit_search_isomorphic(f"d3n{n}.txt",n,d))
# print(len(o))
if __name__ == "__main__":
    multiprocessing.freeze_support()
    asyncio.run(generate_graphs(f"d3n{n}.txt", n, d))
# draw_graphs(o,n,d, cols=6)
#print(len(generate_graphs(f"d3n{n}.txt",n,d)))

#orbit_search_isomorphic_from_file(n,d)
#v,m=find_orbit_isomorphic_no_global_filter(2227370722648,n,d)
#print(len(v))
#draw_graph(m,n,d)
"""
ts=time.time()
print(time.time()-ts)
gg=generate_graphs("d3n4.txt",n,d)
print(len(gg))
print(time.time()-ts)
ts=time.time()
print(time.time()-ts)
cg=create_graphs("d3n4.txt",n,d)
print(len(cg))
print(time.time()-ts)
"""
#for g in o:
#    draw_graph(g,n,d)
#for i in range(12):
   #draw_graph(gg[i],n,d)
   #draw_graph(cg[i],n,d)