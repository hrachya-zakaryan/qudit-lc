import networkx as nx
from collections import deque
<<<<<<< Updated upstream
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph, draw_graphs, draw_representatives
=======
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph, draw_graphs
>>>>>>> Stashed changes
import numpy as np
import time
import multiprocessing
import os
from matplotlib import pyplot as plt
import subprocess
import matplotlib.lines as mlines
import sympy as sp
from itertools import combinations

def total_weight(graph):
    bit_length=int(graph).bit_length()
    pair_sum = []
    for i in range(0, bit_length, 2):
        pair = (graph >> i) & 0b11 
        pair_sum.append(pair)  
    return sum(pair_sum)


<<<<<<< Updated upstream
def call_generate_graphs(n,d, encoded_value,path,index):
    """
    Calls the C code to generate all non-isomorphic weighted graphs
    for a given graph structure and writes the bitpack encoded graphs
    to a single file based on the index.

    Args:
        n(int): Number of vertices
        d(int): Local dimension
        encoded_value(int): Bitpacked graph for the simple graph
        path(str): Path to the C file
        index(int): The index of the simple graph
    """

    result = subprocess.run([path, str(n), str(d), str(encoded_value)], capture_output=True, text=True) #Calls the C code
    output_dir = os.path.join(os.getcwd(),f"d{d}_c",str(n))
    graphs = {int(line) for line in result.stdout.splitlines()} #Parse the C output
=======
def call_generate_graphs(n,d, encoded_value,path,index,ts):
    #print("start:",index,":",time.time()-ts)
    result = subprocess.run(
        [path, str(n), str(d), str(encoded_value)], 
        capture_output=True, text=True
    )
    print("end", index,":",time.time()-ts)
    output_dir = os.path.join(os.getcwd(),f"d{d}_c",str(n))
    graphs = {int(line) for line in result.stdout.splitlines()}
>>>>>>> Stashed changes
    output_file = os.path.join(output_dir, f"{index}.txt")
    with open(output_file, "w") as f:
        for graph in graphs:
            f.write(str(graph) + "\n")
<<<<<<< Updated upstream

def generate_graphs_c(n,d,path, num_workers=multiprocessing.cpu_count()):
    """
    Generates all the weighted non-isomorphic graphs for a given n and d.
    Uses multiprocessing.

    Args:
        n(int): Number of vertices
        d(int): Local dimension
        path(str): Path to the C file
        num_workers(int): Number of CPU threads to use. Uses all of them by default
    """

    filename=f"d3n{n}.txt"
    with open(filename, "r") as file:
        graph6_lines = [line.strip() for line in file if line.strip()]
    iso_graphs = [bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes(line.encode()),dtype=int),d) for line in graph6_lines] #Read the non-isomorphic simple graphs and convert from graph6 to bitpack encoding
    ts=time.time()

    output_dir = os.path.join(os.getcwd(),f"d{d}_c",str(n))
    os.makedirs(output_dir, exist_ok=True)
    with multiprocessing.Pool(num_workers) as pool:
        index=0
        while iso_graphs:
            batch_size = min(len(iso_graphs), num_workers)
            batch = [iso_graphs.pop(0) for _ in range(batch_size)] #Create a batch
            print(f"{len(iso_graphs)}:{time.time()-ts}")
            pool.starmap(call_generate_graphs, [(n, d, g, path,index*num_workers+batch.index(g),ts) for g in batch]) #Run call_generate_graphs in parallel for the graphs in batch
            index+=1
    


def call_complementation_layer(n, d, encoded_value,path):
    """
    Calls the C code that scales and complements a single layer.
    Returns a set of the minimal isomorphic graph encodings.

    Args:
        n(int): Number of vertices
        d(int): Local dimension
        encoded_value(int): Bitpacked graph
        path(str): Path to the C file
    """
    result = subprocess.run([path, str(n), str(d), str(encoded_value)], capture_output=True, text=True) #Calls the C file
    return {int(line) for line in result.stdout.splitlines()}


def orbit_atlas_c(path,n,d, num_workers=multiprocessing.cpu_count()):
    """
    Creates an atlas of local scalings and local complementatinos for
    all the graphs of local dimension d and number of vertices n.
    
    Args:
        path(str): Path to the C file
        n(int): Number of vertices
        d(int): Local dimension
        num_workers(int): Number of CPU threads to use. Uses all of them by default
    """
=======

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
>>>>>>> Stashed changes
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
    G.add_nodes_from(graphs) #Initialize the atlas graph with all the graphs as vertices.

    with multiprocessing.Pool(num_workers) as pool:
        while graphs:
            batch_size = min(len(graphs), num_workers)
            batch = [graphs.pop() for _ in range(batch_size)] #Create a batch
            print(f"{len(graphs)}:{time.time()-ts}")
<<<<<<< Updated upstream
            complements = pool.starmap(call_complementation_layer, [(n, d, g, path) for g in batch]) #Run call_complementation_layer in parallel for the graphs in batch
=======
            complements = pool.starmap(call_complementation_layer, [(n, d, g, path) for g in batch])
>>>>>>> Stashed changes
            for i in range(batch_size):
                for complement in complements[i]:
                    G.add_edge(batch[i],complement) #Add the edges obtained from local scaling and local complementing
    return G


def separate_orbits(G):
    """
    Separates the atlas into orbits.

    Args:
        G(nx.Graph): Atlas graph
    """
    components = list(nx.connected_components(G)) #Separates disjoint parts of the atlas
    subgraphs = [G.subgraph(nodes).copy() for nodes in components]
    return subgraphs

def find_representatives(n,d):
    represenatives=[]
    for i in range(3,n+1):
        directory=f"orbits_d{d}_n{i}_separated"
        for j in range(len(os.listdir(directory))):
            file_path=os.path.join(directory,f"orbit_{j}")
            G=nx.read_edgelist(file_path,data=False)
            graphs=list(G.nodes)
            represenatives.append([int(min(graphs,key=lambda x: (int(x).bit_count(), total_weight(int(x)),x))),(i,j)])
    return represenatives
            

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
        adjacency_matrix = bitpack_decode(int(node), n, d)
        
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

def call_measure(n,d,directory,i,strategy):
<<<<<<< Updated upstream
    """
    Calculates the schmidt measure for orbit i.
    
    Args:
        n(int): Number of vertices
        d(int): Local dimension
        directory(str): Directory of orbit files
        i(int): Index of the orbit
        strategy(str): Coloring strategy for greedy_color
    """

    current_measure=0
    orbit=nx.read_edgelist(directory+f"/orbit_{i}", data=False) #Import orbit
    exists_two_colorable=False
    for g in orbit.nodes: #Check if there is a two colorable graph in the orbit
=======
    current_measure=0
    orbit=nx.read_edgelist(directory+f"/orbit_{i}", data=False)
    exists_two_colorable=False
    for g in orbit.nodes:
>>>>>>> Stashed changes
        adj_matrix=bitpack_decode(int(g),n,d)
        nx_g=nx.from_numpy_array(adj_matrix)
        coloring = nx.coloring.greedy_color(nx_g, strategy=strategy)
        if max(coloring.values())==1:
            exists_two_colorable=True
            break
    max_max_rank=0
    found=False
    for g in orbit.nodes:
        result=calculate_schmidt_measure(g,exists_two_colorable, max_max_rank,n,d, strategy)
        if type(result[1])==int:
            max_max_rank=result[1]
        if result[2]==True:
            current_measure=result[0]
            if type(result[1])==bool:
                found=True
                break
    if type(current_measure) == int and current_measure==0:
        current_measure=[0,0]
    if found==False:
        temp_graphs=list(orbit.nodes)
        for j in range(1,n-2): 
            graphs=temp_graphs.copy()
            temp_graphs=[]
            for g in graphs:
                adj_matrix=bitpack_decode(int(g),n-j+1,d)
                nx_g=nx.from_numpy_array(adj_matrix)
                for v in nx_g.nodes:
                    temp_nx_g=nx_g.copy()
                    temp_nx_g.remove_node(v)
                    if nx.is_connected(temp_nx_g):
                        temp_graphs.append(bitpack_encode(nx.to_numpy_array(temp_nx_g,dtype=int),d))
                        mapping = dict(zip(temp_nx_g, range(n-j)))
                        temp_nx_g=nx.relabel_nodes(temp_nx_g, mapping)
                        result=sub_calculate_schmidt_measure(temp_nx_g,n-j,d, strategy)
                        if result !=0 and (current_measure[1]==0 or result+j<current_measure[1]):
                            current_measure[1]=result+j
                            #print(j)
                        if result !=0 and (result>current_measure[0]):
                            current_measure[0]=result
                        if result !=0 and current_measure[0]==current_measure[1]:
                            break
                if current_measure[0]==current_measure[1]:
                    break
            if current_measure[1]!=0 and current_measure[0]==current_measure[1]:
                current_measure=current_measure[1]
                found = True
                break
        print(f"{i}:{current_measure}:Color")
    if found == False:
        for g in orbit.nodes:
            adj_matrix=bitpack_decode(int(g),n,d)
            nx_g=nx.from_numpy_array(adj_matrix)
            result =  all_bipartitions_measure(adj_matrix,n,d)
            if result != 0  and result > current_measure[0]:
                current_measure[0]=result
<<<<<<< Updated upstream
            if current_measure[1]!=0 and current_measure[0]==current_measure[1]:
                current_measure=current_measure[1]
                found = True
                break
=======
>>>>>>> Stashed changes
        print(f"{i}:{current_measure}")
    else:
        print(f"{i}:{current_measure}:Early")
    return current_measure


def all_bipartitions_measure(adj_matrix,n,d):
    nx_g=nx.from_numpy_array(adj_matrix)
    max_rank=0
    for i in range(1, n // 2 + 1):
        for subset in combinations(list(nx_g.nodes), i):
            submatrix=adj_matrix[np.ix_(sorted(subset),sorted(set(nx_g.nodes)-set(subset)))]
            rank = rank_over_finite_field(submatrix,d)
            if rank>max_rank:
                max_rank=rank
    return max_rank



def orbit_schmidt_measure(n,d, strategy):
    
    directory=f"orbits_d{d}_n{n}_separated"
    num_orbits=len(list(os.listdir(directory)))
    schmidt_measures=[]
    num_workers = multiprocessing.cpu_count()
    orbits=list(range(num_orbits))
    with multiprocessing.Pool(num_workers) as pool:
        while orbits:
            batch_size = min(len(orbits), num_workers)
            batch = [orbits.pop(0) for _ in range(batch_size)]
            schmidt_measures.extend(pool.starmap(call_measure, [(n, d, directory, i, strategy) for i in batch]))
            print(len(orbits))                
    return schmidt_measures


def sub_calculate_schmidt_measure(nx_g,n,d, strategy):
    adj_matrix=nx.to_numpy_array(nx_g,dtype=int)
    coloring = nx.coloring.greedy_color(nx_g, strategy=strategy)
    color_classes = {c: [v for v in coloring if coloring[v] == c] for c in set(coloring.values())}
    if max(coloring.values()) == 1 and rank_over_finite_field(adj_matrix,d)==n:
        return n//2
        
    elif max(coloring.values()) == 1:
        max_rank = 0
        colors=list(color_classes.keys())
        min1=10
        for k in range(len(colors)):
            A = color_classes[colors[k]]  # One part of the bipartition
            if len(A)<min1:
                min1=len(A)
            B = [v for j in range(len(colors)) if j != k for v in color_classes[colors[j]]]  # Other part
            submatrix=adj_matrix[np.ix_(A, B)]
            rank = rank_over_finite_field(submatrix,d)
            max_rank = max(max_rank, rank)
        if max_rank==min1:
            return max_rank
        elif max_rank==n//2:
            return max_rank
    elif max(coloring.values()) == 2:
        max_rank = 0
        colors=list(color_classes.keys())
        min1=10
        min2=10
        for k in range(len(colors)):
            A = color_classes[colors[k]]  # One part of the bipartition
            if len(A)<min1:
                min1=len(A)
            elif len(A)<min2:
                min2=len(A)
            B = [v for j in range(len(colors)) if j != k for v in color_classes[colors[j]]]  # Other part
            submatrix=adj_matrix[np.ix_(A, B)]
            rank = rank_over_finite_field(submatrix,d)
            max_rank = max(max_rank, rank)
        if max_rank==min1+min2:
            return max_rank
    elif max(coloring.values()) == 1:
        return n//2
    return 0

def calculate_schmidt_measure(g, exists_two_colorable, max_max_rank,n,d, strategy):
    adj_matrix=bitpack_decode(int(g),n,d)
    nx_g=nx.from_numpy_array(adj_matrix)
    coloring = nx.coloring.greedy_color(nx_g, strategy=strategy)
    color_classes = {c: [v for v in coloring if coloring[v] == c] for c in set(coloring.values())}
    if max(coloring.values()) == 1 and rank_over_finite_field(adj_matrix,d)==n:
        return n//2, True, True
        
    elif max(coloring.values()) == 1:
        max_rank = 0
        colors=list(color_classes.keys())
        min1=10
        for k in range(len(colors)):
            A = color_classes[colors[k]]  # One part of the bipartition
            if len(A)<min1:
                min1=len(A)
            B = [v for j in range(len(colors)) if j != k for v in color_classes[colors[j]]]  # Other part
            submatrix=adj_matrix[np.ix_(A, B)]
            rank = rank_over_finite_field(submatrix,d)
            max_rank = max(max_rank, rank)
        if max_rank==min1:
            return max_rank, True, True
        elif max_rank==n//2:
            return max_rank, True, True
        elif  max_max_rank<=max_rank:
            max_max_rank=max_rank
            #print(f"{max_rank}<{n//2}")
            # draw_graph(int(g),n,d)
            # plot_graph(orbit,n,d)
            
            return [max_rank,n//2], max_max_rank, True
    elif max(coloring.values()) == 2:
        max_rank = 0
        colors=list(color_classes.keys())
        min1=10
        min2=10
        for k in range(len(colors)):
            A = color_classes[colors[k]]  # One part of the bipartition
            if len(A)<min1:
                min1=len(A)
            elif len(A)<min2:
                min2=len(A)
            B = [v for j in range(len(colors)) if j != k for v in color_classes[colors[j]]]  # Other part
            submatrix=adj_matrix[np.ix_(A, B)]
            rank = rank_over_finite_field(submatrix,d)
            max_rank = max(max_rank, rank)
        if max_rank==min1+min2:
            return max_rank, True, True

        elif max_max_rank<max_rank:
            max_max_rank=max_rank
            if exists_two_colorable:
                #print(f"{max_rank}<{n//2}")
                return [max_rank,n//2], max_max_rank, True
            else:
                #print(f"{max_rank}<")
                return [max_rank,0], max_max_rank, True
            # draw_gfrom itertools import combinationsaph(orbit,n,d)
    return 0, False, False
def mod_inv(a, d):
    """Compute modular inverse of a mod d (assuming d is prime)."""
    return pow(a, -1, d)  # Python 3.8+ supports pow(a, -1, d) for modular inverse

def rank_over_finite_field(matrix, d):
    """Compute the rank of a matrix over the finite field F_d."""
    mat = sp.Matrix(matrix).applyfunc(lambda x: x % d)  # Reduce elements mod d
    rref_matrix, pivot_cols = mat.rref(iszerofunc=lambda x: x % d == 0)
    return len(pivot_cols)


<<<<<<< Updated upstream
n=7
d=3
#directory=f"orbits_d{d}_n{n}_separated"
#call_measure(n,d,directory,46,'independent_set')
#print(orbit_schmidt_measure(n,d,'independent_set'))
#r=find_representatives(n,d)
r=[[5, (3, 0)], [21, (4, 0)], [81, (4, 1)], [344, (4, 2)], [20816, (5, 0)], [1368, (5, 1)], [1348, (5, 2)], [337, (5, 3)], [85, (5, 4)], [279888, (6, 0)], [344452, (6, 1)], [278789, (6, 2)], [1361, (6, 3)], [1335569, (6, 4)], [1380736, (6, 5)], [18433360, (6, 6)], [1380672, (6, 7)], [21316944, (6, 8)], [279824, (6, 9)], [5444, (6, 10)], [21784, (6, 11)], [345488, (6, 12)], [279826, (6, 13)], [341, (6, 14)], [5464, (6, 15)], [21780, (6, 16)], [279814, (6, 17)], [279813, (6, 18)], [5441, (6, 19)], [1401232, (6, 20)], [1146160384, (7, 0)], [21251136, (7, 1)], [4265217, (7, 2)], [1141904656, (7, 3)], [1163150656, (7, 4)], [5453714704, (7, 5)], [4265232, (7, 6)], [1146230080, (7, 7)], [1146160640, (7, 8)], [1146160416, (7, 9)], [22299968, (7, 10)], [5470491920, (7, 11)], [22300224, (7, 12)], [1167345984, (7, 13)], [1146423620, (7, 14)], [89212192, (7, 15)], [1145439297, (7, 16)], [1146177856, (7, 17)], [21251392, (7, 18)], [1146119488, (7, 19)], [4527184, (7, 20)], [21251137, (7, 21)], [21251168, (7, 22)], [21251169, (7, 23)], [22287968, (7, 24)], [4527488, (7, 25)], [25510992, (7, 26)], [21251456, (7, 27)], [21333252, (7, 28)], [1145373716, (7, 29)], [4527504, (7, 30)], [5470492752, (7, 31)], [1141921088, (7, 32)], [1099240784, (7, 33)], [1141920840, (7, 34)], [87300, (7, 35)], [1146423616, (7, 36)], [5655258688, (7, 37)], [1099240528, (7, 38)], [1146115392, (7, 39)], [4265234, (7, 40)], [1145378324, (7, 41)], [1145373720, (7, 42)], [1145382144, (7, 43)], [5655238032, (7, 44)], [21251138, (7, 45)], [1145378128, (7, 46)], [89212176, (7, 47)], [21333392, (7, 48)], [1365, (7, 49)], [4527172, (7, 50)], [4265296, (7, 51)], [21333248, (7, 52)], [21071121, (7, 53)], [21828, (7, 54)], [21825, (7, 55)], [4523076, (7, 56)], [21848, (7, 57)], [4523396, (7, 58)], [22283872, (7, 59)], [21251090, (7, 60)], [4265221, (7, 61)], [4265222, (7, 62)], [4261125, (7, 63)], [4527176, (7, 64)], [21251089, (7, 65)], [21071105, (7, 66)], [21251094, (7, 67)], [5457, (7, 68)], [349272, (7, 69)], [87316, (7, 70)], [87320, (7, 71)], [21251088, (7, 72)]]
draw_representatives(r,d,8)
# generate_graphs_c(n,d,"c/generate_graphs")
# g=orbit_atlas_c("c/complement",n,d)
# if __name__ == "__main__":
#     multiprocessing.freeze_support()
#     g=orbit_atlas_c(n,d)
# g=nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_6", data=False)
# print("Vertices:",g.number_of_nodes())
# print("Edges:", g.number_of_edges())
# encoded_value=bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes("EU~w".encode()),dtype=int),d)
# call_generate_graphs(n,d, encoded_value,"c/generate_graphs",96,time.time())
# orbits=separate_orbits(g)
=======
n=6
d=3
>>>>>>> Stashed changes

# directory = f"orbits_d{d}_n{n}_separated"
# os.makedirs(directory, exist_ok=True)
# for i in range(len(orbits)):
#     nx.write_edgelist(orbits[i],f"{directory}/orbit_{i}", data=False)

<<<<<<< Updated upstream
=======
# generate_graphs_c(n,d,"c/generate_graphs")
# g=orbit_atlas_c("c/complement",n,d)
if __name__ == "__main__":
    multiprocessing.freeze_support()
    # print(orbit_schmidt_measure(n,d,'independent_set'))
    sm_in_og = []
    for cnt in range(3,n+1):
        print(f"------------------------------------Start calculating for n={cnt}------------------------------------")
        sm_in_og.append(orbit_schmidt_measure(cnt,d,'independent_set'))
    
    
    with open("og_data/join_sm_in_og.txt", "w") as f:
        f.write(str(sm_in_og))
    # g=orbit_atlas_c(n,d)
# g=nx.read_edgelist(f"orbits_d{d}_n{n}_separated/orbit_6", data=False)
# print("Vertices:",g.number_of_nodes())
# print("Edges:", g.number_of_edges())
# encoded_value=bitpack_encode(nx.to_numpy_array(nx.from_graph6_bytes("EU~w".encode()),dtype=int),d)
# call_generate_graphs(n,d, encoded_value,"c/generate_graphs",96,time.time())
# orbits=separate_orbits(g)

# directory = f"orbits_d{d}_n{n}_separated"
# os.makedirs(directory, exist_ok=True)
# for i in range(len(orbits)):
#     nx.write_edgelist(orbits[i],f"{directory}/orbit_{i}", data=False)

>>>>>>> Stashed changes
# # print(time.time()-ts)
# for o in orbits:
#   plot_graph(o,n,d)
# nx.draw(g)
# plt.show()
#orbit_atlas_c(n,d)
# strategies=['largest_first','random_sequential','independent_set','connected_sequential_bfs','connected_sequential_dfs','saturation_largest_first']
# for strategy in strategies:
#     print(strategy)
#     print(orbit_schmidt_measure(n,d,strategy))
#g=nx.from_numpy_array(bitpack_decode(703, n, d))
# nx.draw(g)
<<<<<<< Updated upstream
# plt.show()graph
=======
# plt.show()
>>>>>>> Stashed changes
# g.remove_node(1)
# nx.draw(g)
# plt.show()
#calculate_schmidt_measure(703,False,0,n,d)