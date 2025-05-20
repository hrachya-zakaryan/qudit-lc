import networkx as nx
from collections import deque
from auxilary_unparallelized import local_complementation, local_scaling, bitpack_decode, bitpack_encode, draw_graph, draw_graphs, draw_representatives
import numpy as np
import time
import multiprocessing
import os
from matplotlib import pyplot as plt
import subprocess
import matplotlib.lines as mlines
import sympy as sp
from itertools import combinations
import shutil
import matplotlib.colors as mcolors

def total_weight(graph):
    bit_length=int(graph).bit_length()
    pair_sum = []
    for i in range(0, bit_length, 2):
        pair = (graph >> i) & 0b11 
        pair_sum.append(pair)  
    return sum(pair_sum)


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
    output_file = os.path.join(output_dir, f"{index}.txt")
    with open(output_file, "w") as f:
        for graph in graphs:
            f.write(str(graph) + "\n")

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
            complements = pool.starmap(call_complementation_layer, [(n, d, g, path) for g in batch]) #Run call_complementation_layer in parallel for the graphs in batch
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
        directory=f"orbits_d{d}_n{i}_separated_sorted"
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
    edge_colors = ['black' if adjacency_matrix[u, v] == 1 else 'red' for u, v in edges]

    # Draw subgraph edges with the appropriate color
    nx.draw_networkx_edges(subG, pos, ax=ax, edge_color=edge_colors, alpha=1, width=1.3)

    # Draw subgraph nodes (light blue fill, black outline)
    nx.draw_networkx_nodes(subG, pos, ax=ax, node_size=23, node_color='blue', edgecolors='black', linewidths=0.6)

    # Draw edge weights (small font)
    #nx.draw_networkx_edge_labels(subG, pos, edge_labels=edge_weights, ax=ax, font_size=4, font_color='black')

def plot_graph(G, n, d, seed, orbit):
    """
    Plot a graph G using NetworkX with each node as a circle containing a small graph.
    This version draws into an existing matplotlib Axes (ax).
    """
    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(G, seed=seed, k=1/np.sqrt(len(G.nodes)))  # Layout for main graph

    # Draw main graph edges
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='blue', alpha=0.7, width=0.8)

    # Draw main graph nodes as circles with subgraphs inside
    for node in G.nodes():
        x, y = pos[node]
        circle = plt.Circle((x, y), 0.1, color='white', ec='black', lw=1.2)
        ax.add_patch(circle)

        adjacency_matrix = bitpack_decode(int(node), n, d)
        draw_subgraph_inside_circle(ax, (x, y), 0.08, adjacency_matrix)

    # Set axis limits and aspect ratio
    ax.set_xlim(min(x for x, y in pos.values()) - 0.2, max(x for x, y in pos.values()) + 0.2)
    ax.set_ylim(min(y for x, y in pos.values()) - 0.2, max(y for x, y in pos.values()) + 0.2)
    ax.set_aspect('equal')
    ax.axis('off')

    # Legend (optional to skip in each subplot if too cluttered)
    legend_elements = [
        mlines.Line2D([], [], color='black', lw=1.3, label='Edge Weight = 1'),
        mlines.Line2D([], [], color='red', lw=1.3, label='Edge Weight = 2')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=7, frameon=True)
    plt.tight_layout()
    plt.savefig(f"n{n}d{d}_o{orbit}_s{seed}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

def call_measure(n,d,directory,i,strategy):
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
    
    directory=f"orbits_d{d}_n{n}_separated_sorted"
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


def orbit_value(orbit,directory):
    nodes=nx.read_edgelist(os.path.join(directory,orbit)).nodes
    rep=int(min(nodes,key=lambda x: (int(x).bit_count(), total_weight(int(x)),x)))
    return (int(rep).bit_count(), total_weight(int(rep)),int(rep))


def sort_orbits(d):
    for n in range(3,8):
        directory=f"orbits_d{d}_n{n}_separated"
        sorted_orbits=sorted(os.listdir(directory),key= lambda o: orbit_value(o,directory))
        new_directory=f"orbits_d{d}_n{n}_separated_sorted"
        os.makedirs(new_directory, exist_ok=False)
        for new_index, original_filename in enumerate(sorted_orbits):
            src_path = os.path.join(directory, original_filename)
            dst_filename = f'orbit_{new_index}'
            dst_path = os.path.join(new_directory, dst_filename)
            shutil.copy2(src_path, dst_path)  # copy2 preserves metadata


def plot_distance_matrix(G,filename):
    """
    Compute and plot the shortest path distance matrix of a graph G using matplotlib.
    The matrix contains only integer distances.

    Parameters:
        G (networkx.Graph): Input graph.
    """
    # Get the list of nodes
    nodes = sorted(list(G.nodes()),key=lambda x: (int(x).bit_count(), total_weight(int(x)),int(x)))
    n = len(nodes)
    
    edge_changes=[]
    temp_edges=0
    for i in range(n):
        if temp_edges<int(nodes[i]).bit_count():
            edge_changes.append((i,int(nodes[i]).bit_count()))
        temp_edges=int(nodes[i]).bit_count()
    # Initialize an empty distance matrix
    dist_matrix = np.zeros((n, n), dtype=int)
    
    # Compute the shortest path lengths
    for i, node_i in enumerate(nodes):
        lengths = nx.shortest_path_length(G, source=node_i)
        for j, node_j in enumerate(nodes):
            if node_j in lengths:
                dist_matrix[i, j] = int(lengths[node_j])  # Ensure integers
            else:
                dist_matrix[i, j] = np.inf  # Use infinity for unreachable nodes

    # Plot the matrix using imshow
    fig, ax = plt.subplots(figsize=(8, 6))

    cmap = plt.cm.PuBuGn  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    
    cmap=mcolors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.arange(0, np.max(dist_matrix)+2)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    cax = ax.imshow(dist_matrix, cmap=cmap, interpolation="nearest",norm=norm)
    cbar = fig.colorbar(cax, label="Shortest Path Distance")
    cbar.set_ticks(np.arange(0.5, np.max(dist_matrix)+1, 1),labels=np.arange(0, np.max(dist_matrix)+1, 1))
    cbar.ax.tick_params(which='both', length=0)
    ax.set_title("Distance Matrix")
    ax.set_xlabel("Number of Edges in the Graph")
    ax.set_ylabel("Encoded Graph")
    ax.set_xticks([])
    yticks = list(zip(*edge_changes))[0]
    ax.set_yticks([tick - 0.5 for tick in yticks])
    ylabels=[]
    for i in yticks:
        ylabels.append(int(nodes[i]))
    #ax.set_xticklabels(nodes, rotation=90)
    ax.set_yticklabels(ylabels)
    for idx in range(1,len(edge_changes)):
            ax.axhline(y=edge_changes[idx][0]-0.5, color='black', linestyle='solid', linewidth=1)
            ax.axvline(x=edge_changes[idx][0]-0.5, color='black', linestyle='solid', linewidth=1)
    xticks=[]
    for idx in range(len(edge_changes)-1):
            xticks.append((edge_changes[idx+1][0] + edge_changes[idx][0]) / 2-0.5)
    xticks.append((len(nodes) + edge_changes[-1][0]) / 2-0.5)
    ax.set_xticks(xticks, labels=list(zip(*edge_changes))[1])
    ax.tick_params(axis='x', which='both', length=0)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)

n=7
d=3
n=4
g= nx.read_edgelist(f"orbits_d{d}_n{n}_separated_sorted/orbit_2", data=False)
plot_distance_matrix(g,'n4o2.png')
#plot_graph(g, n, d, 49, 2)

#directory=f"orbits_d{d}_n{n}_separated"
#call_measure(n,d,directory,46,'independent_set')
#print(orbit_schmidt_measure(n,d,'independent_set'))
# r=find_representatives(n,d)
# print(r)
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

# directory = f"orbits_d{d}_n{n}_separated"
# os.makedirs(directory, exist_ok=True)
# for i in range(len(orbits)):
#     nx.write_edgelist(orbits[i],f"{directory}/orbit_{i}", data=False)


# directory = f"orbits_d{d}_n{n}_separated"
# os.makedirs(directory, exist_ok=True)
# for i in range(len(orbits)):
#     nx.write_edgelist(orbits[i],f"{directory}/orbit_{i}", data=False)

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
# plt.show()
# g.remove_node(1)
# nx.draw(g)
# plt.show()
#calculate_schmidt_measure(703,False,0,n,d)