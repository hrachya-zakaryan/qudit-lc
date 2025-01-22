import taichi as ti
import numpy as np
from itertools import product

def bitpack_encode(matrix, d):
    n = matrix.shape[0]
    bit_length = d.bit_length()  # Number of bits per weight
    packed = 0
    shift = 0
    
    # Iterate through the upper triangular part (excluding diagonal)
    for i in range(n):
        for j in range(i + 1, n):
            weight = matrix[i, j]
            packed |= (weight << shift)  # Shift and add weight to packed
            shift += bit_length  # Update shift for the next weight
    
    return packed

def bitpack_decode(packed, n, d):
    bit_length = d.bit_length()
    matrix = np.zeros((n, n), dtype=int)
    shift = 0
    
    # Decode weights from packed integer
    for i in range(n):
        for j in range(i + 1, n):
            weight = (packed >> shift) & ((1 << bit_length) - 1)  # Extract bits
            matrix[i, j] = weight
            matrix[j, i] = weight  # Reflect symmetry
            shift += bit_length  # Move to the next weight
    
    return matrix



# Local scaling
@ti.func
def local_scaling(G: ti.types.ndarray(dtype=ti.i32, ndim=2), v: int, k: int, d: int):
    n = G.shape[0]
    for u in range(n):
        if u != v:
            G[v, u] = (k * G[v, u]) % d
            G[u, v] = G[v, u]  # Symmetric

# Local complementation

@ti.func
def local_complementation(G: ti.types.ndarray(dtype=ti.i32, ndim=2), v: int, k: int, d: int):
    n = G.shape[0]
    for u in range(n):
        for w in range(n):
            if u!=w and u != v and w != v and G[v, u] != 0 and G[v, w] != 0:
                G[u, w] = (G[u, w] + k * G[v, u] * G[v, w]) % d


@ti.kernel
def create_graphs(n: int,d: int):
    num_entries = n * (n - 1) // 2
    counter=0
    for combination in product(range(d), repeat=num_entries):
            # Reconstruct the symmetric matrix
        matrix = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                matrix[i, j] = combination[idx]
                matrix[j, i] = combination[idx]
                idx += 1
        if not np.any(np.all(matrix == 0, axis=0)) and not np.any(np.all(matrix == 0, axis=1)):
            print(matrix)
            counter+=1
    print(counter)

ti.init(arch=ti.gpu)
create_graphs(3,3)
