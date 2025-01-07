import taichi as ti
import numpy as np

ti.init(arch=ti.gpu)
n=4
# Adjacency matrix (n x n)
W = ti.ndarray(dtype=ti.i32, shape=(n, n))
for i in range(n):
    for j in range(n):
        W[i,j]=2

W[0,3]=0
W[3,0]=0
W[1,2]=1
W[2,1]=1
print(W.to_numpy())
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
