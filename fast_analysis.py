from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path


# ----------------------------------------------------------------- loading

@dataclass
class OG:
    adj: csr_matrix           # n x n, symmetric, no self-loops
    labels: np.ndarray        # dense index -> original graph code
    n_selfloops: int          # fixed-point edges (g, g) seen in the file

    @property
    def n_nodes(self):
        return self.adj.shape[0]

    @property
    def n_edges(self):
        return self.adj.nnz // 2


def load_orbit(path):
    df = pd.read_csv(path, sep=" ", header=None, dtype=np.uint64, engine="c")
    both = np.concatenate([df[0].to_numpy(), df[1].to_numpy()])
    del df
    codes, labels = pd.factorize(both)
    del both
    m = codes.size // 2
    u = codes[:m].astype(np.int32)
    v = codes[m:].astype(np.int32)
    del codes
    loops = int(np.count_nonzero(u == v))
    keep = u != v
    u, v = u[keep], v[keep]
    n = len(labels)
    adj = csr_matrix(
        (np.ones(2 * len(u), dtype=np.int8),
         (np.concatenate([u, v]), np.concatenate([v, u]))),
        shape=(n, n))
    return OG(adj, np.asarray(labels), loops)


def to_networkx(og):
    import networkx as nx
    g = nx.from_scipy_sparse_array(og.adj)
    return nx.relabel_nodes(g, dict(enumerate(og.labels.tolist())))


def _dists_from(og, sources):
    return shortest_path(og.adj, method="D", unweighted=True,
                         indices=np.asarray(sources, dtype=np.int64))


def _bfs_reduce(og, sources, workers=1, mem_doubles=150_000_000):
    src = np.asarray(sources, dtype=np.int64)
    k, n = len(src), og.n_nodes
    workers = max(1, min(workers, k))
    chunk = max(1, min(k, int(mem_doubles / max(n, 1) / workers)))
    sums = np.empty(k)
    maxs = np.empty(k)

    def do(lo):
        d = _dists_from(og, src[lo:lo + chunk])
        if np.isinf(d).any():
            raise ValueError("orbit graph is not connected")
        sums[lo:lo + d.shape[0]] = d.sum(axis=1)
        maxs[lo:lo + d.shape[0]] = d.max(axis=1)

    starts = range(0, k, chunk)
    if workers > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(do, starts))
    else:
        for lo in starts:
            do(lo)
    return sums, maxs


# ------------------------------------------------------------------- ASPL

def aspl_exact(og, chunk=512, workers=1):
    del chunk
    n = og.n_nodes
    sums, _ = _bfs_reduce(og, np.arange(n), workers=workers)
    return sums.sum() / (n * (n - 1))


def aspl_sampled(og, n_sources=400, seed=0, workers=1):
   
    n = og.n_nodes
    rng = np.random.default_rng(seed)
    k = min(n_sources, n)
    src = rng.choice(n, size=k, replace=False)
    sums, _ = _bfs_reduce(og, src, workers=workers)
    per_source = sums / (n - 1)
    return per_source.mean(), per_source.std(ddof=1) / np.sqrt(k)


# ------------------------------------------- exact stats in one BFS pass

def exact_bfs_stats(og, chunk=1024, workers=1, mem_doubles=150_000_000):
    del chunk
    n = og.n_nodes
    sums, ecc = _bfs_reduce(og, np.arange(n), workers=workers,
                            mem_doubles=mem_doubles)
    return sums.sum() / (n * (n - 1)), int(ecc.max()), ecc


def diameter_sampled(og, n_sources=1000, seed=0, workers=1):
    n = og.n_nodes
    rng = np.random.default_rng(seed)
    src = rng.choice(n, size=min(n_sources, n), replace=False)
    _, ecc = _bfs_reduce(og, src, workers=workers)
    m = int(ecc.max())
    return m, float((ecc == m).mean())


def estimate_bfs_cost(og, probe=64):
    import time
    t = time.time()
    _dists_from(og, np.arange(min(probe, og.n_nodes)))
    return (time.time() - t) / min(probe, og.n_nodes)


# ------------------------------------------------------------- colouring

def greedy_chromatic(og):
    indptr, indices = og.adj.indptr, og.adj.indices
    n = og.n_nodes
    color = np.full(n, -1, dtype=np.int64)
    for v in np.argsort(-np.diff(indptr), kind="stable"):
        used = color[indices[indptr[v]:indptr[v + 1]]]
        used = np.unique(used[used >= 0])
        free = np.flatnonzero(used != np.arange(len(used)))
        color[v] = free[0] if free.size else len(used)
    return int(color.max()) + 1



def orbit_stats(og, budget_seconds=60.0, workers=1, seed=0,
                sample_sources=1000, mem_doubles=150_000_000):
    n = og.n_nodes
    degrees = np.diff(og.adj.indptr)
    base = dict(nodes=n, edges=og.n_edges, selfloops=og.n_selfloops,
                density=2 * og.n_edges / (n * (n - 1)) if n > 1 else 0.0,
                min_degree=int(degrees.min()) if n else 0,
                max_degree=int(degrees.max()) if n else 0,
                chi_og=greedy_chromatic(og))
    cost = estimate_bfs_cost(og) * n / max(workers, 1)
    if cost <= budget_seconds:
        aspl, dia, _ = exact_bfs_stats(og, workers=workers,
                                       mem_doubles=mem_doubles)
        base.update(aspl=aspl, aspl_se=0.0, diameter=dia,
                    diameter_exactly_known=True)
    else:
        rng = np.random.default_rng(seed)
        k = min(sample_sources, n)
        src = rng.choice(n, size=k, replace=False)
        sums, ecc = _bfs_reduce(og, src, workers=workers,
                                mem_doubles=mem_doubles)
        per_source = sums / (n - 1)
        base.update(aspl=per_source.mean(),
                    aspl_se=per_source.std(ddof=1) / np.sqrt(k),
                    diameter=int(ecc.max()),
                    diameter_exactly_known=False)
    return base


def approximate_aspl_pairs(og, sample_size=10000, seed=0):
    n = og.n_nodes
    rng = np.random.default_rng(seed)
    total = 0.0
    for _ in range(sample_size):
        u, v = rng.choice(n, size=2, replace=False)
        total += _dists_from(og, [u])[0][v]
    return total / sample_size
