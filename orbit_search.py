import taichi as ti
import numpy as np
from ls_lc import local_complementation, local_scaling

@ti.kernel
def orbit_bfs(G:ti.types.ndarray(dtype=ti.i32, ndim=2,)):
    return