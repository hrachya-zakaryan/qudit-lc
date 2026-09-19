"""Extract and draw the representatives of the LC orbits."""

import os
from itertools import permutations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Tried in order when no dir_fmt is given.  The "_sorted" variant is what
# sort_orbits() produces; the plain one is what orbits.c writes directly.
DEFAULT_DIR_FMTS = ("orbits_d{d}_n{n}_separated_sorted",
                    "orbits_d{d}_n{n}_separated")

# weight -> colour.  Index 0 is unused (no edge).  Weights 1 and 2 are kept
# as black/red so d = 3 figures are identical to the published ones; higher
# weights use the Okabe-Ito colourblind-safe palette.
WEIGHT_COLORS = ["none",       # 0  no edge
                 "black",      # 1
                 "red",        # 2
                 "#0072B2",    # 3  blue
                 "#009E73",    # 4  bluish green
                 "#CC79A7",    # 5  reddish purple
                 "#E69F00",    # 6  orange
                 "#56B4E9",    # 7  sky blue
                 "#8B4513",    # 8  brown
                 "#7F7F7F"]    # 9  grey

WEIGHT_STYLES = ["-", "-", (0, (5, 1.5)), (0, (1, 1.2)), (0, (5, 1.2, 1, 1.2)),
                 (0, (3, 1, 1, 1, 1, 1)), (0, (8, 2)), (0, (2, 2)),
                 (0, (6, 1, 2, 1)), (0, (4, 1, 4, 3))]

MAX_D = len(WEIGHT_COLORS)


# --------------------------------------------------------------- encoding

def bits_per_weight(d):
    """Bits used per edge weight, matching bitpack_encode."""
    return max(1, (d - 1).bit_length())


def pair_slots(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def decode(code, n, d):
    """Bitpacked integer -> n x n adjacency matrix of weights."""
    b = bits_per_weight(d)
    mask = (1 << b) - 1
    a = np.zeros((n, n), dtype=int)
    for p, (i, j) in enumerate(pair_slots(n)):
        a[i, j] = a[j, i] = (code >> (p * b)) & mask
    return a


def weight_list(code, n, d):
    b = bits_per_weight(d)
    mask = (1 << b) - 1
    m = n * (n - 1) // 2
    return [(code >> (p * b)) & mask for p in range(m)]


def codes_edges_weights(codes, n, d):
    codes = np.asarray(codes, dtype=np.uint64)
    b = bits_per_weight(d)
    fm = np.uint64((1 << b) - 1)
    edges = np.zeros(len(codes), dtype=np.int64)
    weights = np.zeros(len(codes), dtype=np.int64)
    for p in range(n * (n - 1) // 2):
        f = (codes >> np.uint64(b * p)) & fm
        edges += f != 0
        weights += f.astype(np.int64)
    return edges, weights


def sort_key(code, n, d, legacy=False):
    if legacy:
        if d != 3:
            raise ValueError(
                f"legacy=True reproduces the original tie-breaking, which is "
                f"only defined for d = 3 (got d = {d}).  The original "
                f"total_weight() strides 2 bits at a time and bit_count() only "
                f"equals the edge count when d = 3.")
        return (code.bit_count(), _legacy_total_weight(code), str(code))
    w = weight_list(code, n, d)
    return (sum(1 for x in w if x), sum(w), code)


def _legacy_total_weight(code):
    return sum((code >> i) & 0b11 for i in range(0, code.bit_length(), 2))


# ------------------------------------------------------------ file access

def orbit_codes(path):
    import pandas as pd
    a = pd.read_csv(path, sep=" ", header=None, dtype=np.uint64,
                    engine="c").values
    return np.unique(a)


def orbit_representative(path, n, d, legacy=False):
    codes = orbit_codes(path)
    if legacy:
        return min((int(c) for c in codes),
                   key=lambda c: sort_key(c, n, d, legacy))
    edges, weights = codes_edges_weights(codes, n, d)
    return int(codes[np.lexsort((codes, weights, edges))[0]])


def _resolve_dir(root, d, n, dir_fmt):
    fmts = (dir_fmt,) if dir_fmt else DEFAULT_DIR_FMTS
    for fmt in fmts:
        path = os.path.join(root, fmt.format(d=d, n=n))
        if os.path.isdir(path):
            return path
    return None


def _c_engine_reps(d, ns, dir_fmt, root):
    ns = sorted(ns)
    if not ns or ns != list(range(ns[0], ns[-1] + 1)):
        return None
    for fmt in (dir_fmt,) if dir_fmt else DEFAULT_DIR_FMTS:
        if all(os.path.isdir(os.path.join(root, fmt.format(d=d, n=n)))
               for n in ns):
            try:
                from scan_tool import run_scan
                rows = run_scan(d, ns[-1], dir_fmt=fmt, root=root,
                                nmin=ns[0], dp=False)
            except Exception:
                return None
            return [[r["rep"], (r["n"], r["orbit"])] for r in rows
                    if r["n"] >= ns[0]]
    return None


def find_representatives(d, ns=range(3, 8), dir_fmt=None, root=".",
                         verbose=False, legacy=False, engine="auto"):
    if engine not in ("auto", "c", "numpy"):
        raise ValueError(f"unknown engine {engine!r}")
    if not legacy and engine in ("auto", "c"):
        reps = _c_engine_reps(d, ns, dir_fmt, root)
        if reps is not None:
            if verbose:
                for code, (n, j) in reps:
                    print(f"n={n} orbit {j}: {code}")
            return reps
        if engine == "c":
            raise RuntimeError(
                "engine='c' needs the orbit_scan helper, contiguous ns and "
                "a single directory convention -- see scan_tool.ensure_built")

    reps = []
    searched, found = [], []
    for n in ns:
        directory = _resolve_dir(root, d, n, dir_fmt)
        if directory is None:
            fmts = (dir_fmt,) if dir_fmt else DEFAULT_DIR_FMTS
            searched += [os.path.join(root, f.format(d=d, n=n)) for f in fmts]
            continue
        found.append(directory)
        n_orbits = len(os.listdir(directory))
        for j in range(n_orbits):
            path = os.path.join(directory, f"orbit_{j}")
            reps.append([orbit_representative(path, n, d, legacy), (n, j)])
            if verbose:
                print(f"n={n} orbit {j}: {reps[-1][0]}")

    if not reps:
        here = os.path.abspath(root)
        present = sorted(e for e in os.listdir(here)
                         if os.path.isdir(os.path.join(here, e)))
        raise FileNotFoundError(
            "No orbit directories found.\n"
            f"  looked in: {here}\n"
            "  tried:     " + "\n             ".join(searched) + "\n"
            "  present:   " + (", ".join(present) if present else "(no subdirectories)") + "\n"
            "Pass root=... if the orbit folders are elsewhere, or dir_fmt=... "
            "if they are named differently.")
    return reps


def orbit_value(path, n, d, legacy=False):
    return sort_key(orbit_representative(path, n, d, legacy), n, d, legacy)


def sort_orbits(d, ns=range(3, 8), root=".",
                in_fmt="orbits_d{d}_n{n}_separated",
                out_fmt="orbits_d{d}_n{n}_separated_sorted",
                legacy=False):
    import shutil
    written = []
    for n in ns:
        src_dir = os.path.join(root, in_fmt.format(d=d, n=n))
        if not os.path.isdir(src_dir):
            continue
        dst_dir = os.path.join(root, out_fmt.format(d=d, n=n))
        os.makedirs(dst_dir, exist_ok=True)
        names = sorted(os.listdir(src_dir),
                       key=_orbit_sort_keys(src_dir, n, d, in_fmt, root, legacy))
        for new_index, name in enumerate(names):
            shutil.copy2(os.path.join(src_dir, name),
                         os.path.join(dst_dir, f"orbit_{new_index}"))
        written.append((dst_dir, len(names)))
    return written


def _orbit_sort_keys(src_dir, n, d, in_fmt, root, legacy):
    if not legacy:
        try:
            from scan_tool import run_scan
            rows = run_scan(d, n, dir_fmt=in_fmt, root=root, nmin=n, dp=False)
            keys = {f"orbit_{r['orbit']}":
                    (r["rep_edges"], r["rep_weight"], r["rep"])
                    for r in rows if r["n"] == n}
            if keys:
                return keys.__getitem__
        except Exception:
            pass
    return lambda o: orbit_value(os.path.join(src_dir, o), n, d, legacy)


# ------------------------------------------------------- circular layout

def count_crossings(edges, pos):
    chords = []
    for u, v in edges:
        a, b = pos[u], pos[v]
        chords.append((a, b) if a < b else (b, a))
    crossings = 0
    for i in range(len(chords)):
        p1, q1 = chords[i]
        for j in range(i + 1, len(chords)):
            p2, q2 = chords[j]
            if p1 < p2 < q1 < q2 or p2 < p1 < q2 < q1:
                crossings += 1
    return crossings


def best_circular_order(adj):
    n = adj.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i, j]]
    if n < 4:
        return list(range(n)), 0

    best_order, best = list(range(n)), None
    rest = list(range(1, n))
    pos = [0] * n
    for tail in permutations(rest):
        if tail[0] > tail[-1]:
            continue
        order = (0,) + tail
        for idx, v in enumerate(order):
            pos[v] = idx
        c = count_crossings(edges, pos)
        if best is None or c < best:
            best, best_order = c, list(order)
            if best == 0:
                break
    return best_order, best


def circular_positions(order):
    n = len(order)
    step = 2 * np.pi / n
    return {v: (np.cos(i * step), np.sin(i * step)) for i, v in enumerate(order)}


# ------------------------------------------------------------- rendering

def draw_representatives(reps, d, cols=6, filename="Representatives.png",
                         dpi=300, node_size=80, label_fmt="No. {i}",
                         linestyles=None, lw=2):
    if d > MAX_D:
        raise ValueError(
            f"d = {d} needs {d - 1} edge weights but only {MAX_D - 1} are "
            f"defined.  Add entries to WEIGHT_COLORS and WEIGHT_STYLES.")
    if linestyles is None:
        linestyles = d >= 4
    if not reps:
        raise ValueError("draw_representatives() got an empty list of "
                         "representatives -- nothing to draw.")
    num = len(reps)
    cols = max(1, min(cols, num))
    rows = (num + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 1.6, rows * 1.6),
                             squeeze=False)
    axes = axes.flatten()

    for i, (code, (n, _orbit)) in enumerate(reps):
        ax = axes[i]
        adj = decode(code, n, d)
        order, _ = best_circular_order(adj)
        pos = circular_positions(order)

        for u in range(n):
            for v in range(u + 1, n):
                w = adj[u, v]
                if not w:
                    continue
                (x1, y1), (x2, y2) = pos[u], pos[v]
                ax.plot([x1, x2], [y1, y2], color=WEIGHT_COLORS[w], lw=lw,
                        linestyle=WEIGHT_STYLES[w] if linestyles else "-",
                        zorder=1, solid_capstyle="round")

        xs = [pos[v][0] for v in range(n)]
        ys = [pos[v][1] for v in range(n)]
        ax.scatter(xs, ys, s=node_size, c="blue", edgecolors="black",
                   linewidths=1, zorder=2)

        ax.set_title(label_fmt.format(i=i + 1, n=n), fontsize=16)
        ax.set_aspect("equal")
        ax.set_xlim(-1.25, 1.25)
        ax.set_ylim(-1.25, 1.25)
        ax.axis("off")

    for j in range(num, len(axes)):
        fig.delaxes(axes[j])

    legend = [mlines.Line2D([], [], color=WEIGHT_COLORS[w], lw=1.6,
                            linestyle=WEIGHT_STYLES[w] if linestyles else "-",
                            label=f"Edge weight = {w}")
              for w in range(1, d)]
    fig.legend(handles=legend, loc="lower right", fontsize=16, frameon=True,
               bbox_to_anchor=(0.99, 0.005))

    plt.tight_layout()
    plt.savefig(filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return filename


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Draw the orbit representatives.")
    ap.add_argument("d", type=int, nargs="?", default=3, help="local dimension")
    ap.add_argument("n_max", type=int, nargs="?", default=7, help="max vertices")
    ap.add_argument("--n-min", type=int, default=3)
    ap.add_argument("--root", default="data",
                    help="folder containing the orbit directories")
    ap.add_argument("--dir-fmt", default=None,
                    help="directory name pattern, e.g. orbits_d{d}_n{n}_separated")
    ap.add_argument("--cols", type=int, default=9)
    ap.add_argument("--out", default="Representatives.png")
    ap.add_argument("--legacy", action="store_true",
                    help="reproduce the original tie-breaking (d = 3 only)")
    ap.add_argument("--file-order", action="store_true",
                    help="number the panels by on-disk file index instead of "
                         "the representative order the tables use")
    args = ap.parse_args()

    reps = find_representatives(args.d, range(args.n_min, args.n_max + 1),
                                dir_fmt=args.dir_fmt, root=args.root,
                                legacy=args.legacy)
    if not args.file_order:
        # panel i then matches row `no` = i of the make_table/sort_orbits csv
        reps.sort(key=lambda r: (r[1][0],)
                  + sort_key(r[0], r[1][0], args.d, args.legacy))
    print(f"{len(reps)} representatives")
    print(draw_representatives(reps, args.d, cols=args.cols, filename=args.out))
