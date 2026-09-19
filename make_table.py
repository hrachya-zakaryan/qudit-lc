"""Build the per-orbit data table (Tables in App. B)
for any prime d, from the orbit directories written by orbits_mt.c.

    python make_table.py 3 7 --budget 20

Column conventions, matching the paper:
    orbit          global index; n ascending, then representatives ordered
                   by (#edges, total weight, code)
    |V|            number of graphs in the orbit
    |e|            edges of the representative
    chi_OG         greedy (largest-first) colouring of the orbit graph --
                   an upper bound on its chromatic number
    ln(N_L+1)      N_L = orbit-graph self-loops
    chi_i          minimal chromatic number over the orbit's members
    D              orbit-graph density
    <d_OG>         average shortest path length (exact when the all-sources
                   pass fits --budget seconds, else source-sampled, marked *)
    d_OG^max       diameter (sampled values are lower bounds, marked *)
    deg(g)_min     minimal over the orbit's members of the member's own
                   maximum vertex degree
    deg(OG)_max    maximal orbit-graph degree
    E_S            Schmidt measure
The expensive orbit-graph statistics are cached in <out>_stats_cache.csv,
so an interrupted run resumes where it stopped.
"""

import argparse
import csv
import os
import time

import numpy as np
import pandas as pd

from scan_tool import run_scan
from fast_analysis import load_orbit, orbit_stats

CACHE_COLS = ["n", "orbit", "nodes", "edges", "selfloops", "density",
              "min_degree", "max_degree", "aspl", "aspl_se", "diameter",
              "diameter_exactly_known", "chi_og"]


def _append_cache(path, row):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, CACHE_COLS)
        if new:
            w.writeheader()
        w.writerow({k: row[k] for k in CACHE_COLS})


def _stats_one(task):
    (d, root, dir_fmt, n, orbit, members, loops,
     budget, sources, mem_doubles) = task
    t = time.time()
    path = os.path.join(root, dir_fmt.format(d=d, n=n), f"orbit_{orbit}")
    og = load_orbit(path)
    if og.n_nodes != members or og.n_selfloops != loops:
        raise AssertionError(
            f"{path}: scanner and loader disagree "
            f"({og.n_nodes}/{members} nodes, "
            f"{og.n_selfloops}/{loops} loops)")
    s = orbit_stats(og, budget_seconds=budget, workers=1,
                    sample_sources=sources, mem_doubles=mem_doubles)
    s.update(n=n, orbit=orbit, secs=time.time() - t)
    return s


def _available_ram():
    if os.name == "nt":
        import ctypes

        class _MemStat(ctypes.Structure):
            _fields_ = ([("dwLength", ctypes.c_uint32),
                         ("dwMemoryLoad", ctypes.c_uint32)] +
                        [(f, ctypes.c_uint64) for f in
                         ("ullTotalPhys", "ullAvailPhys", "ullTotalPageFile",
                          "ullAvailPageFile", "ullTotalVirtual",
                          "ullAvailVirtual", "ullAvailExtendedVirtual")])

        st = _MemStat(dwLength=ctypes.sizeof(_MemStat))
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
            return int(st.ullAvailPhys)
        return 8 << 30
    try:
        return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (ValueError, OSError):
        return 8 << 30


def gather(d, nmax, dir_fmt, root, budget, workers, sources, cache_path,
           verbose=True, mem_doubles=150_000_000, ram_bytes=None):
    
    rows = run_scan(d, nmax, dir_fmt=dir_fmt, root=root)
    cache = {}
    if os.path.exists(cache_path):
        for r in pd.read_csv(cache_path,
                             float_precision="round_trip").to_dict("records"):
            cache[(int(r["n"]), int(r["orbit"]))] = r
    todo = [r for r in rows if (r["n"], r["orbit"]) not in cache]
    if todo:
        todo.sort(key=lambda r: r["members"], reverse=True)
        procs = max(1, min(workers, len(todo)))
        per_proc_doubles = max(mem_doubles // procs, 1_000_000)
        ram = ram_bytes or int(0.75 * _available_ram())
        pending = []
        for r in todo:
            path = os.path.join(root, dir_fmt.format(d=d, n=r["n"]),
                                f"orbit_{r['orbit']}")
            est = int(1.5 * os.path.getsize(path)) + (400 << 20)
            task = (d, root, dir_fmt, r["n"], r["orbit"], r["members"],
                    r["loops"], budget, sources, per_proc_doubles)
            pending.append([task, min(est, ram)])
        if verbose:
            print(f"{len(rows) - len(todo)} orbits cached, {len(todo)} to "
                  f"compute; {procs} processes, {ram / 2**30:.1f} GiB RAM "
                  "budget", flush=True)
        from concurrent.futures import (ProcessPoolExecutor, wait,
                                        FIRST_COMPLETED)
        futs, inflight, retried = {}, 0, set()
        with ProcessPoolExecutor(procs) as ex:
            try:
                def admit():
                    nonlocal inflight
                    i = 0
                    while i < len(pending) and len(futs) < procs:
                        task, est = pending[i]
                        if futs and inflight + est > ram:
                            i += 1          # too big right now; try a smaller one
                            continue
                        futs[ex.submit(_stats_one, task)] = (task, est)
                        inflight += est
                        pending.pop(i)

                admit()
                while futs:
                    done, _ = wait(list(futs), return_when=FIRST_COMPLETED)
                    for fut in done:
                        task, est = futs.pop(fut)
                        inflight -= est
                        try:
                            s = fut.result()
                        except MemoryError:
                            key = (task[3], task[4])
                            if key in retried:
                                raise
                            retried.add(key)
                            pending.insert(0, [task, ram])  # rerun solo
                            if verbose:
                                print(f"  n={key[0]} orbit_{key[1]}: out of "
                                      "memory, will retry alone", flush=True)
                            continue
                        _append_cache(cache_path, s)
                        cache[(s["n"], s["orbit"])] = s
                        if verbose:
                            note = "" if s["aspl_se"] == 0 else " (sampled)"
                            print(f"  n={s['n']} orbit_{s['orbit']}: "
                                  f"{s['nodes']} nodes  "
                                  f"ASPL {s['aspl']:.4f}{note}  "
                                  f"diam {s['diameter']}  [{s['secs']:.1f}s]",
                                  flush=True)
                    admit()
            except BaseException:
                for f in futs:
                    f.cancel()
                for p in getattr(ex, "_processes", {}).values():
                    try:
                        p.terminate()
                    except OSError:
                        pass
                raise
    out = []
    for r in rows:
        s = cache[(r["n"], r["orbit"])]
        if int(s["nodes"]) != r["members"]:
            raise AssertionError(
                f"stats cache is stale for n={r['n']} orbit_{r['orbit']}: "
                f"cached {s['nodes']} nodes vs {r['members']} members in "
                f"the scan -- delete the cache and rerun")
        out.append({**r, **{k: s[k] for k in CACHE_COLS[2:]}})
    return pd.DataFrame(out)


def number_orbits(df):
    df = df.sort_values(["n", "rep_edges", "rep_weight", "rep"],
                        kind="stable").reset_index(drop=True)
    df.insert(0, "no", np.arange(1, len(df) + 1))
    return df


def _es_cell(lo, up):
    return f"${lo}$" if lo == up else f"$({lo}, {up})$"


def latex_table(df, d, path):
    nmin, nmax = int(df["n"].min()), int(df["n"].max())
    n_sampled = int((df["aspl_se"] > 0).sum())
    se_max = df["aspl_se"].max()
    star_note = ""
    if n_sampled:
        star_note = (
            f" Starred entries are source-sampled ({n_sampled} orbits): "
            f"$\\langle d_{{OG}}\\rangle$ carries a standard error below "
            f"${se_max:.3f}$ and $d_{{OG}}^{{\\text{{max}}}}$ is a sampled "
            "lower bound.")
    caption = (
        f"Per-orbit data for $d={d}$ graph states on $n={nmin},\\dots,{nmax}$ "
        "qudits, columns as in arXiv:2506.05478: orbit size $|V|$, edges "
        "$|e|$ of the representative (minimal under (edges, total weight, "
        "code)), greedy bound $\\chi_{OG}$ on the orbit-graph chromatic "
        "number, self-loops $N_{\\mathcal{L}}$, minimal member chromatic "
        "number $\\chi_i$, orbit-graph density $\\mathcal{D}$, mean distance "
        "$\\langle d_{OG}\\rangle$, diameter $d_{OG}^{\\text{max}}$, "
        "minimal member max-degree $\\deg(g)_{\\text{min}}$, maximal "
        "orbit-graph degree $\\deg(OG)_{\\text{max}}$, and the Schmidt "
        "measure $E_S$ in "
        "units of $\\log_" + str(d) + "$ (exact value where the cut-rank "
        "lower bound meets the measurement upper bound, otherwise the "
        "interval)." + star_note)
    head = (r"orbit & $n$ & $|V|$ & $|e|$ & $\chi_{OG}$ & "
            r"$\ln(N_{\mathcal{L}}+1)$ & $\chi_i$ & $\mathcal{D}$ & "
            r"$\langle d_{OG}\rangle$ & $d_{OG}^{\text{max}}$ & "
            r"$\deg(g)_{\text{min}}$ & $\deg(OG)_{\text{max}}$ & $E_S$\\")
    lines = [
        "% generated by make_table.py -- requires \\usepackage{booktabs,longtable}",
        r"\begingroup",
        r"\setlength{\tabcolsep}{4.5pt}",
        r"\begin{longtable}{rrrrrrrrrrrrc}",
        r"\caption{" + caption + r"}\label{tab:orbit_data_d" + str(d) + r"}\\",
        r"\toprule", head, r"\midrule", r"\endfirsthead",
        r"\toprule", head, r"\midrule", r"\endhead",
        r"\midrule",
        r"\multicolumn{13}{r}{\emph{continued on the next page}}\\",
        r"\endfoot",
        r"\bottomrule", r"\endlastfoot",
    ]
    for r in df.itertuples():
        star = "" if r.aspl_se == 0 else r"$^*$"
        lines.append(
            f"{r.no} & {r.n} & {r.nodes} & {r.rep_edges} & {r.chi_og} & "
            f"{np.log(r.loop_nodes + 1):.2f} & {r.chi_min} & "
            f"{r.density:.5f} & {r.aspl:.2f}{star} & {r.diameter}{star} & "
            f"{r.mindeg_g} & {r.max_degree} & "
            f"{_es_cell(r.es_lo, r.es_up)}\\\\")
    lines += [r"\end{longtable}", r"\endgroup", ""]
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return path


def main():
    ap = argparse.ArgumentParser(
        description="Per-orbit LaTeX table")
    ap.add_argument("d", type=int, nargs="?", default=3)
    ap.add_argument("nmax", type=int, nargs="?", default=7)
    ap.add_argument("--root", default=".")
    ap.add_argument("--dir-fmt", default="orbits_d{d}_n{n}_separated")
    ap.add_argument("--budget", type=float, default=20.0,
                    help="single-core seconds allowed for an exact "
                         "all-sources BFS pass")
    ap.add_argument("--sources", type=int, default=1000,
                    help="BFS sources for sampled ASPL/diameter")
    ap.add_argument("--workers", type=int, default=os.cpu_count(),
                    help="orbits processed in parallel (processes)")
    ap.add_argument("--out", default=None, help="output prefix")
    args = ap.parse_args()

    prefix = args.out or f"table_d{args.d}_n{args.nmax}"
    df = gather(args.d, args.nmax, args.dir_fmt, args.root, args.budget,
                args.workers, args.sources, prefix + "_stats_cache.csv")
    df = number_orbits(df)
    df.to_csv(prefix + ".csv", index=False)
    tex = latex_table(df, args.d, prefix + ".tex")

    exact = int((df["es_lo"] == df["es_up"]).sum())
    print(f"\n{len(df)} orbits, {df['members'].sum()} graphs; "
          f"E_S exact for {exact}, interval for {len(df) - exact}")
    print(f"wrote {prefix}.csv and {tex}")


if __name__ == "__main__":
    main()
