"""Build and run the orbit_scan C helper, and parse its output."""

import os
import shutil
import subprocess

COLS = ("n", "orbit", "members", "loops", "loop_nodes", "rep",
        "rep_edges", "rep_weight", "es_lo", "vc_min", "chi_min", "es_up",
        "mindeg_g")

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "orbit_scan.c")
_EXE = os.path.join(_HERE, "orbit_scan.exe" if os.name == "nt" else "orbit_scan")


def ensure_built(force=False):
    """Path of a usable orbit_scan binary, compiling it if needed.
    Returns None when neither a binary nor a C compiler is available."""
    have = os.path.exists(_EXE)
    fresh = have and (not os.path.exists(_SRC)
                      or os.path.getmtime(_EXE) >= os.path.getmtime(_SRC))
    if fresh and not force:
        return _EXE
    cc = shutil.which("gcc") or shutil.which("cc") or shutil.which("clang")
    if cc is None or not os.path.exists(_SRC):
        return _EXE if have else None
    for extra in (["-march=native", "-fopenmp"], ["-fopenmp"], []):
        r = subprocess.run([cc, "-O3", *extra, "-o", _EXE, _SRC],
                           capture_output=True, text=True)
        if r.returncode == 0:
            return _EXE
    return _EXE if have else None


def run_scan(d, nmax, dir_fmt="orbits_d{d}_n{n}_separated", root=".",
             nmin=3, dp=True, threads=0):
    exe = ensure_built()
    if exe is None:
        raise RuntimeError("orbit_scan is not built and no C compiler was "
                           "found -- use the numpy fallbacks instead")
    cfmt = dir_fmt.format(d=d, n="{n}").replace("{n}", "%d")
    args = [exe, str(d), str(nmax), cfmt]
    if not dp:
        args += ["--no-dp", "--nmin", str(nmin)]
    if threads:
        args += ["--threads", str(threads)]
    p = subprocess.run(args, cwd=root, capture_output=True, text=True)
    if p.returncode:
        raise RuntimeError(f"orbit_scan failed:\n{p.stderr}")
    rows = []
    for line in p.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        rows.append(dict(zip(COLS, map(int, line.split()))))
    return rows
