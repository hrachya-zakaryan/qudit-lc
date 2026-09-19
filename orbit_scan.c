/* orbit_scan.c -- single-pass per-orbit analysis of LC-orbit edge lists.
 *
 * Reads the orbit_* files written by orbits_mt.c (any prime d) and computes,
 * per orbit, everything that needs the member list but not the orbit-graph
 * topology:
 *
 *   members     number of distinct graphs in the orbit (|V| of the OG)
 *   loops       raw self-loop lines "g g" in the file (fixed points, N_L)
 *   rep         representative: minimum under (#edges, total weight, code)
 *   es_lo       Schmidt-measure lower bound: max over bipartitions (A,B) of
 *               rank_{F_d}(Gamma_AB), computed on the representative only --
 *               the Schmidt rank across any cut is d^rank and local Cliffords
 *               cannot change it, so this is an orbit invariant.
 *   vc_min      min over members of the minimal vertex cover.
 *   chi_min     min over members of the chromatic number of the support graph.
 *   es_up       Schmidt-measure upper bound.
 *
 * build:  gcc -O3 -march=native -fopenmp -o orbit_scan orbit_scan.c
 * usage:  orbit_scan <d> <nmax> <dirfmt-with-%d> [--nmin K] [--no-dp]
 *                    [--threads T]
 * output: TSV on stdout, one row per orbit:
 *         n orbit members loops rep rep_edges rep_weight es_lo vc_min chi_min es_up
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#define MAXN 9

static int D, B;
static int PIDX[MAXN + 1][MAXN][MAXN];    /* slot of pair (i,j) on n verts   */
static uint64_t MASK[MAXN + 1][MAXN + 1]; /* decided slots, see orbits_mt.c  */
static int PI[MAXN + 1][MAXN * (MAXN - 1) / 2];  /* slot -> smaller vertex   */
static int PJ[MAXN + 1][MAXN * (MAXN - 1) / 2];  /* slot -> larger vertex    */
static uint64_t INC[MAXN + 1][MAXN];      /* INC[n][v]: slots touching v     */

typedef uint8_t Mat[MAXN][MAXN];

static void init_tables(void) {
    B = 1; while ((1 << B) < D) B++;
    for (int n = 1; n <= MAXN; n++) {
        int p = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                PIDX[n][i][j] = PIDX[n][j][i] = p;
                PI[n][p] = i; PJ[n][p] = j;
                p++;
            }
        uint64_t m = ((uint64_t)1 << B) - 1;
        for (int pos = 0; pos <= n; pos++) {
            uint64_t acc = 0;
            for (int i = pos; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    acc |= m << (PIDX[n][i][j] * B);
            MASK[n][pos] = acc;
        }
        for (int v = 0; v < n; v++) {
            uint64_t acc = 0;
            for (int q = 0; q < p; q++)
                if (PI[n][q] == v || PJ[n][q] == v) acc |= (uint64_t)1 << q;
            INC[n][v] = acc;
        }
    }
}

static inline uint64_t encode(const Mat m, int n) {
    uint64_t c = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            c |= (uint64_t)m[i][j] << (PIDX[n][i][j] * B);
    return c;
}

static inline void decode(uint64_t c, int n, Mat m) {
    uint64_t mask = ((uint64_t)1 << B) - 1;
    memset(m, 0, sizeof(Mat));
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            m[i][j] = m[j][i] = (c >> (PIDX[n][i][j] * B)) & mask;
}

/* ---- canonical form -------- */

typedef struct {
    Mat      m;
    int      n;
    uint64_t best;
    int      perm[MAXN];
    bool     used[MAXN];
} Canon;

static void canon_dfs(Canon *restrict c, int n, int pos, uint64_t cur) {
    if (pos < 0) { if (cur < c->best) c->best = cur; return; }

    const uint64_t mk   = MASK[n][pos];
    const int *restrict pidx = PIDX[n][pos];
    uint64_t best = c->best;

    int      cand[MAXN], ncand = 0;
    uint64_t code[MAXN];

    for (int v = 0; v < n; v++) {
        if (c->used[v]) continue;
        uint64_t nc = cur;
        const uint8_t *restrict row = c->m[v];
        for (int j = pos + 1; j < n; j++)
            nc |= (uint64_t)row[c->perm[j]] << (pidx[j] * B);
        if ((nc & mk) > (best & mk)) continue;         /* prune */
        int k = ncand++;
        while (k > 0 && (code[k - 1] & mk) > (nc & mk)) { code[k] = code[k-1]; cand[k] = cand[k-1]; k--; }
        code[k] = nc; cand[k] = v;
    }
    for (int t = 0; t < ncand; t++) {
        best = c->best;
        if ((code[t] & mk) > (best & mk)) break;
        c->used[cand[t]] = true; c->perm[pos] = cand[t];
        canon_dfs(c, n, pos - 1, code[t]);
        c->used[cand[t]] = false;
    }
}

static uint64_t canon(const Mat m, int n) {
    Canon c;
    memcpy(c.m, m, sizeof(Mat));
    c.n = n;
    c.best = encode(m, n);
    memset(c.used, 0, sizeof(c.used));
    canon_dfs(&c, n, n - 1, 0);
    return c.best;
}

/* ---- growable open-addressing set (single-threaded, per file) ---------- */

static inline uint64_t mix(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

typedef struct { uint64_t *sl; uint64_t mask, count; } Set;

static void set_init(Set *s, int logsize) {
    s->mask = ((uint64_t)1 << logsize) - 1;
    s->sl = calloc(s->mask + 1, sizeof(uint64_t));
    if (!s->sl) { fprintf(stderr, "set alloc failed\n"); exit(1); }
    s->count = 0;
}

static bool set_add(Set *s, uint64_t key);

static void set_grow(Set *s) {
    Set big;
    uint64_t oldmask = s->mask, *old = s->sl;
    int log2 = 0; while (((oldmask + 1) >> log2) > 1) log2++;
    set_init(&big, log2 + 1);
    for (uint64_t i = 0; i <= oldmask; i++)
        if (old[i]) set_add(&big, old[i] - 1);
    free(old);
    *s = big;
}

static bool set_add(Set *s, uint64_t key) {
    if (10 * (s->count + 1) > 7 * (s->mask + 1)) set_grow(s);
    uint64_t k = key + 1, i = mix(key) & s->mask;
    for (;;) {
        if (s->sl[i] == k) return false;
        if (s->sl[i] == 0) { s->sl[i] = k; s->count++; return true; }
        i = (i + 1) & s->mask;
    }
}


typedef struct { uint64_t *sl; uint32_t *val; uint64_t mask; } Map;

static void map_init(Map *m, uint64_t expect) {
    int log2 = 4; while (((uint64_t)1 << log2) < 2 * expect + 16) log2++;
    m->mask = ((uint64_t)1 << log2) - 1;
    m->sl = calloc(m->mask + 1, sizeof(uint64_t));
    m->val = malloc((m->mask + 1) * sizeof(uint32_t));
    if (!m->sl || !m->val) { fprintf(stderr, "map alloc failed\n"); exit(1); }
}

static void map_put(Map *m, uint64_t key, uint32_t v) {
    uint64_t k = key + 1, i = mix(key) & m->mask;
    while (m->sl[i] && m->sl[i] != k) i = (i + 1) & m->mask;
    m->sl[i] = k; m->val[i] = v;
}

static int64_t map_get(const Map *m, uint64_t key) {
    uint64_t k = key + 1, i = mix(key) & m->mask;
    for (;;) {
        if (m->sl[i] == k) return m->val[i];
        if (m->sl[i] == 0) return -1;
        i = (i + 1) & m->mask;
    }
}

static inline void edges_weight(uint64_t code, int P, int *e, int *w) {
    uint64_t fm = ((uint64_t)1 << B) - 1;
    int ec = 0, wc = 0;
    for (int p = 0; p < P; p++) {
        int f = (code >> (p * B)) & fm;
        ec += f != 0; wc += f;
    }
    *e = ec; *w = wc;
}

static inline uint64_t occupancy(uint64_t code, int P) {
    uint64_t fm = ((uint64_t)1 << B) - 1, m = 0;
    for (int p = 0; p < P; p++)
        if ((code >> (p * B)) & fm) m |= (uint64_t)1 << p;
    return m;
}

/* minimal vertex cover of a pair-occupancy mask: branch on an uncovered edge */
static int vc_rec(uint64_t mask, int n) {
    if (!mask) return 0;
    int p = __builtin_ctzll(mask);
    int a = 1 + vc_rec(mask & ~INC[n][PI[n][p]], n);
    if (a == 1) return 1;
    int b = 1 + vc_rec(mask & ~INC[n][PJ[n][p]], n);
    return a < b ? a : b;
}

/* chromatic number of a support graph */
static bool color_try(int v, int n, int k, uint8_t *col, int used,
                      const bool adj[MAXN][MAXN]) {
    if (v == n) return true;
    int lim = used + 1 < k ? used + 1 : k;
    for (int c = 0; c < lim; c++) {
        bool ok = true;
        for (int u = 0; u < v; u++)
            if (adj[u][v] && col[u] == c) { ok = false; break; }
        if (!ok) continue;
        col[v] = c;
        if (color_try(v + 1, n, k, col, c == used ? used + 1 : used, adj))
            return true;
    }
    return false;
}

static int chromatic(uint64_t mask, int n) {
    bool adj[MAXN][MAXN] = {{false}};
    for (uint64_t m = mask; m; m &= m - 1) {
        int p = __builtin_ctzll(m);
        adj[PI[n][p]][PJ[n][p]] = adj[PJ[n][p]][PI[n][p]] = true;
    }
    uint8_t col[MAXN];
    for (int k = 2; k <= n; k++)
        if (color_try(0, n, k, col, 0, adj)) return k;
    return n;
}

/* rank over F_d by Gaussian elimination (entries already reduced mod d) */
static int rank_fd(int rows, int cols, uint8_t a[MAXN][MAXN]) {
    int r = 0;
    for (int c = 0; c < cols && r < rows; c++) {
        int p = -1;
        for (int i = r; i < rows; i++) if (a[i][c]) { p = i; break; }
        if (p < 0) continue;
        if (p != r)
            for (int j = 0; j < cols; j++) {
                uint8_t t = a[p][j]; a[p][j] = a[r][j]; a[r][j] = t;
            }
        int inv = 1;
        for (int x = 1; x < D; x++) if (x * a[r][c] % D == 1) { inv = x; break; }
        for (int j = 0; j < cols; j++) a[r][j] = a[r][j] * inv % D;
        for (int i = 0; i < rows; i++) {
            if (i == r || !a[i][c]) continue;
            int f = a[i][c];
            for (int j = 0; j < cols; j++)
                a[i][j] = (a[i][j] + (D - f * a[r][j] % D)) % D;
        }
        r++;
    }
    return r;
}

/* max over bipartitions of rank_{F_d}(Gamma_AB); vertex 0 kept in A so each
   cut is visited once.  LC-invariant, so the representative suffices. */
static int es_lower(uint64_t code, int n) {
    Mat m; decode(code, n, m);
    int best = 0, half = n / 2;
    for (uint32_t A = 1; A < (uint32_t)(1 << n); A += 2) {
        int ka = __builtin_popcount(A);
        if (ka == n) continue;
        int ra[MAXN], rb[MAXN], na = 0, nb = 0;
        for (int v = 0; v < n; v++)
            if ((A >> v) & 1) ra[na++] = v; else rb[nb++] = v;
        int lo = na < nb ? na : nb;
        if (lo <= best) continue;                     /* rank <= min(|A|,|B|) */
        uint8_t sub[MAXN][MAXN];
        for (int i = 0; i < na; i++)
            for (int j = 0; j < nb; j++)
                sub[i][j] = m[ra[i]][rb[j]];
        int r = rank_fd(na, nb, sub);
        if (r > best) { best = r; if (best == half) return best; }
    }
    return best;
}


typedef struct {
    uint64_t *members; size_t nmem;
    uint64_t loops;                  /* raw "g g" lines (fixed-point moves)  */
    uint64_t loop_nodes;             /* distinct members with a self-loop    */
    uint64_t rep; int rep_e, rep_w;
    int es_lo, vc, chi, es_up;
    int mindeg_g;                    /* min over members of max weighted
                                        (multigraph) vertex degree           */
} Res;

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

static void scan_file(const char *path, int n, Res *r) {
    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno)); exit(1); }
    enum { BUF = 1 << 22 };
    char *buf = malloc(BUF);
    Set set, loopset;
    set_init(&set, 12); set_init(&loopset, 6);
    size_t cap = 1024;
    uint64_t *mem = malloc(cap * sizeof(uint64_t));
    size_t nmem = 0;
    uint64_t loops = 0, val = 0, pair[2];
    int innum = 0, pc = 0;
    for (;;) {
        size_t got = fread(buf, 1, BUF, f);
        if (!got) break;
        for (size_t i = 0; i < got; i++) {
            unsigned c = (unsigned char)buf[i] - '0';
            if (c <= 9) { val = val * 10 + c; innum = 1; continue; }
            if (!innum) continue;
            pair[pc++] = val; val = 0; innum = 0;
            if (pc == 2) {
                pc = 0;
                if (pair[0] == pair[1]) { loops++; set_add(&loopset, pair[0]); }
                for (int e = 0; e < 2; e++)
                    if (set_add(&set, pair[e])) {
                        if (nmem == cap) { cap *= 2; mem = realloc(mem, cap * sizeof(uint64_t)); }
                        mem[nmem++] = pair[e];
                    }
            }
        }
    }
    if (innum) pair[pc++] = val;
    if (pc == 2) {
        if (pair[0] == pair[1]) { loops++; set_add(&loopset, pair[0]); }
        for (int e = 0; e < 2; e++)
            if (set_add(&set, pair[e])) {
                if (nmem == cap) { cap *= 2; mem = realloc(mem, cap * sizeof(uint64_t)); }
                mem[nmem++] = pair[e];
            }
    }
    uint64_t loop_nodes = loopset.count;
    fclose(f); free(buf); free(set.sl); free(loopset.sl);
    if (!nmem) { fprintf(stderr, "empty orbit file: %s\n", path); exit(1); }

    int P = n * (n - 1) / 2;
    uint64_t rep = mem[0]; int re, rw;
    edges_weight(rep, P, &re, &rw);
    uint64_t *masks = malloc(nmem * sizeof(uint64_t));
    uint64_t fm = ((uint64_t)1 << B) - 1;
    int mindeg = n * D;
    for (size_t i = 0; i < nmem; i++) {
        int e, w;
        edges_weight(mem[i], P, &e, &w);
        if (e < re || (e == re && (w < rw || (w == rw && mem[i] < rep)))) {
            rep = mem[i]; re = e; rw = w;
        }
        masks[i] = occupancy(mem[i], P);
        int wdeg[MAXN] = {0}, md = 0;
        for (int p = 0; p < P; p++) {
            int f = (int)((mem[i] >> (p * B)) & fm);
            wdeg[PI[n][p]] += f; wdeg[PJ[n][p]] += f;
        }
        for (int vx = 0; vx < n; vx++) if (wdeg[vx] > md) md = wdeg[vx];
        if (md < mindeg) mindeg = md;
    }
    qsort(masks, nmem, sizeof(uint64_t), cmp_u64);
    int vc = n, chi = n;
    for (size_t i = 0; i < nmem; i++) {
        if (i && masks[i] == masks[i - 1]) continue;
        int v = vc_rec(masks[i], n);
        if (v < vc) vc = v;
        int c = chromatic(masks[i], n);
        if (c < chi) chi = c;
    }
    free(masks);

    r->members = mem; r->nmem = nmem;
    r->loops = loops; r->loop_nodes = loop_nodes;
    r->rep = rep; r->rep_e = re; r->rep_w = rw;
    r->es_lo = es_lower(rep, n);
    r->vc = vc; r->chi = chi; r->es_up = vc;
    r->mindeg_g = mindeg;
}


static Map LVLMAP[MAXN + 1];      
static int *LVLU[MAXN + 1];       


static long long deleted_value(const Mat m, int n, int v) {
    bool seen[MAXN] = {false};
    seen[v] = true;
    long long total = 0;
    for (int s0 = 0; s0 < n; s0++) {
        if (seen[s0]) continue;
        int comp[MAXN], nc = 0, stack[MAXN], sp = 0;
        seen[s0] = true; stack[sp++] = s0;
        while (sp) {
            int x = stack[--sp];
            comp[nc++] = x;
            for (int u = 0; u < n; u++)
                if (!seen[u] && u != v && m[x][u]) { seen[u] = true; stack[sp++] = u; }
        }
        if (nc == 1) continue;
        if (nc == 2) { total += 1; continue; }
        Mat sub; memset(sub, 0, sizeof(Mat));
        for (int a = 0; a < nc; a++)
            for (int b = a + 1; b < nc; b++)
                sub[a][b] = sub[b][a] = m[comp[a]][comp[b]];
        int64_t oid = map_get(&LVLMAP[nc], canon(sub, nc));
        if (oid < 0) return -1;
        total += LVLU[nc][oid];
    }
    return total;
}

static void dp_orbit(Res *r, int n) {
    if (r->es_up <= r->es_lo) return;          
    int best = r->es_up;
    for (size_t i = 0; i < r->nmem; i++) {
        Mat m; decode(r->members[i], n, m);
        for (int v = 0; v < n; v++) {
            long long t = deleted_value(m, n, v);
            if (t < 0) {
                fprintf(stderr, "DP lookup miss (n=%d, member %llu, v=%d)\n",
                        n, (unsigned long long)r->members[i], v);
                exit(1);
            }
            if (1 + t < best) best = (int)(1 + t);
        }
        if (best <= r->es_lo) break;
    }
    if (best < r->es_lo)
        fprintf(stderr, "WARNING: DP bound %d below lower bound %d (n=%d) -- "
                        "this indicates a bug\n", best, r->es_lo, n);
    r->es_up = best;
}



int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr,
            "usage: %s <d> <nmax> <dirfmt-with-%%d> [--nmin K] [--no-dp] [--threads T]\n"
            "e.g.:  %s 3 7 \"orbits_d3_n%%d_separated\"\n", argv[0], argv[0]);
        return 1;
    }
    D = atoi(argv[1]);
    int nmax = atoi(argv[2]);
    const char *fmt = argv[3];
    int nmin = 3, nodp = 0, nthreads = 0;
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "--nmin") && i + 1 < argc) nmin = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no-dp")) nodp = 1;
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else { fprintf(stderr, "unknown option: %s\n", argv[i]); return 1; }
    }
    if (D < 2 || D > 13) { fprintf(stderr, "d out of range\n"); return 1; }
    if (nmax > MAXN) { fprintf(stderr, "nmax > %d\n", MAXN); return 1; }
    init_tables();
    if ((long)(nmax * (nmax - 1) / 2) * B > 64) {
        fprintf(stderr, "code exceeds 64 bits\n"); return 1;
    }
    if (!nodp) nmin = 3;
#ifdef _OPENMP
    if (nthreads > 0) omp_set_num_threads(nthreads);
#else
    (void)nthreads;
#endif

    printf("# orbit_scan d=%d nmin=%d nmax=%d dp=%d\n", D, nmin, nmax, !nodp);
    printf("# n orbit members loops loop_nodes rep rep_edges rep_weight"
           " es_lo vc_min chi_min es_up mindeg_g\n");

    for (int n = nmin; n <= nmax; n++) {
        char dir[1024], path[1200];
        snprintf(dir, sizeof dir, fmt, n);
        int norb = 0;
        for (;;) {
            snprintf(path, sizeof path, "%s/orbit_%d", dir, norb);
            FILE *f = fopen(path, "rb");
            if (!f) break;
            fclose(f); norb++;
        }
        if (!norb) { fprintf(stderr, "no orbit files in %s\n", dir); return 1; }

        Res *res = calloc(norb, sizeof(Res));
        #pragma omp parallel for schedule(dynamic, 1)
        for (int o = 0; o < norb; o++) {
            char p[1200];
            snprintf(p, sizeof p, "%s/orbit_%d", dir, o);
            scan_file(p, n, &res[o]);
        }

        if (!nodp) {
            #pragma omp parallel for schedule(dynamic, 1)
            for (int o = 0; o < norb; o++) dp_orbit(&res[o], n);
        }

        size_t total = 0;
        for (int o = 0; o < norb; o++) total += res[o].nmem;
        if (n < nmax) {
            map_init(&LVLMAP[n], total);
            LVLU[n] = malloc(norb * sizeof(int));
            for (int o = 0; o < norb; o++) {
                LVLU[n][o] = res[o].es_up;
                for (size_t i = 0; i < res[o].nmem; i++)
                    map_put(&LVLMAP[n], res[o].members[i], (uint32_t)o);
            }
        }

        for (int o = 0; o < norb; o++) {
            printf("%d %d %zu %llu %llu %llu %d %d %d %d %d %d %d\n",
                   n, o, res[o].nmem, (unsigned long long)res[o].loops,
                   (unsigned long long)res[o].loop_nodes,
                   (unsigned long long)res[o].rep, res[o].rep_e, res[o].rep_w,
                   res[o].es_lo, res[o].vc, res[o].chi, res[o].es_up,
                   res[o].mindeg_g);
            free(res[o].members);
        }
        fflush(stdout);
        fprintf(stderr, "n=%d: %d orbits, %zu graphs\n", n, norb, total);
        free(res);
    }
    return 0;
}
