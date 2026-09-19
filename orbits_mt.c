/* orbits.c -- LC/LS orbit enumeration for qudit graph states.

 * usage:  ./orbits <n> <d> [--outdir-fmt FMT] [--upto] [--loghash K]
 *
 * --outdir-fmt writes orbit_0, orbit_1, ... edge lists in exactly the format
 * read by nx.read_edgelist().  It is a printf format taking one %d (the
 * vertex count), so one run with --upto fills every directory at once:
 *
 *   ./orbits 7 3 --upto --outdir-fmt "orbits_d3_n%d_separated" --loghash 23
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <errno.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#if defined(_WIN32)
  #include <direct.h>
  #define MKDIR(p)  _mkdir(p)
  #define IS_SEP(c) ((c) == '/' || (c) == '\\')
#else
  #include <sys/stat.h>
  #define MKDIR(p)  mkdir((p), 0777)
  #define IS_SEP(c) ((c) == '/')
#endif

#define MAXN 9

static int N, D, B;
static int PIDX[MAXN + 1][MAXN][MAXN];   /* PIDX[n][i][j] = slot index      */
static uint64_t MASK[MAXN + 1][MAXN + 1];/* MASK[n][pos]  = decided slots   */
static uint8_t  ADD[16][16], MUL[16][16];/* mod-d arithmetic tables         */

typedef uint8_t Mat[MAXN][MAXN];

static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

/* ---------------------------------------------------------------- tables */

static void init_tables(void) {
    B = 1; while ((1 << B) < D) B++;
    for (int n = 1; n <= MAXN; n++) {
        int p = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                PIDX[n][i][j] = PIDX[n][j][i] = p++;
        uint64_t m = ((uint64_t)1 << B) - 1;
        for (int pos = 0; pos <= n; pos++) {
            uint64_t acc = 0;
            for (int i = pos; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    acc |= m << (PIDX[n][i][j] * B);
            MASK[n][pos] = acc;
        }
    }
    for (int a = 0; a < D; a++)
        for (int b = 0; b < D; b++) { ADD[a][b] = (a + b) % D; MUL[a][b] = (a * b) % D; }
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

/* ------------------------------------------------- canonical form (B&B) */

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

/* ------------------------------------------- lock-free hash (no growth) */

static uint64_t *H;
static uint64_t  HMASK;

static inline uint64_t mix(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

static void hinit(int logsize) {
    HMASK = ((uint64_t)1 << logsize) - 1;
    H = calloc(HMASK + 1, sizeof(uint64_t));
    if (!H) { fprintf(stderr, "hash alloc failed (2^%d slots)\n", logsize); exit(1); }
}

/* insert-or-find; returns the slot, sets *is_new */
static inline uint64_t hput(uint64_t key, bool *is_new) {
    uint64_t k = key + 1, i = mix(key) & HMASK, probes = 0;
    for (;;) {
        uint64_t cur = __atomic_load_n(&H[i], __ATOMIC_ACQUIRE);
        if (cur == k) { *is_new = false; return i; }
        if (cur == 0) {
            uint64_t exp = 0;
            if (__atomic_compare_exchange_n(&H[i], &exp, k, false,
                    __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) { *is_new = true; return i; }
            if (exp == k) { *is_new = false; return i; }
        }
        i = (i + 1) & HMASK;
        if (++probes > HMASK) {
            fprintf(stderr, "hash full -- raise --loghash\n"); exit(1);
        }
    }
}

static inline uint64_t hget(uint64_t key) {
    uint64_t k = key + 1, i = mix(key) & HMASK;
    for (;;) {
        uint64_t cur = H[i];
        if (cur == k) return i;
        if (cur == 0) return UINT64_MAX;
        i = (i + 1) & HMASK;
    }
}

/* --------------------------------------------- lock-free union-find */

static uint32_t *UF;

static inline uint32_t uf_find(uint32_t x) {
    for (;;) {
        uint32_t p = __atomic_load_n(&UF[x], __ATOMIC_RELAXED);
        if (p == x) return x;
        uint32_t gp = __atomic_load_n(&UF[p], __ATOMIC_RELAXED);
        uint32_t exp = p;
        __atomic_compare_exchange_n(&UF[x], &exp, gp, false,
                                    __ATOMIC_RELAXED, __ATOMIC_RELAXED);
        x = gp;
    }
}

static inline void uf_union(uint32_t a, uint32_t b) {
    for (;;) {
        a = uf_find(a); b = uf_find(b);
        if (a == b) return;
        if (a > b) { uint32_t t = a; a = b; b = t; }
        uint32_t exp = b;
        if (__atomic_compare_exchange_n(&UF[b], &exp, a, false,
                                        __ATOMIC_ACQ_REL, __ATOMIC_RELAXED)) return;
    }
}

static inline void local_scaling(Mat m, int v, int k) {
    for (int u = 0; u < N; u++)
        if (u != v) m[v][u] = m[u][v] = MUL[k][m[v][u]];
}

static inline void local_complementation(Mat m, int v, int k) {
    int nb[MAXN], nn = 0;
    for (int u = 0; u < N; u++) if (u != v && m[v][u]) nb[nn++] = u;
    for (int a = 0; a < nn; a++) {
        int u = nb[a], t = MUL[k][m[v][u]];
        if (!t) continue;
        for (int b = a + 1; b < nn; b++) {
            int w = nb[b];
            m[u][w] = m[w][u] = ADD[m[u][w]][MUL[t][m[v][w]]];
        }
    }
}

static bool connected(const Mat m, int n) {
    int stack[MAXN], sp = 0, seen = 1;
    bool vis[MAXN] = {0};
    vis[0] = true; stack[sp++] = 0;
    while (sp) {
        int v = stack[--sp];
        for (int u = 0; u < n; u++)
            if (!vis[u] && m[v][u]) { vis[u] = true; seen++; stack[sp++] = u; }
    }
    return seen == n;
}


static inline int neighbours(uint64_t g, int n, uint64_t *out) {
    Mat m, tmp;
    decode(g, n, m);
    int c = 0;
    for (int v = 0; v < n; v++) {
        for (int k = 2; k < D; k++) {
            memcpy(tmp, m, sizeof(Mat));
            local_scaling(tmp, v, k);
            out[c++] = canon(tmp, n);
        }
        for (int k = 1; k < D; k++) {
            memcpy(tmp, m, sizeof(Mat));
            local_complementation(tmp, v, k);
            out[c++] = canon(tmp, n);
        }
    }
    return c;
}

static void mkdir_p(const char *p) {
    char buf[1024];
    snprintf(buf, sizeof buf, "%s", p);
    size_t start = 0;
#if defined(_WIN32)
    if (buf[0] && buf[1] == ':') start = 2;
    #define IS_SEP2(c) ((c) == '/' || (c) == '\\')
#else
    #define IS_SEP2(c) ((c) == '/')
#endif
    if (IS_SEP2(buf[start])) start++;
    for (char *q = buf + start; *q; q++) {
        if (!IS_SEP2(*q)) continue;
        char save = *q; *q = 0;
        if (*buf && MKDIR(buf) != 0 && errno != EEXIST) { perror(buf); exit(1); }
        *q = save;
    }
    if (*buf && MKDIR(buf) != 0 && errno != EEXIST) { perror(buf); exit(1); }
}


static void link_collect_output(int n, const char *outdir, double t_gen) {
    (void)t_gen;
    double t2 = now();
    UF = malloc((HMASK + 1) * sizeof(uint32_t));
    if (!UF) { fprintf(stderr, "union-find alloc failed\n"); exit(1); }
    #pragma omp parallel for
    for (long long i = 0; i <= (long long)HMASK; i++) UF[i] = (uint32_t)i;

    #pragma omp parallel for schedule(dynamic, 1024)
    for (long long i = 0; i <= (long long)HMASK; i++) {
        if (!H[i]) continue;
        uint64_t g = H[i] - 1;
        Mat m; decode(g, n, m);
        if (!connected(m, n)) continue;
        uint64_t nb[MAXN * 8];
        int c = neighbours(g, n, nb);
        for (int t = 0; t < c; t++) {
            uint64_t j = hget(nb[t]);
            if (j != UINT64_MAX) uf_union((uint32_t)i, (uint32_t)j);
        }
    }
    double t_link = now() - t2;

    double t3 = now();
    uint32_t *orbit_of = malloc((HMASK + 1) * sizeof(uint32_t));
    for (uint64_t i = 0; i <= HMASK; i++) orbit_of[i] = UINT32_MAX;
    size_t norbits = 0, ngraphs = 0;
    for (uint64_t i = 0; i <= HMASK; i++) {
        if (!H[i]) continue;
        Mat m; decode(H[i] - 1, n, m);
        if (!connected(m, n)) continue;
        uint32_t r = uf_find((uint32_t)i);
        if (orbit_of[r] == UINT32_MAX) orbit_of[r] = (uint32_t)norbits++;
        ngraphs++;
    }
    fprintf(stderr, "n=%d: %llu orbits, %llu graphs (link %.2fs, collect %.2fs)\n",
            n, (unsigned long long)norbits, (unsigned long long)ngraphs,
            t_link, now() - t3);
    printf("%d %llu %llu\n", n, (unsigned long long)norbits,
           (unsigned long long)ngraphs);

    if (outdir) {
        double t4 = now();
        mkdir_p(outdir);
        size_t *cnt = calloc(norbits + 1, sizeof(size_t));
        for (uint64_t i = 0; i <= HMASK; i++) {
            if (!H[i]) continue;
            Mat m; decode(H[i] - 1, n, m);
            if (connected(m, n)) cnt[orbit_of[uf_find((uint32_t)i)]]++;
        }
        size_t *off = malloc((norbits + 1) * sizeof(size_t));
        off[0] = 0;
        for (size_t k = 0; k < norbits; k++) off[k + 1] = off[k] + cnt[k];
        uint64_t *members = malloc(ngraphs * sizeof(uint64_t));
        size_t *fill = calloc(norbits, sizeof(size_t));
        for (uint64_t i = 0; i <= HMASK; i++) {
            if (!H[i]) continue;
            Mat m; decode(H[i] - 1, n, m);
            if (!connected(m, n)) continue;
            uint32_t o = orbit_of[uf_find((uint32_t)i)];
            members[off[o] + fill[o]++] = H[i] - 1;
        }
        #pragma omp parallel for schedule(dynamic, 1)
        for (long long o = 0; o < (long long)norbits; o++) {
            char path[1152];
            snprintf(path, sizeof path, "%s/orbit_%lld", outdir, o);
            FILE *f = fopen(path, "w");
            if (!f) { perror(path); exit(1); }
            uint64_t nb[MAXN * 8];
            for (size_t t = off[o]; t < off[o + 1]; t++) {
                uint64_t g = members[t];
                int c = neighbours(g, n, nb);
                uint64_t uniq[MAXN * 8]; int nu = 0;
                for (int q = 0; q < c; q++) {
                    if (nb[q] < g) continue;
                    int r = 0; while (r < nu && uniq[r] != nb[q]) r++;
                    if (r == nu) uniq[nu++] = nb[q];
                }
                for (int r = 0; r < nu; r++)
                    fprintf(f, "%llu %llu\n", (unsigned long long)g,
                            (unsigned long long)uniq[r]);
            }
            fclose(f);
        }
        free(cnt); free(off); free(members); free(fill);
        fprintf(stderr, "  edge lists: %.2fs -> %s\n", now() - t4, outdir);
    }
    free(orbit_of); free(UF); UF = NULL;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr,
          "usage: %s n d [--outdir-fmt FMT] [--upto] [--loghash K] [--threads T]\n"
          "  --outdir-fmt  printf format for the orbit directory, one %%d for n,\n"
          "                e.g. \"orbits_d3_n%%d_separated\"\n"
          "  --upto        also separate orbits for every 3 <= j <= n\n"
          "  --loghash     log2 of the hash table for the top level; it never\n"
          "                grows, so size it under ~70%% load.  n=7,d=3 -> 23,\n"
          "                n=8,d=3 -> 30.\n", argv[0]);
        return 1;
    }
    int n_target = atoi(argv[1]); D = atoi(argv[2]);
    const char *fmt = NULL;
    bool upto = false;
    int loghash = 23, nthreads = 0;
    for (int i = 3; i < argc; i++) {
        if ((!strcmp(argv[i], "--outdir-fmt") || !strcmp(argv[i], "--outdir"))
            && i + 1 < argc) fmt = argv[++i];
        else if (!strcmp(argv[i], "--upto")) upto = true;
        else if (!strcmp(argv[i], "--loghash") && i + 1 < argc) loghash = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else { fprintf(stderr, "unknown option: %s\n", argv[i]); return 1; }
    }
    if (n_target > MAXN) { fprintf(stderr, "n too large\n"); return 1; }
    if (fmt && upto && !strstr(fmt, "%d")) {
        fprintf(stderr, "--upto needs a %%d in --outdir-fmt, else every level "
                        "would overwrite the same directory\n");
        return 1;
    }
    N = n_target; init_tables();
    if ((long)(n_target * (n_target - 1) / 2) * B > 64) {
        fprintf(stderr, "code exceeds 64 bits\n"); return 1;
    }
#ifdef _OPENMP
    if (nthreads > 0) omp_set_num_threads(nthreads);
    fprintf(stderr, "threads: %d\n", omp_get_max_threads());
#else
    fprintf(stderr, "threads: 1 (built without -fopenmp)\n");
#endif

    double t0 = now();
    char dir[1024];

    uint64_t *cur = malloc(sizeof(uint64_t)); cur[0] = 0;
    size_t cur_n = 1;
    for (int j = 1; j < n_target - 1; j++) {
        int nn = j + 1;
        uint64_t *next = NULL; size_t next_n = 0, next_cap = 0;
        uint64_t *tbl = calloc(1 << 22, sizeof(uint64_t)); uint64_t tm = (1 << 22) - 1;
        Mat m;
        for (size_t s = 0; s < cur_n; s++) {
            decode(cur[s], j, m);
            for (int i = 0; i < nn; i++) m[i][j] = m[j][i] = 0;
            uint64_t total = 1; for (int i = 0; i < j; i++) total *= D;
            for (uint64_t x = 0, col; x < total; x++) {
                col = x;
                for (int i = 0; i < j; i++) { m[i][j] = m[j][i] = col % D; col /= D; }
                uint64_t c = canon(m, nn), k = c + 1, i2 = mix(c) & tm;
                while (tbl[i2] && tbl[i2] != k) i2 = (i2 + 1) & tm;
                if (tbl[i2] == k) continue;
                tbl[i2] = k;
                if (next_n == next_cap) { next_cap = next_cap ? next_cap * 2 : 1024;
                    next = realloc(next, next_cap * sizeof(uint64_t)); }
                next[next_n++] = c;
            }
        }
        free(tbl); free(cur); cur = next; cur_n = next_n;
        fprintf(stderr, "  level n=%d : %llu canonical graphs (%.2fs)\n",
                nn, (unsigned long long)cur_n, now() - t0);

        /* orbits at this level too, if asked */
        if (upto && nn >= 3) {
            int lg = 12; while (((size_t)1 << lg) < cur_n * 2) lg++;
            N = nn; hinit(lg);
            for (size_t s = 0; s < cur_n; s++) { bool nw; hput(cur[s], &nw); }
            const char *out = NULL;
            if (fmt) { snprintf(dir, sizeof dir, fmt, nn); out = dir; }
            link_collect_output(nn, out, 0);
            free(H); H = NULL;
            N = n_target;
        }
    }
    fprintf(stderr, "parents on %d vertices: %llu (%.2fs)\n",
            n_target - 1, (unsigned long long)cur_n, now() - t0);

    double t1 = now();
    N = n_target;
    hinit(loghash);
    uint64_t dpow = 1; for (int i = 0; i < n_target - 1; i++) dpow *= D;
    #pragma omp parallel for schedule(dynamic, 32)
    for (long long s = 0; s < (long long)cur_n; s++) {
        Mat m;
        decode(cur[s], n_target - 1, m);
        for (int i = 0; i < n_target; i++) m[i][n_target-1] = m[n_target-1][i] = 0;
        for (uint64_t x = 0; x < dpow; x++) {
            uint64_t col = x;
            for (int i = 0; i < n_target - 1; i++) {
                m[i][n_target-1] = m[n_target-1][i] = col % D; col /= D;
            }
            bool isnew;
            hput(canon(m, n_target), &isnew);
        }
    }
    free(cur);

    size_t ngen = 0;
    #pragma omp parallel for reduction(+:ngen)
    for (long long i = 0; i <= (long long)HMASK; i++) if (H[i]) ngen++;
    fprintf(stderr, "generation: %llu graphs on %d vertices (%.2fs)\n",
            (unsigned long long)ngen, n_target, now() - t1);
    if (ngen * 10 > (HMASK + 1) * 7)
        fprintf(stderr, "WARNING: table at %.0f%% load, raise --loghash\n",
                100.0 * ngen / (HMASK + 1));

    const char *out = NULL;
    if (fmt) { snprintf(dir, sizeof dir, fmt, n_target); out = dir; }
    link_collect_output(n_target, out, 0);

    fprintf(stderr, "total %.2fs\n", now() - t0);
    return 0;
}
