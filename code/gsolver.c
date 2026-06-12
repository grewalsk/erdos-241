/* gsolver.c — exact g(m) = max size of a B_3 set in Z/mZ (multiset convention).
 *
 * Same incremental scheme as fsolver.c, with all sums reduced mod m.
 * WLOG 0 in A (translation by t shifts triple sums by 3t, a bijection on
 * Z_m, preserving multiset-uniqueness), so DFS starts from A={0} and adds
 * elements of {1..m-1} in increasing order.
 *
 * Every B_3 set in Z_m is Sidon in Z_m (same padding argument as over Z),
 * so maintaining the Sidon invariant is lossless and S2 is faithful.
 *
 * Legality of x (> all current nonzero elements, x not in A):
 *   (P)  (x+a) mod m not in S2 for all a in A; (2x) mod m not in S2
 *   (T)  (S2 + x mod m) disjoint from S3; (2x+a) mod m not in S3; (3x) mod m not in S3
 * Cross-collision coverage identical to the integer case:
 *   x+a+b = 2x+c (mod m) <=> x+c in S2;  x+a+b = 3x <=> 2x in S2;
 *   2x+a = 3x <=> a = x impossible (x fresh);  values in S2+x distinct mod m.
 *
 * Pruning: HONEST ONLY — remaining-candidates bound and the counting bound
 * C(k+2,3) <= m (all triple-sum multisets distinct in Z_m).  We deliberately
 * do NOT use the in-house spectral bound, since g(m) data is meant to test it.
 *
 * Output: CSV  m,gm,ratio6,ratio2,seconds,nodes,witness
 *   ratio6 = gm/(6m)^(1/3)   (counting bound scale)
 *   ratio2 = gm/(2m)^(1/3)   (in-house candidate bound scale)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

typedef uint64_t u64;
#define MAXM 4096
#define MAXWM (MAXM / 64 + 2)
#define MAXK 100

static int M, W, best, kcap;
static int Ael[MAXK];
static u64 S2s[MAXK][MAXWM], S3s[MAXK][MAXWM];
static int bestset[MAXK], bestlen;
static unsigned long long nodes;

static inline int GET(const u64 *B, int i) { return (int)((B[i >> 6] >> (i & 63)) & 1ull); }
static inline void SETB(u64 *B, int i) { B[i >> 6] |= 1ull << (i & 63); }
static inline int MD(int v) { return v >= M ? (v >= 2 * M ? v - 2 * M : v - M) : v; }

static int legal(int s, int x) {
    const u64 *s2 = S2s[s], *s3 = S3s[s];
    for (int i = 0; i < s; i++)
        if (GET(s2, MD(x + Ael[i]))) return 0;
    if (GET(s2, MD(2 * x))) return 0;
    int x2 = MD(2 * x);
    for (int i = 0; i < s; i++)
        if (GET(s3, MD(x2 + Ael[i]))) return 0;
    if (GET(s3, MD(x2 + x))) return 0;
    /* rotate S2 by x and intersect S3: iterate set bits of S2 (sparse: ~s^2/2 bits) */
    for (int w = 0; w < W; w++) {
        u64 v = s2[w];
        while (v) {
            int b = __builtin_ctzll(v);
            v &= v - 1;
            int p = (w << 6) + b;
            if (GET(s3, MD(p + x))) return 0;
        }
    }
    return 1;
}

static void apply_x(int s, int x) {
    u64 *s2 = S2s[s + 1], *s3 = S3s[s + 1];
    const u64 *o2 = S2s[s];
    memcpy(s2, S2s[s], (size_t)W * 8);
    memcpy(s3, S3s[s], (size_t)W * 8);
    for (int w = 0; w < W; w++) {
        u64 v = o2[w];
        while (v) {
            int b = __builtin_ctzll(v);
            v &= v - 1;
            int p = (w << 6) + b;
            SETB(s3, MD(p + x));
        }
    }
    int x2 = MD(2 * x);
    for (int i = 0; i < s; i++) SETB(s3, MD(x2 + Ael[i]));
    SETB(s3, MD(x2 + x));
    for (int i = 0; i < s; i++) SETB(s2, MD(x + Ael[i]));
    SETB(s2, x2);
    Ael[s] = x;
}

static void dfs(int s, int next) {
    nodes++;
    if (s >= kcap) return;
    for (int x = next; x < M; x++) {
        if (s + 1 + (M - 1 - x) <= best) break;
        if (legal(s, x)) {
            apply_x(s, x);
            if (s + 1 > best) {
                best = s + 1;
                bestlen = s + 1;
                memcpy(bestset, Ael, (size_t)(s + 1) * sizeof(int));
            }
            dfs(s + 1, x + 1);
        }
    }
}

int main(int argc, char **argv) {
    int mlo = (argc > 1) ? atoi(argv[1]) : 3;
    int mhi = (argc > 2) ? atoi(argv[2]) : 80;
    double percap = (argc > 3) ? atof(argv[3]) : 1e18;
    printf("m,gm,ratio6,ratio2,seconds,nodes,witness\n");
    for (M = mlo; M <= mhi; M++) {
        W = M / 64 + 1;
        /* counting cap: C(k+2,3) <= m */
        kcap = 0;
        while ((long long)(kcap + 3) * (kcap + 2) * (kcap + 1) / 6 <= M) kcap++;
        memset(S2s[0], 0, (size_t)W * 8);
        memset(S3s[0], 0, (size_t)W * 8);
        /* level 1 = {0}: pair sum 0+0=0, triple sum 0 */
        Ael[0] = 0;
        memcpy(S2s[1], S2s[0], (size_t)W * 8);
        memcpy(S3s[1], S3s[0], (size_t)W * 8);
        SETB(S2s[1], 0);
        SETB(S3s[1], 0);
        best = 1;
        bestlen = 1;
        bestset[0] = 0;
        nodes = 0;
        clock_t t0 = clock();
        dfs(1, 1);
        double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
        printf("%d,%d,%.6f,%.6f,%.3f,%llu,\"", M, best,
               best / cbrt(6.0 * M), best / cbrt(2.0 * M), dt, nodes);
        for (int i = 0; i < bestlen; i++) printf("%d%s", bestset[i], i + 1 < bestlen ? " " : "");
        printf("\"\n");
        fflush(stdout);
        if (dt > percap) {
            fprintf(stderr, "time cap exceeded at m=%d (%.1fs), stopping\n", M, dt);
            break;
        }
    }
    return 0;
}
