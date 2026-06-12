/* fsolver.c — exact f(N) for B_3 subsets of {1..N}, multiset convention.
 *
 * Method: for N = 1..NMAX, depth-first branch and bound adding elements in
 * increasing order.  State per level: bitset S2 of pairwise sums a+b (a<=b),
 * bitset S3 of triple sums.  Invariant: the current set is Sidon (B_2) and
 * B_3; requiring Sidon loses nothing because every B_3 set is Sidon
 * (a+b = c+d with {a,b} != {c,d} plus any x in A gives two distinct multiset
 * representations of a+b+x).  Sidon-ness makes S2 a faithful set of pair
 * sums (one bit = one pair), so:
 *
 * x (> all current elements) is addable iff
 *   (P)  x+a not in S2 for all a in A, and 2x not in S2
 *   (T)  (S2 + x) disjoint from S3; 2x+a not in S3 for all a; 3x not in S3
 *
 * (P) also rules out all cross-collisions among the NEW triples:
 *   x+a+b = 2x+c  <=>  x+c in S2;   x+a+b = 3x  <=>  2x in S2.
 * 2x+a = 3x is impossible (a < x); 2x+a = 2x+b forces a=b; values in S2+x
 * are pairwise distinct by Sidon-ness.
 *
 * Pruning: search only for sets beating best = f(N-1) (valid since
 * f(N) <= f(N-1)+1); inside DFS, break the candidate loop at x once
 * size + 1 + f(N-x) <= best  (future elements live in [x+1,N], which
 * translates into [1,N-x], so at most f(N-x) of them; f is nondecreasing,
 * so the break is safe for all larger x).
 *
 * Output (stdout): CSV  N,f(N),ratio,seconds,witness
 * Cross-checked against the verifier-driven brute force (brute_small.py).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

typedef uint64_t u64;
#define MAXN 4000
#define MAXW ((3*MAXN)/64 + 2)
#define MAXK 100

static int N, W, best;
static int ftab[MAXN + 1];
static int Ael[MAXK];
static u64 S2s[MAXK][MAXW], S3s[MAXK][MAXW];
static int bestset[MAXK], bestlen;
static unsigned long long nodes;

static inline int GET(const u64 *B, int i) { return (int)((B[i >> 6] >> (i & 63)) & 1ull); }
static inline void SETB(u64 *B, int i) { B[i >> 6] |= 1ull << (i & 63); }

/* is x addable at level s? */
static int legal(int s, int x) {
    const u64 *s2 = S2s[s], *s3 = S3s[s];
    for (int i = 0; i < s; i++)
        if (GET(s2, x + Ael[i])) return 0;
    if (GET(s2, 2 * x)) return 0;
    for (int i = 0; i < s; i++)
        if (GET(s3, 2 * x + Ael[i])) return 0;
    if (GET(s3, 3 * x)) return 0;
    /* (S2 << x) & S3 == 0 ? */
    int wsh = x >> 6, bsh = x & 63;
    if (bsh) {
        for (int w = W - 1; w >= wsh; w--) {
            u64 v = (s2[w - wsh] << bsh) |
                    ((w - wsh - 1) >= 0 ? (s2[w - wsh - 1] >> (64 - bsh)) : 0ull);
            if (v & s3[w]) return 0;
        }
    } else {
        for (int w = W - 1; w >= wsh; w--)
            if (s2[w - wsh] & s3[w]) return 0;
    }
    return 1;
}

/* build level s+1 from level s by adding x (assumes legal) */
static void apply_x(int s, int x) {
    u64 *s2 = S2s[s + 1], *s3 = S3s[s + 1];
    const u64 *o2 = S2s[s];
    memcpy(s2, S2s[s], (size_t)W * 8);
    memcpy(s3, S3s[s], (size_t)W * 8);
    int wsh = x >> 6, bsh = x & 63;
    if (bsh) {
        for (int w = W - 1; w >= wsh; w--) {
            u64 v = (o2[w - wsh] << bsh) |
                    ((w - wsh - 1) >= 0 ? (o2[w - wsh - 1] >> (64 - bsh)) : 0ull);
            s3[w] |= v;
        }
    } else {
        for (int w = W - 1; w >= wsh; w--) s3[w] |= o2[w - wsh];
    }
    for (int i = 0; i < s; i++) SETB(s3, 2 * x + Ael[i]);
    SETB(s3, 3 * x);
    for (int i = 0; i < s; i++) SETB(s2, x + Ael[i]);
    SETB(s2, 2 * x);
    Ael[s] = x;
}

static void dfs(int s, int next) {
    nodes++;
    for (int x = next; x <= N; x++) {
        if (s + 1 + ftab[N - x] <= best) break;
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
    int nmax = (argc > 1) ? atoi(argv[1]) : 200;
    double percap = (argc > 2) ? atof(argv[2]) : 1e18; /* stop entirely if one N exceeds this many seconds */
    if (nmax > MAXN) nmax = MAXN;
    ftab[0] = 0;
    int prev = 0;
    bestlen = 0;
    printf("N,fN,ratio,seconds,nodes,witness\n");
    for (N = 1; N <= nmax; N++) {
        W = (3 * N) / 64 + 1;
        memset(S2s[0], 0, (size_t)W * 8);
        memset(S3s[0], 0, (size_t)W * 8);
        best = prev;
        nodes = 0;
        clock_t t0 = clock();
        dfs(0, 1);
        double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
        ftab[N] = best;
        printf("%d,%d,%.6f,%.3f,%llu,\"", N, best, best / cbrt((double)N), dt, nodes);
        for (int i = 0; i < bestlen; i++) printf("%d%s", bestset[i], i + 1 < bestlen ? " " : "");
        printf("\"\n");
        fflush(stdout);
        prev = best;
        if (dt > percap) {
            fprintf(stderr, "per-N time cap exceeded at N=%d (%.1fs), stopping\n", N, dt);
            break;
        }
    }
    return 0;
}
