/**
 * @file walkcounts.c
 * @brief Molecular Walk and Path Count Descriptors
 *
 * Walk-based topological descriptors:
 * - MWC: Molecular Walk Counts at lengths 1-10
 * - SRW: Self-returning walks (closed walks)
 * - Path counts at various lengths
 * - Weighted path descriptors
 *
 * Uses adjacency matrix powers for efficient walk counting.
 * Total: 36 descriptors
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

#define NUM_WALK_DESCRIPTORS 36
#define MAX_WALK_ATOMS 128
#define MAX_PATH_LENGTH 10

/* ============================================================================
 * Adjacency Matrix Operations
 * ============================================================================ */

/**
 * Build adjacency matrix for heavy atoms only
 */
static void build_adjacency(const molecule_t* mol, double* adj,
                            int* heavy_idx, int* n_heavy) {
    int nh = 0;
    for (int i = 0; i < mol->num_atoms && nh < MAX_WALK_ATOMS; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_idx[nh++] = i;
        }
    }
    *n_heavy = nh;

    /* Initialize to zero */
    memset(adj, 0, nh * nh * sizeof(double));

    /* Create reverse mapping */
    int reverse[512];
    memset(reverse, -1, sizeof(reverse));
    for (int i = 0; i < nh; i++) {
        reverse[heavy_idx[i]] = i;
    }

    /* Fill adjacency from bonds */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int i = reverse[bond->atom1];
        int j = reverse[bond->atom2];
        if (i >= 0 && j >= 0) {
            adj[i * nh + j] = 1.0;
            adj[j * nh + i] = 1.0;
        }
    }
}

/**
 * Matrix multiplication: C = A * B
 */
static void matrix_multiply(const double* A, const double* B, double* C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

/**
 * Sum all elements of matrix
 */
static double matrix_sum(const double* M, int n) {
    double sum = 0.0;
    for (int i = 0; i < n * n; i++) {
        sum += M[i];
    }
    return sum;
}

/**
 * Sum diagonal elements (trace)
 */
static double matrix_trace(const double* M, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += M[i * n + i];
    }
    return sum;
}

/* ============================================================================
 * Path Counting Using DFS
 * ============================================================================ */

/**
 * Count simple paths of exactly length k starting from atom start
 */
static int count_paths_dfs(const double* adj, int n, int start, int length,
                           bool* visited, int current, int depth) {
    if (depth == length) return 1;

    int count = 0;
    for (int next = 0; next < n; next++) {
        if (adj[current * n + next] > 0 && !visited[next]) {
            visited[next] = true;
            count += count_paths_dfs(adj, n, start, length, visited, next, depth + 1);
            visited[next] = false;
        }
    }
    return count;
}

/**
 * Count all simple paths of length k in the molecule
 */
static int count_all_paths(const double* adj, int n, int length) {
    if (length == 0) return n;
    if (length == 1) return (int)(matrix_sum(adj, n) / 2);

    bool* visited = (bool*)alloca(n * sizeof(bool));
    int total = 0;

    for (int start = 0; start < n; start++) {
        memset(visited, 0, n * sizeof(bool));
        visited[start] = true;
        total += count_paths_dfs(adj, n, start, length, visited, start, 0);
    }

    /* Each path counted twice (forward and backward) for length > 0 */
    return total / 2;
}

/* ============================================================================
 * Batch Computation
 * ============================================================================ */

int descriptors_compute_walk_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_WALK_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    int n_heavy;
    int heavy_idx[MAX_WALK_ATOMS];

    /* Check molecule size */
    int n_check = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_check++;
    }
    if (n_check == 0 || n_check > MAX_WALK_ATOMS) return NUM_WALK_DESCRIPTORS;

    /* Build adjacency matrix */
    double* adj = (double*)alloca(MAX_WALK_ATOMS * MAX_WALK_ATOMS * sizeof(double));
    double* power = (double*)alloca(MAX_WALK_ATOMS * MAX_WALK_ATOMS * sizeof(double));
    double* temp = (double*)alloca(MAX_WALK_ATOMS * MAX_WALK_ATOMS * sizeof(double));

    build_adjacency(mol, adj, heavy_idx, &n_heavy);

    if (n_heavy == 0) return NUM_WALK_DESCRIPTORS;

    /* Molecular Walk Counts: A^k gives walks of length k */
    double MWC[MAX_PATH_LENGTH + 1];
    double SRW[MAX_PATH_LENGTH + 1];  /* Self-returning walks (trace of A^k) */

    /* A^0 = I */
    MWC[0] = n_heavy;
    SRW[0] = n_heavy;

    /* A^1 = A */
    memcpy(power, adj, n_heavy * n_heavy * sizeof(double));
    MWC[1] = matrix_sum(power, n_heavy);
    SRW[1] = matrix_trace(power, n_heavy);

    /* A^k = A^(k-1) * A for k >= 2 */
    for (int k = 2; k <= MAX_PATH_LENGTH; k++) {
        matrix_multiply(power, adj, temp, n_heavy);
        memcpy(power, temp, n_heavy * n_heavy * sizeof(double));
        MWC[k] = matrix_sum(power, n_heavy);
        SRW[k] = matrix_trace(power, n_heavy);
    }

    /* Simple path counts (more expensive, limit to small lengths) */
    double PC[8];
    for (int k = 0; k <= 7; k++) {
        if (n_heavy <= 64) {  /* Only compute for smaller molecules */
            PC[k] = count_all_paths(adj, n_heavy, k);
        } else {
            PC[k] = 0.0;
        }
    }

    /* Weighted walk counts */
    double total_mwc = 0.0;
    double weighted_mwc = 0.0;
    for (int k = 1; k <= MAX_PATH_LENGTH; k++) {
        total_mwc += MWC[k];
        weighted_mwc += MWC[k] / k;  /* Weight by 1/length */
    }

    /* Average walk length (weighted mean) */
    double avg_walk = 0.0;
    double walk_sum = 0.0;
    for (int k = 1; k <= MAX_PATH_LENGTH; k++) {
        avg_walk += k * MWC[k];
        walk_sum += MWC[k];
    }
    if (walk_sum > 0) avg_walk /= walk_sum;

    /* Store results */
    int idx = 0;

    /* Walk counts at lengths 1-10 */
    for (int k = 1; k <= 10; k++) {
        values[idx++].d = MWC[k];  /* 0-9: MWC1-MWC10 */
    }

    /* Self-returning walks at lengths 2,4,6,8,10 (odd lengths are 0 for bipartite) */
    values[idx++].d = SRW[2];   /* 10: SRW2 */
    values[idx++].d = SRW[3];   /* 11: SRW3 */
    values[idx++].d = SRW[4];   /* 12: SRW4 */
    values[idx++].d = SRW[5];   /* 13: SRW5 */
    values[idx++].d = SRW[6];   /* 14: SRW6 */

    /* Path counts at lengths 1-7 */
    for (int k = 1; k <= 7; k++) {
        values[idx++].d = PC[k];  /* 15-21: PC1-PC7 */
    }

    /* Summary statistics */
    values[idx++].d = total_mwc;     /* 22: TotalMWC */
    values[idx++].d = weighted_mwc;  /* 23: WeightedMWC */
    values[idx++].d = avg_walk;      /* 24: AvgWalkLength */
    values[idx++].d = MWC[2] > 0 ? MWC[3] / MWC[2] : 0.0;  /* 25: MWC3/MWC2 ratio */
    values[idx++].d = MWC[1] > 0 ? sqrt(MWC[2] / MWC[1]) : 0.0;  /* 26: sqrt(MWC2/MWC1) */

    /* Cyclic walk indicators */
    values[idx++].d = SRW[3] > 0 ? 1.0 : 0.0;   /* 27: Has 3-cycles */
    values[idx++].d = SRW[4] > 0 ? 1.0 : 0.0;   /* 28: Has 4-cycles */
    values[idx++].d = SRW[5] > 0 ? 1.0 : 0.0;   /* 29: Has 5-cycles */
    values[idx++].d = SRW[6] > 0 ? 1.0 : 0.0;   /* 30: Has 6-cycles */

    /* Normalized walks */
    values[idx++].d = n_heavy > 0 ? MWC[2] / n_heavy : 0.0;  /* 31: MWC2 per atom */
    values[idx++].d = n_heavy > 0 ? MWC[3] / n_heavy : 0.0;  /* 32: MWC3 per atom */
    values[idx++].d = n_heavy > 0 ? total_mwc / n_heavy : 0.0;  /* 33: Total MWC per atom */

    /* Log walks (for scale invariance) */
    values[idx++].d = MWC[2] > 0 ? log(MWC[2] + 1) : 0.0;  /* 34: log(MWC2+1) */
    values[idx++].d = total_mwc > 0 ? log(total_mwc + 1) : 0.0;  /* 35: log(TotalMWC+1) */

    return NUM_WALK_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* walk_cached_mol = NULL;
static _Thread_local descriptor_value_t walk_cached_values[NUM_WALK_DESCRIPTORS];

static inline void ensure_walk_computed(const molecule_t* mol) {
    if (walk_cached_mol != mol) {
        descriptors_compute_walk_all(mol, walk_cached_values);
        walk_cached_mol = mol;
    }
}

#define DEFINE_WALK_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_walk_computed(mol); \
    value->d = walk_cached_values[idx].d; \
    return CCHEM_OK; \
}

DEFINE_WALK_FUNC(mwc1, 0)
DEFINE_WALK_FUNC(mwc2, 1)
DEFINE_WALK_FUNC(mwc3, 2)
DEFINE_WALK_FUNC(mwc4, 3)
DEFINE_WALK_FUNC(mwc5, 4)
DEFINE_WALK_FUNC(mwc6, 5)
DEFINE_WALK_FUNC(mwc7, 6)
DEFINE_WALK_FUNC(mwc8, 7)
DEFINE_WALK_FUNC(mwc9, 8)
DEFINE_WALK_FUNC(mwc10, 9)
DEFINE_WALK_FUNC(srw2, 10)
DEFINE_WALK_FUNC(srw3, 11)
DEFINE_WALK_FUNC(srw4, 12)
DEFINE_WALK_FUNC(srw5, 13)
DEFINE_WALK_FUNC(srw6, 14)
DEFINE_WALK_FUNC(pc1, 15)
DEFINE_WALK_FUNC(pc2, 16)
DEFINE_WALK_FUNC(pc3, 17)
DEFINE_WALK_FUNC(pc4, 18)
DEFINE_WALK_FUNC(pc5, 19)
DEFINE_WALK_FUNC(pc6, 20)
DEFINE_WALK_FUNC(pc7, 21)
DEFINE_WALK_FUNC(total_mwc, 22)
DEFINE_WALK_FUNC(weighted_mwc, 23)
DEFINE_WALK_FUNC(avg_walk, 24)
DEFINE_WALK_FUNC(mwc_ratio, 25)
DEFINE_WALK_FUNC(mwc_sqrt, 26)
DEFINE_WALK_FUNC(has_3cyc, 27)
DEFINE_WALK_FUNC(has_4cyc, 28)
DEFINE_WALK_FUNC(has_5cyc, 29)
DEFINE_WALK_FUNC(has_6cyc, 30)
DEFINE_WALK_FUNC(mwc2_per_atom, 31)
DEFINE_WALK_FUNC(mwc3_per_atom, 32)
DEFINE_WALK_FUNC(mwc_per_atom, 33)
DEFINE_WALK_FUNC(log_mwc2, 34)
DEFINE_WALK_FUNC(log_mwc, 35)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_WALK(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_walkcounts(void) {
    /* Molecular walk counts */
    REGISTER_WALK("MWC1", "Molecular walk count length 1", desc_mwc1);
    REGISTER_WALK("MWC2", "Molecular walk count length 2", desc_mwc2);
    REGISTER_WALK("MWC3", "Molecular walk count length 3", desc_mwc3);
    REGISTER_WALK("MWC4", "Molecular walk count length 4", desc_mwc4);
    REGISTER_WALK("MWC5", "Molecular walk count length 5", desc_mwc5);
    REGISTER_WALK("MWC6", "Molecular walk count length 6", desc_mwc6);
    REGISTER_WALK("MWC7", "Molecular walk count length 7", desc_mwc7);
    REGISTER_WALK("MWC8", "Molecular walk count length 8", desc_mwc8);
    REGISTER_WALK("MWC9", "Molecular walk count length 9", desc_mwc9);
    REGISTER_WALK("MWC10", "Molecular walk count length 10", desc_mwc10);

    /* Self-returning walks */
    REGISTER_WALK("SRW2", "Self-returning walks length 2", desc_srw2);
    REGISTER_WALK("SRW3", "Self-returning walks length 3", desc_srw3);
    REGISTER_WALK("SRW4", "Self-returning walks length 4", desc_srw4);
    REGISTER_WALK("SRW5", "Self-returning walks length 5", desc_srw5);
    REGISTER_WALK("SRW6", "Self-returning walks length 6", desc_srw6);

    /* Path counts */
    REGISTER_WALK("PathCount1", "Simple paths of length 1", desc_pc1);
    REGISTER_WALK("PathCount2", "Simple paths of length 2", desc_pc2);
    REGISTER_WALK("PathCount3", "Simple paths of length 3", desc_pc3);
    REGISTER_WALK("PathCount4", "Simple paths of length 4", desc_pc4);
    REGISTER_WALK("PathCount5", "Simple paths of length 5", desc_pc5);
    REGISTER_WALK("PathCount6", "Simple paths of length 6", desc_pc6);
    REGISTER_WALK("PathCount7", "Simple paths of length 7", desc_pc7);

    /* Summary statistics */
    REGISTER_WALK("TotalMWC", "Total molecular walk count", desc_total_mwc);
    REGISTER_WALK("WeightedMWC", "Length-weighted walk count", desc_weighted_mwc);
    REGISTER_WALK("AvgWalkLength", "Average walk length", desc_avg_walk);
    REGISTER_WALK("MWC_Ratio", "MWC3/MWC2 ratio", desc_mwc_ratio);
    REGISTER_WALK("MWC_Sqrt", "sqrt(MWC2/MWC1)", desc_mwc_sqrt);

    /* Cycle indicators */
    REGISTER_WALK("Has3Cycle", "Contains 3-membered cycle", desc_has_3cyc);
    REGISTER_WALK("Has4Cycle", "Contains 4-membered cycle", desc_has_4cyc);
    REGISTER_WALK("Has5Cycle", "Contains 5-membered cycle", desc_has_5cyc);
    REGISTER_WALK("Has6Cycle", "Contains 6-membered cycle", desc_has_6cyc);

    /* Normalized walks */
    REGISTER_WALK("MWC2_PerAtom", "MWC2 per heavy atom", desc_mwc2_per_atom);
    REGISTER_WALK("MWC3_PerAtom", "MWC3 per heavy atom", desc_mwc3_per_atom);
    REGISTER_WALK("MWC_PerAtom", "Total MWC per heavy atom", desc_mwc_per_atom);

    /* Log-transformed */
    REGISTER_WALK("LogMWC2", "Log of MWC2+1", desc_log_mwc2);
    REGISTER_WALK("LogTotalMWC", "Log of total MWC+1", desc_log_mwc);
}
