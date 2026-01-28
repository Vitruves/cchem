/**
 * @file estate_sums.c
 * @brief E-State Index Sums and Statistics
 *
 * Computes actual Kier-Hall electrotopological state (E-State) indices and
 * provides sums by atom type, min/max/mean statistics.
 *
 * E-State index: S_i = I_i + ΔI_i
 * where I_i = intrinsic state, ΔI_i = perturbation from neighbors
 *
 * Descriptors:
 * - Sum of E-States by atom type (C, N, O, S, etc.)
 * - Min/Max/Mean E-State values
 * - E-State range and variance
 * - Hydrogen E-States (attached to heteroatoms)
 *
 * Total: 32 descriptors
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

#define NUM_ESTATE_SUM_DESCRIPTORS 32
#define MAX_ESTATE_ATOMS 512

/* ============================================================================
 * Intrinsic State Values (Kier-Hall)
 * I = (2/N)^2 * delta_v + 1
 * where delta_v = valence electrons - H count
 * ============================================================================ */

static const int PRINCIPAL_QUANTUM[] = {
    [ELEM_H]  = 1, [ELEM_C]  = 2, [ELEM_N]  = 2, [ELEM_O]  = 2,
    [ELEM_F]  = 2, [ELEM_Si] = 3, [ELEM_P]  = 3, [ELEM_S]  = 3,
    [ELEM_Cl] = 3, [ELEM_Br] = 4, [ELEM_I]  = 5, [ELEM_B]  = 2,
    [ELEM_Se] = 4, [ELEM_As] = 4,
};

static const int VALENCE_ELECTRONS[] = {
    [ELEM_H]  = 1, [ELEM_C]  = 4, [ELEM_N]  = 5, [ELEM_O]  = 6,
    [ELEM_F]  = 7, [ELEM_Si] = 4, [ELEM_P]  = 5, [ELEM_S]  = 6,
    [ELEM_Cl] = 7, [ELEM_Br] = 7, [ELEM_I]  = 7, [ELEM_B]  = 3,
    [ELEM_Se] = 6, [ELEM_As] = 5,
};

static inline int get_principal_quantum(element_t elem) {
    if (elem <= 0 || elem >= 128) return 2;
    int n = PRINCIPAL_QUANTUM[elem];
    return n > 0 ? n : 2;
}

static inline int get_valence(element_t elem) {
    if (elem <= 0 || elem >= 128) return 4;
    int v = VALENCE_ELECTRONS[elem];
    return v > 0 ? v : 4;
}

/* ============================================================================
 * E-State Computation
 * ============================================================================ */

/**
 * Count hydrogens attached to atom (explicit + implicit)
 */
static int count_hydrogens(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h++;
    }
    return h;
}

/**
 * Get heavy atom degree
 */
static int get_heavy_degree(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int d = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) d++;
    }
    return d;
}

/**
 * Compute intrinsic state for an atom
 * I = (2/N)^2 * (delta_v + 1)
 */
static double compute_intrinsic_state(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    if (atom->element == ELEM_H) return 0.0;

    int n = get_principal_quantum(atom->element);
    int zv = get_valence(atom->element);
    int h = count_hydrogens(mol, atom_idx);
    int delta = get_heavy_degree(mol, atom_idx);

    /* delta_v = (Zv - h) for sigma electrons */
    double delta_v = (double)(zv - h);
    if (delta <= 0) delta = 1;

    /* I = (2/N)^2 * delta_v / delta + 1 */
    double factor = (2.0 / n) * (2.0 / n);
    return factor * delta_v / delta + 1.0;
}

/**
 * Compute topological distance between two atoms using BFS
 */
static int compute_distance(const molecule_t* mol, int from, int to) {
    if (from == to) return 0;

    bool visited[MAX_ESTATE_ATOMS];
    int dist[MAX_ESTATE_ATOMS];
    int queue[MAX_ESTATE_ATOMS];
    int head = 0, tail = 0;

    memset(visited, 0, mol->num_atoms * sizeof(bool));
    memset(dist, -1, mol->num_atoms * sizeof(int));

    visited[from] = true;
    dist[from] = 0;
    queue[tail++] = from;

    while (head < tail) {
        int curr = queue[head++];
        if (curr == to) return dist[curr];

        const atom_t* atom = &mol->atoms[curr];
        for (int i = 0; i < atom->num_neighbors; i++) {
            int next = atom->neighbors[i];
            if (!visited[next]) {
                visited[next] = true;
                dist[next] = dist[curr] + 1;
                queue[tail++] = next;
            }
        }
    }

    return -1;  /* Not connected */
}

/**
 * Compute all E-State indices for a molecule
 */
static void compute_all_estates(const molecule_t* mol, double* estates) {
    int n = mol->num_atoms;

    /* First pass: compute intrinsic states */
    double* intrinsic = (double*)alloca(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        intrinsic[i] = compute_intrinsic_state(mol, i);
    }

    /* Second pass: add perturbations */
    for (int i = 0; i < n; i++) {
        if (mol->atoms[i].element == ELEM_H) {
            estates[i] = 0.0;
            continue;
        }

        double delta_I = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j || mol->atoms[j].element == ELEM_H) continue;

            int d = compute_distance(mol, i, j);
            if (d > 0 && d < 100) {
                /* ΔI = (I_i - I_j) / (d + 1)^2 */
                delta_I += (intrinsic[i] - intrinsic[j]) / ((d + 1) * (d + 1));
            }
        }

        estates[i] = intrinsic[i] + delta_I;
    }
}

/* ============================================================================
 * Batch Computation
 * ============================================================================ */

int descriptors_compute_estate_sums_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_ESTATE_SUM_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    int n = mol->num_atoms;
    if (n == 0 || n > MAX_ESTATE_ATOMS) return NUM_ESTATE_SUM_DESCRIPTORS;

    /* Compute E-States */
    double* estates = (double*)alloca(n * sizeof(double));
    compute_all_estates(mol, estates);

    /* Collect statistics */
    double sum_C = 0.0, sum_N = 0.0, sum_O = 0.0, sum_S = 0.0;
    double sum_P = 0.0, sum_F = 0.0, sum_Cl = 0.0, sum_Br = 0.0;
    double sum_hal = 0.0, sum_hetero = 0.0;
    double sum_all = 0.0, sum_sq = 0.0;
    double min_e = 1e10, max_e = -1e10;
    int n_heavy = 0;
    int n_C = 0, n_N = 0, n_O = 0, n_S = 0;

    for (int i = 0; i < n; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_H) continue;

        double e = estates[i];
        n_heavy++;
        sum_all += e;
        sum_sq += e * e;
        if (e < min_e) min_e = e;
        if (e > max_e) max_e = e;

        switch (elem) {
            case ELEM_C:  sum_C += e; n_C++; break;
            case ELEM_N:  sum_N += e; sum_hetero += e; n_N++; break;
            case ELEM_O:  sum_O += e; sum_hetero += e; n_O++; break;
            case ELEM_S:  sum_S += e; sum_hetero += e; n_S++; break;
            case ELEM_P:  sum_P += e; sum_hetero += e; break;
            case ELEM_F:  sum_F += e; sum_hal += e; break;
            case ELEM_Cl: sum_Cl += e; sum_hal += e; break;
            case ELEM_Br: sum_Br += e; sum_hal += e; break;
            case ELEM_I:  sum_hal += e; break;
            default: sum_hetero += e; break;
        }
    }

    /* Compute statistics */
    double mean_e = n_heavy > 0 ? sum_all / n_heavy : 0.0;
    double var_e = n_heavy > 1 ? (sum_sq - sum_all * sum_all / n_heavy) / (n_heavy - 1) : 0.0;
    double range_e = max_e - min_e;

    /* E-State of hydrogens on heteroatoms (H-bond related) */
    double sum_h_on_het = 0.0;
    int n_h_on_het = 0;
    for (int i = 0; i < n; i++) {
        element_t elem = mol->atoms[i].element;
        if (elem == ELEM_N || elem == ELEM_O || elem == ELEM_S) {
            int h = count_hydrogens(mol, i);
            /* H E-state approximation: fraction of parent */
            sum_h_on_het += h * estates[i] * 0.3;
            n_h_on_het += h;
        }
    }

    /* Max positive and negative E-states */
    double max_pos = 0.0, max_neg = 0.0;
    for (int i = 0; i < n; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        if (estates[i] > max_pos) max_pos = estates[i];
        if (estates[i] < max_neg) max_neg = estates[i];
    }

    /* Store results */
    int idx = 0;
    values[idx++].d = sum_C;      /* 0: Sum E-State carbons */
    values[idx++].d = sum_N;      /* 1: Sum E-State nitrogens */
    values[idx++].d = sum_O;      /* 2: Sum E-State oxygens */
    values[idx++].d = sum_S;      /* 3: Sum E-State sulfurs */
    values[idx++].d = sum_P;      /* 4: Sum E-State phosphorus */
    values[idx++].d = sum_F;      /* 5: Sum E-State fluorines */
    values[idx++].d = sum_Cl;     /* 6: Sum E-State chlorines */
    values[idx++].d = sum_Br;     /* 7: Sum E-State bromines */
    values[idx++].d = sum_hal;    /* 8: Sum E-State halogens */
    values[idx++].d = sum_hetero; /* 9: Sum E-State heteroatoms */
    values[idx++].d = sum_all;    /* 10: Total E-State sum */

    values[idx++].d = min_e < 1e9 ? min_e : 0.0;  /* 11: Min E-State */
    values[idx++].d = max_e > -1e9 ? max_e : 0.0; /* 12: Max E-State */
    values[idx++].d = mean_e;     /* 13: Mean E-State */
    values[idx++].d = var_e;      /* 14: E-State variance */
    values[idx++].d = range_e > 0 && range_e < 1e9 ? range_e : 0.0; /* 15: E-State range */

    values[idx++].d = max_pos;    /* 16: Max positive E-State */
    values[idx++].d = max_neg;    /* 17: Max negative E-State */
    values[idx++].d = max_pos - max_neg;  /* 18: E-State polarity */

    values[idx++].d = sum_h_on_het;  /* 19: H E-State on heteroatoms */
    values[idx++].d = (double)n_h_on_het;  /* 20: Count H on heteroatoms */

    /* Average E-States by atom type */
    values[idx++].d = n_C > 0 ? sum_C / n_C : 0.0;  /* 21: Avg E-State C */
    values[idx++].d = n_N > 0 ? sum_N / n_N : 0.0;  /* 22: Avg E-State N */
    values[idx++].d = n_O > 0 ? sum_O / n_O : 0.0;  /* 23: Avg E-State O */
    values[idx++].d = n_S > 0 ? sum_S / n_S : 0.0;  /* 24: Avg E-State S */

    /* E-State ratios */
    values[idx++].d = sum_all != 0 ? sum_C / sum_all : 0.0;  /* 25: C E-State fraction */
    values[idx++].d = sum_all != 0 ? sum_hetero / sum_all : 0.0;  /* 26: Hetero E-State fraction */

    /* Standard deviation */
    values[idx++].d = var_e > 0 ? sqrt(var_e) : 0.0;  /* 27: E-State std dev */

    /* Skewness (simplified) */
    double skew = 0.0;
    if (var_e > 0 && n_heavy > 2) {
        for (int i = 0; i < n; i++) {
            if (mol->atoms[i].element != ELEM_H) {
                double z = (estates[i] - mean_e) / sqrt(var_e);
                skew += z * z * z;
            }
        }
        skew /= n_heavy;
    }
    values[idx++].d = skew;  /* 28: E-State skewness */

    /* Absolute sum (ignoring sign) */
    double abs_sum = 0.0;
    for (int i = 0; i < n; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            abs_sum += fabs(estates[i]);
        }
    }
    values[idx++].d = abs_sum;  /* 29: Absolute E-State sum */

    /* E-State index counts */
    values[idx++].d = (double)n_heavy;  /* 30: Heavy atom count */
    values[idx++].d = n_heavy > 0 ? abs_sum / n_heavy : 0.0;  /* 31: Mean abs E-State */

    return NUM_ESTATE_SUM_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* es_cached_mol = NULL;
static _Thread_local descriptor_value_t es_cached_values[NUM_ESTATE_SUM_DESCRIPTORS];

static inline void ensure_es_computed(const molecule_t* mol) {
    if (es_cached_mol != mol) {
        descriptors_compute_estate_sums_all(mol, es_cached_values);
        es_cached_mol = mol;
    }
}

#define DEFINE_ES_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_es_computed(mol); \
    value->d = es_cached_values[idx].d; \
    return CCHEM_OK; \
}

DEFINE_ES_FUNC(es_sum_c, 0)
DEFINE_ES_FUNC(es_sum_n, 1)
DEFINE_ES_FUNC(es_sum_o, 2)
DEFINE_ES_FUNC(es_sum_s, 3)
DEFINE_ES_FUNC(es_sum_p, 4)
DEFINE_ES_FUNC(es_sum_f, 5)
DEFINE_ES_FUNC(es_sum_cl, 6)
DEFINE_ES_FUNC(es_sum_br, 7)
DEFINE_ES_FUNC(es_sum_hal, 8)
DEFINE_ES_FUNC(es_sum_het, 9)
DEFINE_ES_FUNC(es_sum_all, 10)
DEFINE_ES_FUNC(es_min, 11)
DEFINE_ES_FUNC(es_max, 12)
DEFINE_ES_FUNC(es_mean, 13)
DEFINE_ES_FUNC(es_var, 14)
DEFINE_ES_FUNC(es_range, 15)
DEFINE_ES_FUNC(es_max_pos, 16)
DEFINE_ES_FUNC(es_max_neg, 17)
DEFINE_ES_FUNC(es_polarity, 18)
DEFINE_ES_FUNC(es_h_het, 19)
DEFINE_ES_FUNC(es_n_h_het, 20)
DEFINE_ES_FUNC(es_avg_c, 21)
DEFINE_ES_FUNC(es_avg_n, 22)
DEFINE_ES_FUNC(es_avg_o, 23)
DEFINE_ES_FUNC(es_avg_s, 24)
DEFINE_ES_FUNC(es_frac_c, 25)
DEFINE_ES_FUNC(es_frac_het, 26)
DEFINE_ES_FUNC(es_std, 27)
DEFINE_ES_FUNC(es_skew, 28)
DEFINE_ES_FUNC(es_abs_sum, 29)
DEFINE_ES_FUNC(es_n_heavy, 30)
DEFINE_ES_FUNC(es_mean_abs, 31)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_ES(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_ELECTRONIC; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_estate_sums(void) {
    /* E-State sums by element */
    REGISTER_ES("EState_Sum_C", "Sum of E-State indices for C", desc_es_sum_c);
    REGISTER_ES("EState_Sum_N", "Sum of E-State indices for N", desc_es_sum_n);
    REGISTER_ES("EState_Sum_O", "Sum of E-State indices for O", desc_es_sum_o);
    REGISTER_ES("EState_Sum_S", "Sum of E-State indices for S", desc_es_sum_s);
    REGISTER_ES("EState_Sum_P", "Sum of E-State indices for P", desc_es_sum_p);
    REGISTER_ES("EState_Sum_F", "Sum of E-State indices for F", desc_es_sum_f);
    REGISTER_ES("EState_Sum_Cl", "Sum of E-State indices for Cl", desc_es_sum_cl);
    REGISTER_ES("EState_Sum_Br", "Sum of E-State indices for Br", desc_es_sum_br);
    REGISTER_ES("EState_Sum_Hal", "Sum of E-State halogen", desc_es_sum_hal);
    REGISTER_ES("EState_Sum_Het", "Sum of E-State heteroatoms", desc_es_sum_het);
    REGISTER_ES("EState_Sum_All", "Total E-State sum", desc_es_sum_all);

    /* Statistics */
    REGISTER_ES("EState_Min", "Minimum E-State index", desc_es_min);
    REGISTER_ES("EState_Max", "Maximum E-State index", desc_es_max);
    REGISTER_ES("EState_Mean", "Mean E-State index", desc_es_mean);
    REGISTER_ES("EState_Var", "E-State variance", desc_es_var);
    REGISTER_ES("EState_Range", "E-State range (max-min)", desc_es_range);
    REGISTER_ES("EState_MaxPos", "Max positive E-State", desc_es_max_pos);
    REGISTER_ES("EState_MaxNeg", "Max negative E-State", desc_es_max_neg);
    REGISTER_ES("EState_Polarity", "E-State polarity", desc_es_polarity);

    /* Hydrogen E-states */
    REGISTER_ES("EState_H_Het", "E-State of H on heteroatoms", desc_es_h_het);
    REGISTER_ES("EState_nH_Het", "Count H on heteroatoms", desc_es_n_h_het);

    /* Averages */
    REGISTER_ES("EState_Avg_C", "Average E-State for C", desc_es_avg_c);
    REGISTER_ES("EState_Avg_N", "Average E-State for N", desc_es_avg_n);
    REGISTER_ES("EState_Avg_O", "Average E-State for O", desc_es_avg_o);
    REGISTER_ES("EState_Avg_S", "Average E-State for S", desc_es_avg_s);

    /* Fractions and derived */
    REGISTER_ES("EState_Frac_C", "Carbon E-State fraction", desc_es_frac_c);
    REGISTER_ES("EState_Frac_Het", "Heteroatom E-State fraction", desc_es_frac_het);
    REGISTER_ES("EState_Std", "E-State std deviation", desc_es_std);
    REGISTER_ES("EState_Skew", "E-State skewness", desc_es_skew);
    REGISTER_ES("EState_AbsSum", "Absolute E-State sum", desc_es_abs_sum);
    REGISTER_ES("EState_nHeavy", "Heavy atom count", desc_es_n_heavy);
    REGISTER_ES("EState_MeanAbs", "Mean absolute E-State", desc_es_mean_abs);
}
