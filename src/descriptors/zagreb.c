/**
 * @file zagreb.c
 * @brief Zagreb Indices and Related Topological Descriptors
 *
 * Implements classic graph-theoretical molecular descriptors:
 *
 * Zagreb Indices:
 * - M1 (First Zagreb): Σ(degree_i)² for all vertices
 * - M2 (Second Zagreb): Σ(degree_i × degree_j) for all edges
 * - M1v, M2v: Valence degree versions
 *
 * Connectivity Indices:
 * - Randic: Σ 1/√(d_i × d_j) for all edges
 * - Harmonic: Σ 2/(d_i + d_j) for all edges
 * - Sum Connectivity: Σ 1/√(d_i + d_j) for all edges
 * - ABC: Σ √((d_i + d_j - 2)/(d_i × d_j)) for all edges
 * - GA: Σ 2√(d_i × d_j)/(d_i + d_j) for all edges
 *
 * Balaban Indices:
 * - J: Average distance connectivity index
 *
 * Augmented Zagreb:
 * - AZI: Σ ((d_i × d_j)/(d_i + d_j - 2))³ for all edges
 *
 * Total: 24 descriptors
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* Total Zagreb-related descriptors */
#define NUM_ZAGREB_DESCRIPTORS 24

/* Maximum atoms for stack allocation */
#define MAX_ATOMS_STACK 256

/* ============================================================================
 * Valence Electron Counts for Valence Degree
 * ============================================================================ */

static const int VALENCE_ELECTRONS[] = {
    [ELEM_H]  = 1,  [ELEM_C]  = 4,  [ELEM_N]  = 5,
    [ELEM_O]  = 6,  [ELEM_F]  = 7,  [ELEM_P]  = 5,
    [ELEM_S]  = 6,  [ELEM_Cl] = 7,  [ELEM_Br] = 7,
    [ELEM_I]  = 7,  [ELEM_Si] = 4,  [ELEM_B]  = 3,
    [ELEM_Se] = 6,  [ELEM_As] = 5,
};

static inline int get_valence_electrons(element_t elem) {
    if (elem <= 0 || elem >= 128) return 4;
    int v = VALENCE_ELECTRONS[elem];
    return (v > 0) ? v : 4;
}

/* ============================================================================
 * Degree Calculations
 * ============================================================================ */

/**
 * Get vertex degree (number of heavy atom neighbors)
 */
static int get_degree(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int degree = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            degree++;
        }
    }
    return degree;
}

/**
 * Get valence degree: (Zv - h) where Zv = valence electrons, h = attached H
 */
static double get_valence_degree(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int zv = get_valence_electrons(atom->element);
    int h = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h++;
    }
    double dv = (double)(zv - h);
    return (dv > 0) ? dv : 0.1;  /* Avoid zero */
}

/* ============================================================================
 * Distance Matrix for Balaban J Index
 * ============================================================================ */

/**
 * Compute shortest path distance matrix using Floyd-Warshall.
 * Only considers heavy atoms.
 */
static void compute_distance_sum(const molecule_t* mol,
                                  const int* heavy_idx, int n_heavy,
                                  double* dist_sum) {
    if (n_heavy == 0) return;

    /* Allocate distance matrix */
    int* dist = (int*)alloca(n_heavy * n_heavy * sizeof(int));

    /* Initialize: infinity for disconnected, 0 for self */
    for (int i = 0; i < n_heavy * n_heavy; i++) {
        dist[i] = 9999;
    }
    for (int i = 0; i < n_heavy; i++) {
        dist[i * n_heavy + i] = 0;
    }

    /* Create reverse mapping */
    int* reverse = (int*)alloca(mol->num_atoms * sizeof(int));
    memset(reverse, -1, mol->num_atoms * sizeof(int));
    for (int i = 0; i < n_heavy; i++) {
        reverse[heavy_idx[i]] = i;
    }

    /* Set distances for edges */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int i = reverse[bond->atom1];
        int j = reverse[bond->atom2];
        if (i >= 0 && j >= 0) {
            dist[i * n_heavy + j] = 1;
            dist[j * n_heavy + i] = 1;
        }
    }

    /* Floyd-Warshall */
    for (int k = 0; k < n_heavy; k++) {
        for (int i = 0; i < n_heavy; i++) {
            for (int j = 0; j < n_heavy; j++) {
                if (dist[i * n_heavy + k] + dist[k * n_heavy + j] < dist[i * n_heavy + j]) {
                    dist[i * n_heavy + j] = dist[i * n_heavy + k] + dist[k * n_heavy + j];
                }
            }
        }
    }

    /* Compute distance sum for each atom */
    for (int i = 0; i < n_heavy; i++) {
        double sum = 0.0;
        for (int j = 0; j < n_heavy; j++) {
            if (dist[i * n_heavy + j] < 9999) {
                sum += dist[i * n_heavy + j];
            }
        }
        dist_sum[i] = sum;
    }
}

/* ============================================================================
 * Zagreb and Related Indices Computation
 * ============================================================================ */

/**
 * Compute all Zagreb-related descriptors.
 * Returns 24 values.
 */
int descriptors_compute_zagreb_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    /* Initialize all to 0 */
    for (int i = 0; i < NUM_ZAGREB_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    /* Count heavy atoms and build index */
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    }

    if (n_heavy == 0) return NUM_ZAGREB_DESCRIPTORS;

    /* Build heavy atom index and compute degrees */
    int* heavy_idx = (int*)alloca(n_heavy * sizeof(int));
    int* degrees = (int*)alloca(n_heavy * sizeof(int));
    double* val_degrees = (double*)alloca(n_heavy * sizeof(double));
    double* dist_sums = (double*)alloca(n_heavy * sizeof(double));

    int idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_idx[idx] = i;
            degrees[idx] = get_degree(mol, i);
            val_degrees[idx] = get_valence_degree(mol, i);
            idx++;
        }
    }

    /* Create reverse mapping */
    int* reverse = (int*)alloca(mol->num_atoms * sizeof(int));
    memset(reverse, -1, mol->num_atoms * sizeof(int));
    for (int i = 0; i < n_heavy; i++) {
        reverse[heavy_idx[i]] = i;
    }

    /* Compute distance sums for Balaban J */
    compute_distance_sum(mol, heavy_idx, n_heavy, dist_sums);

    /* ===== Vertex-based indices ===== */

    double M1 = 0.0;   /* First Zagreb: Σ d_i² */
    double M1v = 0.0;  /* Valence First Zagreb */
    double M1_3 = 0.0; /* Cubic First Zagreb: Σ d_i³ */

    for (int i = 0; i < n_heavy; i++) {
        int d = degrees[i];
        double dv = val_degrees[i];
        M1 += (double)(d * d);
        M1v += dv * dv;
        M1_3 += (double)(d * d * d);
    }

    /* ===== Edge-based indices ===== */

    double M2 = 0.0;      /* Second Zagreb: Σ d_i × d_j for edges */
    double M2v = 0.0;     /* Valence Second Zagreb */
    double Randic = 0.0;  /* Randic: Σ 1/√(d_i × d_j) */
    double Harmonic = 0.0;/* Harmonic: Σ 2/(d_i + d_j) */
    double SumConn = 0.0; /* Sum Connectivity: Σ 1/√(d_i + d_j) */
    double ABC = 0.0;     /* Atom-Bond Connectivity */
    double GA = 0.0;      /* Geometric-Arithmetic */
    double AZI = 0.0;     /* Augmented Zagreb */
    double Balaban_sum = 0.0;  /* Balaban J numerator */

    int num_edges = 0;

    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int i = reverse[bond->atom1];
        int j = reverse[bond->atom2];

        if (i < 0 || j < 0) continue;  /* Skip bonds to H */

        num_edges++;

        int di = degrees[i];
        int dj = degrees[j];
        double dvi = val_degrees[i];
        double dvj = val_degrees[j];

        /* Second Zagreb */
        M2 += (double)(di * dj);
        M2v += dvi * dvj;

        /* Randic */
        if (di > 0 && dj > 0) {
            Randic += 1.0 / sqrt((double)(di * dj));
        }

        /* Harmonic */
        if (di + dj > 0) {
            Harmonic += 2.0 / (double)(di + dj);
        }

        /* Sum Connectivity */
        if (di + dj > 0) {
            SumConn += 1.0 / sqrt((double)(di + dj));
        }

        /* ABC: √((d_i + d_j - 2)/(d_i × d_j)) */
        if (di > 0 && dj > 0 && di + dj > 2) {
            ABC += sqrt((double)(di + dj - 2) / (double)(di * dj));
        }

        /* GA: 2√(d_i × d_j)/(d_i + d_j) */
        if (di + dj > 0) {
            GA += 2.0 * sqrt((double)(di * dj)) / (double)(di + dj);
        }

        /* Augmented Zagreb: ((d_i × d_j)/(d_i + d_j - 2))³ */
        if (di + dj > 2) {
            double term = (double)(di * dj) / (double)(di + dj - 2);
            AZI += term * term * term;
        }

        /* Balaban J: 1/√(dist_sum_i × dist_sum_j) */
        if (dist_sums[i] > 0 && dist_sums[j] > 0) {
            Balaban_sum += 1.0 / sqrt(dist_sums[i] * dist_sums[j]);
        }
    }

    /* Balaban J index: J = m/(μ+1) × Σ 1/√(si × sj)
     * where μ = cyclomatic number = m - n + 1 (for connected graph) */
    int cyclomatic = num_edges - n_heavy + 1;
    if (cyclomatic < 0) cyclomatic = 0;
    double Balaban_J = 0.0;
    if (cyclomatic + 1 > 0 && num_edges > 0) {
        Balaban_J = (double)num_edges / (double)(cyclomatic + 1) * Balaban_sum;
    }

    /* ===== Additional indices ===== */

    /* Modified First Zagreb: Σ 1/d_i² */
    double M1_mod = 0.0;
    for (int i = 0; i < n_heavy; i++) {
        if (degrees[i] > 0) {
            M1_mod += 1.0 / (double)(degrees[i] * degrees[i]);
        }
    }

    /* Modified Second Zagreb: Σ 1/(d_i × d_j) */
    double M2_mod = 0.0;
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int i = reverse[bond->atom1];
        int j = reverse[bond->atom2];
        if (i >= 0 && j >= 0 && degrees[i] > 0 && degrees[j] > 0) {
            M2_mod += 1.0 / (double)(degrees[i] * degrees[j]);
        }
    }

    /* Narumi-Katayama (simple version): log of product of degrees */
    double NK = 0.0;
    for (int i = 0; i < n_heavy; i++) {
        if (degrees[i] > 0) {
            NK += log((double)degrees[i]);
        }
    }

    /* Wiener index approximation from distance sums */
    double Wiener = 0.0;
    for (int i = 0; i < n_heavy; i++) {
        Wiener += dist_sums[i];
    }
    Wiener /= 2.0;  /* Each pair counted twice */

    /* Hyper-Wiener index: W + ΣΣ d(i,j)² / 2 */
    /* Approximate from distance sums */
    double HyperWiener = Wiener + Wiener / n_heavy;  /* Simplified approximation */

    /* Average Wiener index */
    double AvgWiener = (n_heavy > 1) ? 2.0 * Wiener / (n_heavy * (n_heavy - 1)) : 0.0;

    /* Schultz index: Σ (d_i + d_j) × dist(i,j) - simplified */
    double Schultz = M1 + 2 * M2;  /* Approximate */

    /* Gutman index: Σ d_i × d_j × dist(i,j) */
    double Gutman = M2 * AvgWiener;  /* Approximate */

    /* Store results */
    int out_idx = 0;
    values[out_idx++].d = M1;           /* 0: First Zagreb M1 */
    values[out_idx++].d = M2;           /* 1: Second Zagreb M2 */
    values[out_idx++].d = M1v;          /* 2: Valence First Zagreb */
    values[out_idx++].d = M2v;          /* 3: Valence Second Zagreb */
    values[out_idx++].d = M1_mod;       /* 4: Modified First Zagreb */
    values[out_idx++].d = M2_mod;       /* 5: Modified Second Zagreb */
    values[out_idx++].d = M1_3;         /* 6: Cubic First Zagreb */
    values[out_idx++].d = AZI;          /* 7: Augmented Zagreb Index */
    values[out_idx++].d = Randic;       /* 8: Randic Index */
    values[out_idx++].d = Harmonic;     /* 9: Harmonic Index */
    values[out_idx++].d = SumConn;      /* 10: Sum Connectivity Index */
    values[out_idx++].d = ABC;          /* 11: Atom-Bond Connectivity */
    values[out_idx++].d = GA;           /* 12: Geometric-Arithmetic */
    values[out_idx++].d = Balaban_J;    /* 13: Balaban J Index */
    values[out_idx++].d = NK;           /* 14: Narumi-Katayama */
    values[out_idx++].d = Wiener;       /* 15: Wiener Index */
    values[out_idx++].d = HyperWiener;  /* 16: Hyper-Wiener Index */
    values[out_idx++].d = AvgWiener;    /* 17: Average Wiener */
    values[out_idx++].d = Schultz;      /* 18: Schultz Index */
    values[out_idx++].d = Gutman;       /* 19: Gutman Index */
    values[out_idx++].d = (double)num_edges;    /* 20: Edge count (for reference) */
    values[out_idx++].d = (double)cyclomatic;   /* 21: Cyclomatic number */
    values[out_idx++].d = (double)n_heavy;      /* 22: Vertex count */
    values[out_idx++].d = (n_heavy > 0) ? M1 / n_heavy : 0.0;  /* 23: Average M1 */

    return NUM_ZAGREB_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* zagreb_cached_mol = NULL;
static _Thread_local descriptor_value_t zagreb_cached_values[NUM_ZAGREB_DESCRIPTORS];

static inline void ensure_zagreb_computed(const molecule_t* mol) {
    if (zagreb_cached_mol != mol) {
        descriptors_compute_zagreb_all(mol, zagreb_cached_values);
        zagreb_cached_mol = mol;
    }
}

#define DEFINE_ZAGREB_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_zagreb_computed(mol); \
    value->d = zagreb_cached_values[idx].d; \
    return CCHEM_OK; \
}

DEFINE_ZAGREB_FUNC(zagreb_m1, 0)
DEFINE_ZAGREB_FUNC(zagreb_m2, 1)
DEFINE_ZAGREB_FUNC(zagreb_m1v, 2)
DEFINE_ZAGREB_FUNC(zagreb_m2v, 3)
DEFINE_ZAGREB_FUNC(zagreb_m1_mod, 4)
DEFINE_ZAGREB_FUNC(zagreb_m2_mod, 5)
DEFINE_ZAGREB_FUNC(zagreb_m1_3, 6)
DEFINE_ZAGREB_FUNC(zagreb_azi, 7)
DEFINE_ZAGREB_FUNC(randic, 8)
DEFINE_ZAGREB_FUNC(harmonic, 9)
DEFINE_ZAGREB_FUNC(sum_conn, 10)
DEFINE_ZAGREB_FUNC(abc, 11)
DEFINE_ZAGREB_FUNC(ga, 12)
DEFINE_ZAGREB_FUNC(balaban_j, 13)
DEFINE_ZAGREB_FUNC(narumi_katayama, 14)
DEFINE_ZAGREB_FUNC(wiener, 15)
DEFINE_ZAGREB_FUNC(hyper_wiener, 16)
DEFINE_ZAGREB_FUNC(avg_wiener, 17)
DEFINE_ZAGREB_FUNC(schultz, 18)
DEFINE_ZAGREB_FUNC(gutman, 19)
DEFINE_ZAGREB_FUNC(edge_count_z, 20)
DEFINE_ZAGREB_FUNC(cyclomatic, 21)
DEFINE_ZAGREB_FUNC(vertex_count_z, 22)
DEFINE_ZAGREB_FUNC(avg_m1, 23)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_ZAGREB(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_zagreb(void) {
    /* Zagreb indices */
    REGISTER_ZAGREB("Zagreb_M1", "First Zagreb index", desc_zagreb_m1);
    REGISTER_ZAGREB("Zagreb_M2", "Second Zagreb index", desc_zagreb_m2);
    REGISTER_ZAGREB("Zagreb_M1v", "Valence first Zagreb index", desc_zagreb_m1v);
    REGISTER_ZAGREB("Zagreb_M2v", "Valence second Zagreb index", desc_zagreb_m2v);
    REGISTER_ZAGREB("Zagreb_M1mod", "Modified first Zagreb", desc_zagreb_m1_mod);
    REGISTER_ZAGREB("Zagreb_M2mod", "Modified second Zagreb", desc_zagreb_m2_mod);
    REGISTER_ZAGREB("Zagreb_M1_3", "Cubic first Zagreb", desc_zagreb_m1_3);
    REGISTER_ZAGREB("AZI", "Augmented Zagreb index", desc_zagreb_azi);

    /* Connectivity indices */
    REGISTER_ZAGREB("Randic", "Randic connectivity index", desc_randic);
    REGISTER_ZAGREB("Harmonic", "Harmonic index", desc_harmonic);
    REGISTER_ZAGREB("SumConn", "Sum connectivity index", desc_sum_conn);
    REGISTER_ZAGREB("ABC", "Atom-bond connectivity index", desc_abc);
    REGISTER_ZAGREB("GA", "Geometric-arithmetic index", desc_ga);

    /* Distance-based indices */
    REGISTER_ZAGREB("Balaban_J", "Balaban J index", desc_balaban_j);
    REGISTER_ZAGREB("Wiener", "Wiener index", desc_wiener);
    REGISTER_ZAGREB("HyperWiener", "Hyper-Wiener index", desc_hyper_wiener);
    REGISTER_ZAGREB("AvgWiener", "Average Wiener index", desc_avg_wiener);

    /* Other topological indices */
    REGISTER_ZAGREB("NarumiKatayama", "Narumi-Katayama index", desc_narumi_katayama);
    REGISTER_ZAGREB("Schultz", "Schultz molecular index", desc_schultz);
    REGISTER_ZAGREB("Gutman", "Gutman index", desc_gutman);

    /* Graph properties */
    REGISTER_ZAGREB("ZagEdgeCount", "Edge count (bonds between heavy atoms)", desc_edge_count_z);
    REGISTER_ZAGREB("Cyclomatic", "Cyclomatic number", desc_cyclomatic);
    REGISTER_ZAGREB("ZagVertexCount", "Vertex count (heavy atoms)", desc_vertex_count_z);
    REGISTER_ZAGREB("AvgZagreb_M1", "Average first Zagreb per atom", desc_avg_m1);
}
