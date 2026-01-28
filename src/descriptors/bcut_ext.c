/**
 * @file bcut_ext.c
 * @brief Extended BCUT Descriptors - Additional Burden eigenvalue descriptors
 *
 * Extension of BCUT descriptors with additional atomic properties:
 * - Electron Affinity (ea)
 * - Covalent Radius (r)
 * - Valence Electrons (v)
 * - VdW Volume (vdw)
 * - Oxidation State (ox)
 * - Lone Pairs (lp)
 * - Aromatic Flag (ar)
 * - Ring Count (rc)
 *
 * Total: 8 properties × 6 eigenvalues (3 high, 3 low) = 48 descriptors
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* Number of eigenvalues to extract per property (3 high + 3 low) */
#define NUM_EIGENVALUES 6

/* Number of properties in this extension */
#define NUM_BCUT_EXT_PROPERTIES 8

/* Total BCUT extension descriptors */
#define NUM_BCUT_EXT_DESCRIPTORS (NUM_BCUT_EXT_PROPERTIES * NUM_EIGENVALUES)

/* Maximum atoms for stack allocation */
#define MAX_ATOMS_STACK 128

/* ============================================================================
 * Atomic Property Tables
 * ============================================================================ */

/* Electron Affinity (eV) - Mulliken values */
static const double BCUT_EA[] = {
    [ELEM_H]  = 0.754,  [ELEM_C]  = 1.262,  [ELEM_N]  = -0.07,
    [ELEM_O]  = 1.461,  [ELEM_F]  = 3.401,  [ELEM_P]  = 0.746,
    [ELEM_S]  = 2.077,  [ELEM_Cl] = 3.612,  [ELEM_Br] = 3.364,
    [ELEM_I]  = 3.059,  [ELEM_Si] = 1.389,  [ELEM_B]  = 0.277,
    [ELEM_Se] = 2.021,  [ELEM_As] = 0.804,
};

/* Covalent Radius (Å) */
static const double BCUT_COVALENT_R[] = {
    [ELEM_H]  = 0.31,  [ELEM_C]  = 0.77,  [ELEM_N]  = 0.71,
    [ELEM_O]  = 0.66,  [ELEM_F]  = 0.57,  [ELEM_P]  = 1.07,
    [ELEM_S]  = 1.05,  [ELEM_Cl] = 1.02,  [ELEM_Br] = 1.20,
    [ELEM_I]  = 1.39,  [ELEM_Si] = 1.11,  [ELEM_B]  = 0.84,
    [ELEM_Se] = 1.20,  [ELEM_As] = 1.19,
};

/* Number of valence electrons */
static const double BCUT_VALENCE[] = {
    [ELEM_H]  = 1,  [ELEM_C]  = 4,  [ELEM_N]  = 5,
    [ELEM_O]  = 6,  [ELEM_F]  = 7,  [ELEM_P]  = 5,
    [ELEM_S]  = 6,  [ELEM_Cl] = 7,  [ELEM_Br] = 7,
    [ELEM_I]  = 7,  [ELEM_Si] = 4,  [ELEM_B]  = 3,
    [ELEM_Se] = 6,  [ELEM_As] = 5,
};

/* Van der Waals Volume (Å³) - Bondi values */
static const double BCUT_VDW_VOL[] = {
    [ELEM_H]  = 7.24,   [ELEM_C]  = 20.58,  [ELEM_N]  = 15.60,
    [ELEM_O]  = 14.71,  [ELEM_F]  = 13.31,  [ELEM_P]  = 24.43,
    [ELEM_S]  = 24.43,  [ELEM_Cl] = 22.45,  [ELEM_Br] = 26.52,
    [ELEM_I]  = 32.52,  [ELEM_Si] = 38.79,  [ELEM_B]  = 17.88,
    [ELEM_Se] = 28.73,  [ELEM_As] = 26.52,
};

/* Common oxidation states (most common positive) */
static const double BCUT_OXIDATION[] = {
    [ELEM_H]  = 1,   [ELEM_C]  = 4,   [ELEM_N]  = 3,
    [ELEM_O]  = -2,  [ELEM_F]  = -1,  [ELEM_P]  = 5,
    [ELEM_S]  = 6,   [ELEM_Cl] = -1,  [ELEM_Br] = -1,
    [ELEM_I]  = -1,  [ELEM_Si] = 4,   [ELEM_B]  = 3,
    [ELEM_Se] = 4,   [ELEM_As] = 5,
};

/* Expected lone pairs at neutral state */
static const double BCUT_LONE_PAIRS[] = {
    [ELEM_H]  = 0,  [ELEM_C]  = 0,  [ELEM_N]  = 1,
    [ELEM_O]  = 2,  [ELEM_F]  = 3,  [ELEM_P]  = 1,
    [ELEM_S]  = 2,  [ELEM_Cl] = 3,  [ELEM_Br] = 3,
    [ELEM_I]  = 3,  [ELEM_Si] = 0,  [ELEM_B]  = 0,
    [ELEM_Se] = 2,  [ELEM_As] = 1,
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_property(element_t elem, const double* table, double default_val) {
    if (elem <= 0 || elem >= 128) return default_val;
    double v = table[elem];
    return (v != 0.0 || elem == ELEM_C) ? v : default_val;
}

/* Get aromatic flag (1.0 if aromatic, 0.0 otherwise) */
static double get_aromatic_flag(const molecule_t* mol, int atom_idx) {
    /* Check if atom is in any aromatic ring */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if ((bond->atom1 == atom_idx || bond->atom2 == atom_idx) && bond->aromatic) {
            return 1.0;
        }
    }
    return 0.0;
}

/* Count rings containing this atom */
static double get_ring_count(const molecule_t* mol, int atom_idx) {
    int ring_count = 0;
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        for (int i = 0; i < ring->size; i++) {
            if (ring->atoms[i] == atom_idx) {
                ring_count++;
                break;
            }
        }
    }
    return (double)ring_count;
}

/* ============================================================================
 * Jacobi Eigenvalue Algorithm for Symmetric Matrices
 * ============================================================================ */

#define JACOBI_MAX_ITER 50
#define JACOBI_EPS 1e-10

static int jacobi_eigenvalues(double* A, int n, double* eigenvalues) {
    if (n <= 0) return 0;
    if (n == 1) {
        eigenvalues[0] = A[0];
        return 1;
    }

    double* b = (double*)alloca(n * sizeof(double));
    double* z = (double*)alloca(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        eigenvalues[i] = A[i * n + i];
        b[i] = eigenvalues[i];
        z[i] = 0.0;
    }

    for (int iter = 0; iter < JACOBI_MAX_ITER; iter++) {
        double sm = 0.0;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                sm += fabs(A[i * n + j]);
            }
        }

        if (sm < JACOBI_EPS) break;

        double thresh = (iter < 4) ? 0.2 * sm / (n * n) : 0.0;

        for (int ip = 0; ip < n - 1; ip++) {
            for (int iq = ip + 1; iq < n; iq++) {
                double g = 100.0 * fabs(A[ip * n + iq]);

                if (iter > 4 &&
                    fabs(eigenvalues[ip]) + g == fabs(eigenvalues[ip]) &&
                    fabs(eigenvalues[iq]) + g == fabs(eigenvalues[iq])) {
                    A[ip * n + iq] = 0.0;
                } else if (fabs(A[ip * n + iq]) > thresh) {
                    double h = eigenvalues[iq] - eigenvalues[ip];
                    double t;

                    if (fabs(h) + g == fabs(h)) {
                        t = A[ip * n + iq] / h;
                    } else {
                        double theta = 0.5 * h / A[ip * n + iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }

                    double c = 1.0 / sqrt(1.0 + t * t);
                    double s = t * c;
                    double tau = s / (1.0 + c);
                    h = t * A[ip * n + iq];

                    z[ip] -= h;
                    z[iq] += h;
                    eigenvalues[ip] -= h;
                    eigenvalues[iq] += h;
                    A[ip * n + iq] = 0.0;

                    for (int j = 0; j < ip; j++) {
                        g = A[j * n + ip];
                        h = A[j * n + iq];
                        A[j * n + ip] = g - s * (h + g * tau);
                        A[j * n + iq] = h + s * (g - h * tau);
                    }
                    for (int j = ip + 1; j < iq; j++) {
                        g = A[ip * n + j];
                        h = A[j * n + iq];
                        A[ip * n + j] = g - s * (h + g * tau);
                        A[j * n + iq] = h + s * (g - h * tau);
                    }
                    for (int j = iq + 1; j < n; j++) {
                        g = A[ip * n + j];
                        h = A[iq * n + j];
                        A[ip * n + j] = g - s * (h + g * tau);
                        A[iq * n + j] = h + s * (g - h * tau);
                    }
                }
            }
        }

        for (int i = 0; i < n; i++) {
            b[i] += z[i];
            eigenvalues[i] = b[i];
            z[i] = 0.0;
        }
    }

    /* Sort eigenvalues descending */
    for (int i = 0; i < n - 1; i++) {
        int max_idx = i;
        for (int j = i + 1; j < n; j++) {
            if (eigenvalues[j] > eigenvalues[max_idx]) max_idx = j;
        }
        if (max_idx != i) {
            double tmp = eigenvalues[i];
            eigenvalues[i] = eigenvalues[max_idx];
            eigenvalues[max_idx] = tmp;
        }
    }

    return n;
}

/* ============================================================================
 * BCUT Extension Computation
 * ============================================================================ */

static void build_burden_matrix(const molecule_t* mol,
                                 const int* heavy_idx, int n_heavy,
                                 const double* props,
                                 double* matrix) {
    for (int i = 0; i < n_heavy; i++) {
        for (int j = 0; j < n_heavy; j++) {
            if (i == j) {
                matrix[i * n_heavy + j] = props[i];
            } else {
                matrix[i * n_heavy + j] = 0.001;
            }
        }
    }

    int* reverse_map = (int*)alloca(mol->num_atoms * sizeof(int));
    memset(reverse_map, -1, mol->num_atoms * sizeof(int));
    for (int i = 0; i < n_heavy; i++) {
        reverse_map[heavy_idx[i]] = i;
    }

    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int i = reverse_map[bond->atom1];
        int j = reverse_map[bond->atom2];

        if (i >= 0 && j >= 0) {
            double bond_order = 1.0;
            switch (bond->type) {
                case BOND_SINGLE: bond_order = 1.0; break;
                case BOND_DOUBLE: bond_order = 2.0; break;
                case BOND_TRIPLE: bond_order = 3.0; break;
                default: bond_order = 1.0; break;
            }
            if (bond->aromatic) bond_order = 1.5;
            double val = 0.1 * sqrt(bond_order);
            matrix[i * n_heavy + j] = val;
            matrix[j * n_heavy + i] = val;
        }
    }
}

/* Store eigenvalues for one property: 3 highest + 3 lowest */
static inline void store_eigenvalues(double* eigenvalues, int n_heavy,
                                      descriptor_value_t* out, int* idx) {
    out[(*idx)++].d = eigenvalues[0];
    out[(*idx)++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    out[(*idx)++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    out[(*idx)++].d = eigenvalues[n_heavy - 1];
    out[(*idx)++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    out[(*idx)++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;
}

int descriptors_compute_bcut_ext_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_BCUT_EXT_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    }

    if (n_heavy <= 1) return NUM_BCUT_EXT_DESCRIPTORS;

    bool heap_alloc = (n_heavy > MAX_ATOMS_STACK);
    int* heavy_idx;
    double* props;
    double* matrix;
    double* eigenvalues;

    if (heap_alloc) {
        heavy_idx = (int*)malloc(n_heavy * sizeof(int));
        props = (double*)malloc(n_heavy * sizeof(double));
        matrix = (double*)malloc(n_heavy * n_heavy * sizeof(double));
        eigenvalues = (double*)malloc(n_heavy * sizeof(double));
        if (!heavy_idx || !props || !matrix || !eigenvalues) {
            free(heavy_idx); free(props); free(matrix); free(eigenvalues);
            return -1;
        }
    } else {
        heavy_idx = (int*)alloca(n_heavy * sizeof(int));
        props = (double*)alloca(n_heavy * sizeof(double));
        matrix = (double*)alloca(n_heavy * n_heavy * sizeof(double));
        eigenvalues = (double*)alloca(n_heavy * sizeof(double));
    }

    int idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_idx[idx++] = i;
        }
    }

    int out_idx = 0;

    /* Property 0: Electron Affinity */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_EA, 1.0);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    /* Property 1: Covalent Radius */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_COVALENT_R, 0.77);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    /* Property 2: Valence Electrons */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_VALENCE, 4.0);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    /* Property 3: VdW Volume */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_VDW_VOL, 20.58);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    /* Property 4: Oxidation State */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_OXIDATION, 0.0);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    /* Property 5: Lone Pairs */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_LONE_PAIRS, 0.0);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    /* Property 6: Aromatic Flag */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_aromatic_flag(mol, heavy_idx[i]);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    /* Property 7: Ring Count */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_ring_count(mol, heavy_idx[i]);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    store_eigenvalues(eigenvalues, n_heavy, values, &out_idx);

    if (heap_alloc) {
        free(heavy_idx); free(props); free(matrix); free(eigenvalues);
    }

    return NUM_BCUT_EXT_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* bcut_ext_cached_mol = NULL;
static _Thread_local uint64_t bcut_ext_cached_gen = 0;
static _Thread_local descriptor_value_t bcut_ext_cached_values[NUM_BCUT_EXT_DESCRIPTORS];

static inline void ensure_bcut_ext_computed(const molecule_t* mol) {
    uint64_t current_gen = descriptor_cache_generation();
    if (bcut_ext_cached_mol != mol || bcut_ext_cached_gen != current_gen) {
        descriptors_compute_bcut_ext_all(mol, bcut_ext_cached_values);
        bcut_ext_cached_mol = mol;
        bcut_ext_cached_gen = current_gen;
    }
}

#define DEFINE_BCUT_EXT_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_bcut_ext_computed(mol); \
    value->d = bcut_ext_cached_values[idx].d; \
    return CCHEM_OK; \
}

/* Electron Affinity BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_ea_hi1, 0)
DEFINE_BCUT_EXT_FUNC(bcut_ea_hi2, 1)
DEFINE_BCUT_EXT_FUNC(bcut_ea_hi3, 2)
DEFINE_BCUT_EXT_FUNC(bcut_ea_lo1, 3)
DEFINE_BCUT_EXT_FUNC(bcut_ea_lo2, 4)
DEFINE_BCUT_EXT_FUNC(bcut_ea_lo3, 5)

/* Covalent Radius BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_r_hi1, 6)
DEFINE_BCUT_EXT_FUNC(bcut_r_hi2, 7)
DEFINE_BCUT_EXT_FUNC(bcut_r_hi3, 8)
DEFINE_BCUT_EXT_FUNC(bcut_r_lo1, 9)
DEFINE_BCUT_EXT_FUNC(bcut_r_lo2, 10)
DEFINE_BCUT_EXT_FUNC(bcut_r_lo3, 11)

/* Valence Electron BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_v_hi1, 12)
DEFINE_BCUT_EXT_FUNC(bcut_v_hi2, 13)
DEFINE_BCUT_EXT_FUNC(bcut_v_hi3, 14)
DEFINE_BCUT_EXT_FUNC(bcut_v_lo1, 15)
DEFINE_BCUT_EXT_FUNC(bcut_v_lo2, 16)
DEFINE_BCUT_EXT_FUNC(bcut_v_lo3, 17)

/* VdW Volume BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_vdw_hi1, 18)
DEFINE_BCUT_EXT_FUNC(bcut_vdw_hi2, 19)
DEFINE_BCUT_EXT_FUNC(bcut_vdw_hi3, 20)
DEFINE_BCUT_EXT_FUNC(bcut_vdw_lo1, 21)
DEFINE_BCUT_EXT_FUNC(bcut_vdw_lo2, 22)
DEFINE_BCUT_EXT_FUNC(bcut_vdw_lo3, 23)

/* Oxidation State BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_ox_hi1, 24)
DEFINE_BCUT_EXT_FUNC(bcut_ox_hi2, 25)
DEFINE_BCUT_EXT_FUNC(bcut_ox_hi3, 26)
DEFINE_BCUT_EXT_FUNC(bcut_ox_lo1, 27)
DEFINE_BCUT_EXT_FUNC(bcut_ox_lo2, 28)
DEFINE_BCUT_EXT_FUNC(bcut_ox_lo3, 29)

/* Lone Pair BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_lp_hi1, 30)
DEFINE_BCUT_EXT_FUNC(bcut_lp_hi2, 31)
DEFINE_BCUT_EXT_FUNC(bcut_lp_hi3, 32)
DEFINE_BCUT_EXT_FUNC(bcut_lp_lo1, 33)
DEFINE_BCUT_EXT_FUNC(bcut_lp_lo2, 34)
DEFINE_BCUT_EXT_FUNC(bcut_lp_lo3, 35)

/* Aromatic Flag BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_ar_hi1, 36)
DEFINE_BCUT_EXT_FUNC(bcut_ar_hi2, 37)
DEFINE_BCUT_EXT_FUNC(bcut_ar_hi3, 38)
DEFINE_BCUT_EXT_FUNC(bcut_ar_lo1, 39)
DEFINE_BCUT_EXT_FUNC(bcut_ar_lo2, 40)
DEFINE_BCUT_EXT_FUNC(bcut_ar_lo3, 41)

/* Ring Count BCUT */
DEFINE_BCUT_EXT_FUNC(bcut_rc_hi1, 42)
DEFINE_BCUT_EXT_FUNC(bcut_rc_hi2, 43)
DEFINE_BCUT_EXT_FUNC(bcut_rc_hi3, 44)
DEFINE_BCUT_EXT_FUNC(bcut_rc_lo1, 45)
DEFINE_BCUT_EXT_FUNC(bcut_rc_lo2, 46)
DEFINE_BCUT_EXT_FUNC(bcut_rc_lo3, 47)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_BCUT_EXT(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_bcut_ext(void) {
    /* Electron Affinity BCUT */
    REGISTER_BCUT_EXT("BCUTea_hi1", "EA BCUT highest eigenvalue", desc_bcut_ea_hi1);
    REGISTER_BCUT_EXT("BCUTea_hi2", "EA BCUT 2nd highest eigenvalue", desc_bcut_ea_hi2);
    REGISTER_BCUT_EXT("BCUTea_hi3", "EA BCUT 3rd highest eigenvalue", desc_bcut_ea_hi3);
    REGISTER_BCUT_EXT("BCUTea_lo1", "EA BCUT lowest eigenvalue", desc_bcut_ea_lo1);
    REGISTER_BCUT_EXT("BCUTea_lo2", "EA BCUT 2nd lowest eigenvalue", desc_bcut_ea_lo2);
    REGISTER_BCUT_EXT("BCUTea_lo3", "EA BCUT 3rd lowest eigenvalue", desc_bcut_ea_lo3);

    /* Covalent Radius BCUT */
    REGISTER_BCUT_EXT("BCUTr_hi1", "Radius BCUT highest eigenvalue", desc_bcut_r_hi1);
    REGISTER_BCUT_EXT("BCUTr_hi2", "Radius BCUT 2nd highest eigenvalue", desc_bcut_r_hi2);
    REGISTER_BCUT_EXT("BCUTr_hi3", "Radius BCUT 3rd highest eigenvalue", desc_bcut_r_hi3);
    REGISTER_BCUT_EXT("BCUTr_lo1", "Radius BCUT lowest eigenvalue", desc_bcut_r_lo1);
    REGISTER_BCUT_EXT("BCUTr_lo2", "Radius BCUT 2nd lowest eigenvalue", desc_bcut_r_lo2);
    REGISTER_BCUT_EXT("BCUTr_lo3", "Radius BCUT 3rd lowest eigenvalue", desc_bcut_r_lo3);

    /* Valence Electron BCUT */
    REGISTER_BCUT_EXT("BCUTv_hi1", "Valence BCUT highest eigenvalue", desc_bcut_v_hi1);
    REGISTER_BCUT_EXT("BCUTv_hi2", "Valence BCUT 2nd highest eigenvalue", desc_bcut_v_hi2);
    REGISTER_BCUT_EXT("BCUTv_hi3", "Valence BCUT 3rd highest eigenvalue", desc_bcut_v_hi3);
    REGISTER_BCUT_EXT("BCUTv_lo1", "Valence BCUT lowest eigenvalue", desc_bcut_v_lo1);
    REGISTER_BCUT_EXT("BCUTv_lo2", "Valence BCUT 2nd lowest eigenvalue", desc_bcut_v_lo2);
    REGISTER_BCUT_EXT("BCUTv_lo3", "Valence BCUT 3rd lowest eigenvalue", desc_bcut_v_lo3);

    /* VdW Volume BCUT */
    REGISTER_BCUT_EXT("BCUTvdw_hi1", "VdW Volume BCUT highest eigenvalue", desc_bcut_vdw_hi1);
    REGISTER_BCUT_EXT("BCUTvdw_hi2", "VdW Volume BCUT 2nd highest eigenvalue", desc_bcut_vdw_hi2);
    REGISTER_BCUT_EXT("BCUTvdw_hi3", "VdW Volume BCUT 3rd highest eigenvalue", desc_bcut_vdw_hi3);
    REGISTER_BCUT_EXT("BCUTvdw_lo1", "VdW Volume BCUT lowest eigenvalue", desc_bcut_vdw_lo1);
    REGISTER_BCUT_EXT("BCUTvdw_lo2", "VdW Volume BCUT 2nd lowest eigenvalue", desc_bcut_vdw_lo2);
    REGISTER_BCUT_EXT("BCUTvdw_lo3", "VdW Volume BCUT 3rd lowest eigenvalue", desc_bcut_vdw_lo3);

    /* Oxidation State BCUT */
    REGISTER_BCUT_EXT("BCUTox_hi1", "Oxidation BCUT highest eigenvalue", desc_bcut_ox_hi1);
    REGISTER_BCUT_EXT("BCUTox_hi2", "Oxidation BCUT 2nd highest eigenvalue", desc_bcut_ox_hi2);
    REGISTER_BCUT_EXT("BCUTox_hi3", "Oxidation BCUT 3rd highest eigenvalue", desc_bcut_ox_hi3);
    REGISTER_BCUT_EXT("BCUTox_lo1", "Oxidation BCUT lowest eigenvalue", desc_bcut_ox_lo1);
    REGISTER_BCUT_EXT("BCUTox_lo2", "Oxidation BCUT 2nd lowest eigenvalue", desc_bcut_ox_lo2);
    REGISTER_BCUT_EXT("BCUTox_lo3", "Oxidation BCUT 3rd lowest eigenvalue", desc_bcut_ox_lo3);

    /* Lone Pair BCUT */
    REGISTER_BCUT_EXT("BCUTlp_hi1", "Lone Pair BCUT highest eigenvalue", desc_bcut_lp_hi1);
    REGISTER_BCUT_EXT("BCUTlp_hi2", "Lone Pair BCUT 2nd highest eigenvalue", desc_bcut_lp_hi2);
    REGISTER_BCUT_EXT("BCUTlp_hi3", "Lone Pair BCUT 3rd highest eigenvalue", desc_bcut_lp_hi3);
    REGISTER_BCUT_EXT("BCUTlp_lo1", "Lone Pair BCUT lowest eigenvalue", desc_bcut_lp_lo1);
    REGISTER_BCUT_EXT("BCUTlp_lo2", "Lone Pair BCUT 2nd lowest eigenvalue", desc_bcut_lp_lo2);
    REGISTER_BCUT_EXT("BCUTlp_lo3", "Lone Pair BCUT 3rd lowest eigenvalue", desc_bcut_lp_lo3);

    /* Aromatic Flag BCUT */
    REGISTER_BCUT_EXT("BCUTar_hi1", "Aromatic BCUT highest eigenvalue", desc_bcut_ar_hi1);
    REGISTER_BCUT_EXT("BCUTar_hi2", "Aromatic BCUT 2nd highest eigenvalue", desc_bcut_ar_hi2);
    REGISTER_BCUT_EXT("BCUTar_hi3", "Aromatic BCUT 3rd highest eigenvalue", desc_bcut_ar_hi3);
    REGISTER_BCUT_EXT("BCUTar_lo1", "Aromatic BCUT lowest eigenvalue", desc_bcut_ar_lo1);
    REGISTER_BCUT_EXT("BCUTar_lo2", "Aromatic BCUT 2nd lowest eigenvalue", desc_bcut_ar_lo2);
    REGISTER_BCUT_EXT("BCUTar_lo3", "Aromatic BCUT 3rd lowest eigenvalue", desc_bcut_ar_lo3);

    /* Ring Count BCUT */
    REGISTER_BCUT_EXT("BCUTrc_hi1", "Ring Count BCUT highest eigenvalue", desc_bcut_rc_hi1);
    REGISTER_BCUT_EXT("BCUTrc_hi2", "Ring Count BCUT 2nd highest eigenvalue", desc_bcut_rc_hi2);
    REGISTER_BCUT_EXT("BCUTrc_hi3", "Ring Count BCUT 3rd highest eigenvalue", desc_bcut_rc_hi3);
    REGISTER_BCUT_EXT("BCUTrc_lo1", "Ring Count BCUT lowest eigenvalue", desc_bcut_rc_lo1);
    REGISTER_BCUT_EXT("BCUTrc_lo2", "Ring Count BCUT 2nd lowest eigenvalue", desc_bcut_rc_lo2);
    REGISTER_BCUT_EXT("BCUTrc_lo3", "Ring Count BCUT 3rd lowest eigenvalue", desc_bcut_rc_lo3);
}
