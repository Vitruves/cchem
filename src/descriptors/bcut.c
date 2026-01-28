/**
 * @file bcut.c
 * @brief BCUT Descriptors - Burden-CAS-University of Texas eigenvalue descriptors
 *
 * BCUT descriptors are eigenvalues of the Burden matrix, a modified adjacency
 * matrix where:
 * - Diagonal: atomic property values (mass, charge, EN, polarizability, etc.)
 * - Off-diagonal: 0.1 * sqrt(conventional bond order) for bonded atoms, 0.001 otherwise
 *
 * For each property, we compute eigenvalues and take:
 * - Highest eigenvalue (BCUT_X_hi)
 * - Lowest eigenvalue (BCUT_X_lo)
 * - Plus additional high/low pairs for larger molecules
 *
 * Properties: mass (m), charge (c), electronegativity (e), polarizability (p),
 *             ionization potential (i), hydrogen bond donors (d), acceptors (a), logP (l)
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

/* Number of properties */
#define NUM_BCUT_PROPERTIES 8

/* Total BCUT descriptors */
#define NUM_BCUT_DESCRIPTORS (NUM_BCUT_PROPERTIES * NUM_EIGENVALUES)

/* Maximum atoms for stack allocation */
#define MAX_ATOMS_STACK 128

/* ============================================================================
 * Atomic Property Tables
 * ============================================================================ */

/* Atomic masses (Da) */
static const double BCUT_MASS[] = {
    [ELEM_H]  = 1.008,   [ELEM_C]  = 12.011,  [ELEM_N]  = 14.007,
    [ELEM_O]  = 15.999,  [ELEM_F]  = 18.998,  [ELEM_P]  = 30.974,
    [ELEM_S]  = 32.065,  [ELEM_Cl] = 35.453,  [ELEM_Br] = 79.904,
    [ELEM_I]  = 126.904, [ELEM_Si] = 28.086,  [ELEM_B]  = 10.811,
    [ELEM_Se] = 78.96,   [ELEM_As] = 74.922,
};

/* Sanderson electronegativity */
static const double BCUT_EN[] = {
    [ELEM_H]  = 2.59,  [ELEM_C]  = 2.75,  [ELEM_N]  = 3.19,
    [ELEM_O]  = 3.65,  [ELEM_F]  = 4.00,  [ELEM_P]  = 2.52,
    [ELEM_S]  = 2.96,  [ELEM_Cl] = 3.48,  [ELEM_Br] = 3.22,
    [ELEM_I]  = 2.78,  [ELEM_Si] = 2.14,  [ELEM_B]  = 2.28,
    [ELEM_Se] = 2.76,  [ELEM_As] = 2.26,
};

/* Atomic polarizability (Å³) */
static const double BCUT_POL[] = {
    [ELEM_H]  = 0.387, [ELEM_C]  = 1.76,  [ELEM_N]  = 1.10,
    [ELEM_O]  = 0.802, [ELEM_F]  = 0.557, [ELEM_P]  = 3.63,
    [ELEM_S]  = 2.90,  [ELEM_Cl] = 2.18,  [ELEM_Br] = 3.05,
    [ELEM_I]  = 4.7,   [ELEM_Si] = 5.38,  [ELEM_B]  = 3.03,
    [ELEM_Se] = 3.77,  [ELEM_As] = 4.31,
};

/* First ionization potential (eV) */
static const double BCUT_IP[] = {
    [ELEM_H]  = 13.598, [ELEM_C]  = 11.260, [ELEM_N]  = 14.534,
    [ELEM_O]  = 13.618, [ELEM_F]  = 17.422, [ELEM_P]  = 10.486,
    [ELEM_S]  = 10.360, [ELEM_Cl] = 12.967, [ELEM_Br] = 11.814,
    [ELEM_I]  = 10.451, [ELEM_Si] = 8.151,  [ELEM_B]  = 8.298,
    [ELEM_Se] = 9.752,  [ELEM_As] = 9.815,
};

/* Wildman-Crippen LogP contributions */
static const double BCUT_LOGP[] = {
    [ELEM_H]  = 0.123,  [ELEM_C]  = 0.291,  [ELEM_N]  = -0.534,
    [ELEM_O]  = -0.467, [ELEM_F]  = 0.375,  [ELEM_P]  = 0.8,
    [ELEM_S]  = 0.643,  [ELEM_Cl] = 0.713,  [ELEM_Br] = 1.013,
    [ELEM_I]  = 1.293,  [ELEM_Si] = 0.0,    [ELEM_B]  = 0.0,
    [ELEM_Se] = 0.7,    [ELEM_As] = 0.5,
};

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static inline double get_property(element_t elem, const double* table, double default_val) {
    if (elem <= 0 || elem >= 128) return default_val;
    double v = table[elem];
    return (v != 0.0) ? v : default_val;
}

/* Check if atom is H-bond donor (has H attached to N, O, S) */
static int is_hbd(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    element_t elem = atom->element;
    if (elem != ELEM_N && elem != ELEM_O && elem != ELEM_S) return 0;

    /* Count attached hydrogens */
    int h_count = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
    }
    return h_count > 0 ? 1 : 0;
}

/* Check if atom is H-bond acceptor (N, O with lone pairs) */
static int is_hba(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    element_t elem = atom->element;
    if (elem == ELEM_O) return 1;
    if (elem == ELEM_N) {
        /* Exclude positively charged N and N in nitro */
        if (atom->charge > 0) return 0;
        return 1;
    }
    if (elem == ELEM_S) return 1;
    if (elem == ELEM_F) return 1;
    return 0;
}

/* ============================================================================
 * Jacobi Eigenvalue Algorithm for Symmetric Matrices
 * ============================================================================ */

#define JACOBI_MAX_ITER 50
#define JACOBI_EPS 1e-10

/**
 * Compute eigenvalues of a symmetric matrix using Jacobi method.
 * Matrix is destroyed in process. Eigenvalues returned sorted descending.
 */
static int jacobi_eigenvalues(double* A, int n, double* eigenvalues) {
    if (n <= 0) return 0;
    if (n == 1) {
        eigenvalues[0] = A[0];
        return 1;
    }

    /* Work arrays */
    double* b = (double*)alloca(n * sizeof(double));
    double* z = (double*)alloca(n * sizeof(double));

    /* Initialize */
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = A[i * n + i];  /* Diagonal */
        b[i] = eigenvalues[i];
        z[i] = 0.0;
    }

    for (int iter = 0; iter < JACOBI_MAX_ITER; iter++) {
        /* Sum off-diagonal elements */
        double sm = 0.0;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                sm += fabs(A[i * n + j]);
            }
        }

        if (sm < JACOBI_EPS) break;  /* Converged */

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

                    /* Rotations */
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
 * BCUT Computation
 * ============================================================================ */

/**
 * Build Burden matrix for given property.
 * Diagonal: property values
 * Off-diagonal: 0.1 * sqrt(bond_order) for bonded, 0.001 for non-bonded
 */
static void build_burden_matrix(const molecule_t* mol,
                                 const int* heavy_idx, int n_heavy,
                                 const double* props,
                                 double* matrix) {
    /* Initialize with small off-diagonal values */
    for (int i = 0; i < n_heavy; i++) {
        for (int j = 0; j < n_heavy; j++) {
            if (i == j) {
                matrix[i * n_heavy + j] = props[i];
            } else {
                matrix[i * n_heavy + j] = 0.001;
            }
        }
    }

    /* Create reverse mapping */
    int* reverse_map = (int*)alloca(mol->num_atoms * sizeof(int));
    memset(reverse_map, -1, mol->num_atoms * sizeof(int));
    for (int i = 0; i < n_heavy; i++) {
        reverse_map[heavy_idx[i]] = i;
    }

    /* Set bond entries */
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

/**
 * Compute all BCUT descriptors for a molecule.
 * Returns 48 values: 8 properties × 6 eigenvalues (3 high + 3 low)
 */
int descriptors_compute_bcut_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    /* Initialize all to 0 */
    for (int i = 0; i < NUM_BCUT_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    /* Count and index heavy atoms */
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    }

    if (n_heavy == 0) return NUM_BCUT_DESCRIPTORS;
    if (n_heavy == 1) {
        /* Single atom: eigenvalue = property value, no variation */
        return NUM_BCUT_DESCRIPTORS;
    }

    /* Allocate arrays */
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

    /* Build heavy atom index */
    int idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_idx[idx++] = i;
        }
    }

    int out_idx = 0;

    /* Property 0: Mass */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_MASS, 12.0);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    /* Store 3 highest and 3 lowest */
    values[out_idx++].d = eigenvalues[0];  /* hi1 */
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;  /* hi2 */
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;  /* hi3 */
    values[out_idx++].d = eigenvalues[n_heavy - 1];  /* lo1 */
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;  /* lo2 */
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;  /* lo3 */

    /* Property 1: Gasteiger charge (simplified - use formal charge) */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = mol->atoms[heavy_idx[i]].charge;
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    values[out_idx++].d = eigenvalues[0];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    values[out_idx++].d = eigenvalues[n_heavy - 1];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;

    /* Property 2: Electronegativity */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_EN, 2.75);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    values[out_idx++].d = eigenvalues[0];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    values[out_idx++].d = eigenvalues[n_heavy - 1];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;

    /* Property 3: Polarizability */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_POL, 1.76);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    values[out_idx++].d = eigenvalues[0];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    values[out_idx++].d = eigenvalues[n_heavy - 1];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;

    /* Property 4: Ionization potential */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_IP, 11.26);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    values[out_idx++].d = eigenvalues[0];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    values[out_idx++].d = eigenvalues[n_heavy - 1];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;

    /* Property 5: H-bond donors */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = is_hbd(mol, heavy_idx[i]) ? 1.0 : 0.0;
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    values[out_idx++].d = eigenvalues[0];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    values[out_idx++].d = eigenvalues[n_heavy - 1];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;

    /* Property 6: H-bond acceptors */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = is_hba(mol, heavy_idx[i]) ? 1.0 : 0.0;
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    values[out_idx++].d = eigenvalues[0];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    values[out_idx++].d = eigenvalues[n_heavy - 1];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;

    /* Property 7: LogP contribution */
    for (int i = 0; i < n_heavy; i++) {
        props[i] = get_property(mol->atoms[heavy_idx[i]].element, BCUT_LOGP, 0.0);
    }
    build_burden_matrix(mol, heavy_idx, n_heavy, props, matrix);
    jacobi_eigenvalues(matrix, n_heavy, eigenvalues);
    values[out_idx++].d = eigenvalues[0];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[1] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[2] : 0.0;
    values[out_idx++].d = eigenvalues[n_heavy - 1];
    values[out_idx++].d = (n_heavy > 1) ? eigenvalues[n_heavy - 2] : 0.0;
    values[out_idx++].d = (n_heavy > 2) ? eigenvalues[n_heavy - 3] : 0.0;

    if (heap_alloc) {
        free(heavy_idx); free(props); free(matrix); free(eigenvalues);
    }

    return NUM_BCUT_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* bcut_cached_mol = NULL;
static _Thread_local descriptor_value_t bcut_cached_values[NUM_BCUT_DESCRIPTORS];

static inline void ensure_bcut_computed(const molecule_t* mol) {
    if (bcut_cached_mol != mol) {
        descriptors_compute_bcut_all(mol, bcut_cached_values);
        bcut_cached_mol = mol;
    }
}

#define DEFINE_BCUT_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_bcut_computed(mol); \
    value->d = bcut_cached_values[idx].d; \
    return CCHEM_OK; \
}

/* Mass-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_m_hi1, 0)
DEFINE_BCUT_FUNC(bcut_m_hi2, 1)
DEFINE_BCUT_FUNC(bcut_m_hi3, 2)
DEFINE_BCUT_FUNC(bcut_m_lo1, 3)
DEFINE_BCUT_FUNC(bcut_m_lo2, 4)
DEFINE_BCUT_FUNC(bcut_m_lo3, 5)

/* Charge-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_c_hi1, 6)
DEFINE_BCUT_FUNC(bcut_c_hi2, 7)
DEFINE_BCUT_FUNC(bcut_c_hi3, 8)
DEFINE_BCUT_FUNC(bcut_c_lo1, 9)
DEFINE_BCUT_FUNC(bcut_c_lo2, 10)
DEFINE_BCUT_FUNC(bcut_c_lo3, 11)

/* Electronegativity-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_e_hi1, 12)
DEFINE_BCUT_FUNC(bcut_e_hi2, 13)
DEFINE_BCUT_FUNC(bcut_e_hi3, 14)
DEFINE_BCUT_FUNC(bcut_e_lo1, 15)
DEFINE_BCUT_FUNC(bcut_e_lo2, 16)
DEFINE_BCUT_FUNC(bcut_e_lo3, 17)

/* Polarizability-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_p_hi1, 18)
DEFINE_BCUT_FUNC(bcut_p_hi2, 19)
DEFINE_BCUT_FUNC(bcut_p_hi3, 20)
DEFINE_BCUT_FUNC(bcut_p_lo1, 21)
DEFINE_BCUT_FUNC(bcut_p_lo2, 22)
DEFINE_BCUT_FUNC(bcut_p_lo3, 23)

/* Ionization potential-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_i_hi1, 24)
DEFINE_BCUT_FUNC(bcut_i_hi2, 25)
DEFINE_BCUT_FUNC(bcut_i_hi3, 26)
DEFINE_BCUT_FUNC(bcut_i_lo1, 27)
DEFINE_BCUT_FUNC(bcut_i_lo2, 28)
DEFINE_BCUT_FUNC(bcut_i_lo3, 29)

/* H-bond donor-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_d_hi1, 30)
DEFINE_BCUT_FUNC(bcut_d_hi2, 31)
DEFINE_BCUT_FUNC(bcut_d_hi3, 32)
DEFINE_BCUT_FUNC(bcut_d_lo1, 33)
DEFINE_BCUT_FUNC(bcut_d_lo2, 34)
DEFINE_BCUT_FUNC(bcut_d_lo3, 35)

/* H-bond acceptor-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_a_hi1, 36)
DEFINE_BCUT_FUNC(bcut_a_hi2, 37)
DEFINE_BCUT_FUNC(bcut_a_hi3, 38)
DEFINE_BCUT_FUNC(bcut_a_lo1, 39)
DEFINE_BCUT_FUNC(bcut_a_lo2, 40)
DEFINE_BCUT_FUNC(bcut_a_lo3, 41)

/* LogP-weighted BCUT */
DEFINE_BCUT_FUNC(bcut_l_hi1, 42)
DEFINE_BCUT_FUNC(bcut_l_hi2, 43)
DEFINE_BCUT_FUNC(bcut_l_hi3, 44)
DEFINE_BCUT_FUNC(bcut_l_lo1, 45)
DEFINE_BCUT_FUNC(bcut_l_lo2, 46)
DEFINE_BCUT_FUNC(bcut_l_lo3, 47)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_BCUT(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_bcut(void) {
    /* Mass-weighted BCUT */
    REGISTER_BCUT("BCUTm_hi1", "Mass BCUT highest eigenvalue", desc_bcut_m_hi1);
    REGISTER_BCUT("BCUTm_hi2", "Mass BCUT 2nd highest eigenvalue", desc_bcut_m_hi2);
    REGISTER_BCUT("BCUTm_hi3", "Mass BCUT 3rd highest eigenvalue", desc_bcut_m_hi3);
    REGISTER_BCUT("BCUTm_lo1", "Mass BCUT lowest eigenvalue", desc_bcut_m_lo1);
    REGISTER_BCUT("BCUTm_lo2", "Mass BCUT 2nd lowest eigenvalue", desc_bcut_m_lo2);
    REGISTER_BCUT("BCUTm_lo3", "Mass BCUT 3rd lowest eigenvalue", desc_bcut_m_lo3);

    /* Charge-weighted BCUT */
    REGISTER_BCUT("BCUTc_hi1", "Charge BCUT highest eigenvalue", desc_bcut_c_hi1);
    REGISTER_BCUT("BCUTc_hi2", "Charge BCUT 2nd highest eigenvalue", desc_bcut_c_hi2);
    REGISTER_BCUT("BCUTc_hi3", "Charge BCUT 3rd highest eigenvalue", desc_bcut_c_hi3);
    REGISTER_BCUT("BCUTc_lo1", "Charge BCUT lowest eigenvalue", desc_bcut_c_lo1);
    REGISTER_BCUT("BCUTc_lo2", "Charge BCUT 2nd lowest eigenvalue", desc_bcut_c_lo2);
    REGISTER_BCUT("BCUTc_lo3", "Charge BCUT 3rd lowest eigenvalue", desc_bcut_c_lo3);

    /* Electronegativity-weighted BCUT */
    REGISTER_BCUT("BCUTe_hi1", "EN BCUT highest eigenvalue", desc_bcut_e_hi1);
    REGISTER_BCUT("BCUTe_hi2", "EN BCUT 2nd highest eigenvalue", desc_bcut_e_hi2);
    REGISTER_BCUT("BCUTe_hi3", "EN BCUT 3rd highest eigenvalue", desc_bcut_e_hi3);
    REGISTER_BCUT("BCUTe_lo1", "EN BCUT lowest eigenvalue", desc_bcut_e_lo1);
    REGISTER_BCUT("BCUTe_lo2", "EN BCUT 2nd lowest eigenvalue", desc_bcut_e_lo2);
    REGISTER_BCUT("BCUTe_lo3", "EN BCUT 3rd lowest eigenvalue", desc_bcut_e_lo3);

    /* Polarizability-weighted BCUT */
    REGISTER_BCUT("BCUTp_hi1", "Polarizability BCUT highest eigenvalue", desc_bcut_p_hi1);
    REGISTER_BCUT("BCUTp_hi2", "Polarizability BCUT 2nd highest eigenvalue", desc_bcut_p_hi2);
    REGISTER_BCUT("BCUTp_hi3", "Polarizability BCUT 3rd highest eigenvalue", desc_bcut_p_hi3);
    REGISTER_BCUT("BCUTp_lo1", "Polarizability BCUT lowest eigenvalue", desc_bcut_p_lo1);
    REGISTER_BCUT("BCUTp_lo2", "Polarizability BCUT 2nd lowest eigenvalue", desc_bcut_p_lo2);
    REGISTER_BCUT("BCUTp_lo3", "Polarizability BCUT 3rd lowest eigenvalue", desc_bcut_p_lo3);

    /* Ionization potential-weighted BCUT */
    REGISTER_BCUT("BCUTi_hi1", "IP BCUT highest eigenvalue", desc_bcut_i_hi1);
    REGISTER_BCUT("BCUTi_hi2", "IP BCUT 2nd highest eigenvalue", desc_bcut_i_hi2);
    REGISTER_BCUT("BCUTi_hi3", "IP BCUT 3rd highest eigenvalue", desc_bcut_i_hi3);
    REGISTER_BCUT("BCUTi_lo1", "IP BCUT lowest eigenvalue", desc_bcut_i_lo1);
    REGISTER_BCUT("BCUTi_lo2", "IP BCUT 2nd lowest eigenvalue", desc_bcut_i_lo2);
    REGISTER_BCUT("BCUTi_lo3", "IP BCUT 3rd lowest eigenvalue", desc_bcut_i_lo3);

    /* H-bond donor-weighted BCUT */
    REGISTER_BCUT("BCUTd_hi1", "HBD BCUT highest eigenvalue", desc_bcut_d_hi1);
    REGISTER_BCUT("BCUTd_hi2", "HBD BCUT 2nd highest eigenvalue", desc_bcut_d_hi2);
    REGISTER_BCUT("BCUTd_hi3", "HBD BCUT 3rd highest eigenvalue", desc_bcut_d_hi3);
    REGISTER_BCUT("BCUTd_lo1", "HBD BCUT lowest eigenvalue", desc_bcut_d_lo1);
    REGISTER_BCUT("BCUTd_lo2", "HBD BCUT 2nd lowest eigenvalue", desc_bcut_d_lo2);
    REGISTER_BCUT("BCUTd_lo3", "HBD BCUT 3rd lowest eigenvalue", desc_bcut_d_lo3);

    /* H-bond acceptor-weighted BCUT */
    REGISTER_BCUT("BCUTa_hi1", "HBA BCUT highest eigenvalue", desc_bcut_a_hi1);
    REGISTER_BCUT("BCUTa_hi2", "HBA BCUT 2nd highest eigenvalue", desc_bcut_a_hi2);
    REGISTER_BCUT("BCUTa_hi3", "HBA BCUT 3rd highest eigenvalue", desc_bcut_a_hi3);
    REGISTER_BCUT("BCUTa_lo1", "HBA BCUT lowest eigenvalue", desc_bcut_a_lo1);
    REGISTER_BCUT("BCUTa_lo2", "HBA BCUT 2nd lowest eigenvalue", desc_bcut_a_lo2);
    REGISTER_BCUT("BCUTa_lo3", "HBA BCUT 3rd lowest eigenvalue", desc_bcut_a_lo3);

    /* LogP-weighted BCUT */
    REGISTER_BCUT("BCUTl_hi1", "LogP BCUT highest eigenvalue", desc_bcut_l_hi1);
    REGISTER_BCUT("BCUTl_hi2", "LogP BCUT 2nd highest eigenvalue", desc_bcut_l_hi2);
    REGISTER_BCUT("BCUTl_hi3", "LogP BCUT 3rd highest eigenvalue", desc_bcut_l_hi3);
    REGISTER_BCUT("BCUTl_lo1", "LogP BCUT lowest eigenvalue", desc_bcut_l_lo1);
    REGISTER_BCUT("BCUTl_lo2", "LogP BCUT 2nd lowest eigenvalue", desc_bcut_l_lo2);
    REGISTER_BCUT("BCUTl_lo3", "LogP BCUT 3rd lowest eigenvalue", desc_bcut_l_lo3);
}
