/**
 * @file vsa.c
 * @brief Van der Waals Surface Area (VSA) descriptors
 *
 * Partitions molecular surface area into bins based on atomic properties:
 * - SlogP_VSA: Surface area by LogP contribution bins
 * - SMR_VSA: Surface area by molar refractivity bins
 * - PEOE_VSA: Surface area by partial charge bins
 * - EState_VSA: Surface area by E-state bins
 *
 * Reference: Labute, J. Mol. Graph. Model. 2000
 *
 * All O(n) complexity.
 */

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Atomic Properties
 * ============================================================================ */

/* Van der Waals radii (Bondi) in Angstroms */
static double get_vdw_radius(element_t elem) {
    switch (elem) {
        case ELEM_H:  return 1.20;
        case ELEM_C:  return 1.70;
        case ELEM_N:  return 1.55;
        case ELEM_O:  return 1.52;
        case ELEM_F:  return 1.47;
        case ELEM_P:  return 1.80;
        case ELEM_S:  return 1.80;
        case ELEM_Cl: return 1.75;
        case ELEM_Br: return 1.85;
        case ELEM_I:  return 1.98;
        default:      return 1.70;
    }
}

/* Approximate VdW surface area contribution (simplified Bondi) */
static double get_vsa(element_t elem) {
    double r = get_vdw_radius(elem);
    return 4.0 * 3.14159 * r * r;  /* Sphere surface area */
}

/* Wildman-Crippen LogP atomic contributions */
static double get_logp_contrib(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Count neighbors */
    int heavy = 0, h_count = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
        else heavy++;
    }

    /* Simplified Wildman-Crippen contributions */
    switch (atom->element) {
        case ELEM_C:
            if (atom->aromatic) return 0.13;
            if (heavy == 1) return 0.29;  /* CH3 */
            if (heavy == 2) return 0.19;  /* CH2 */
            if (heavy == 3) return 0.09;  /* CH */
            return -0.03;  /* quaternary */

        case ELEM_N:
            if (atom->aromatic) return -0.55;
            if (h_count >= 2) return -1.03;  /* NH2 */
            if (h_count == 1) return -0.54;  /* NH */
            return -0.28;  /* tertiary */

        case ELEM_O:
            if (h_count >= 1) return -0.47;  /* OH */
            return -0.27;  /* ether */

        case ELEM_S:
            if (h_count >= 1) return 0.42;  /* SH */
            return 0.29;  /* thioether */

        case ELEM_F:  return 0.38;
        case ELEM_Cl: return 0.71;
        case ELEM_Br: return 0.87;
        case ELEM_I:  return 1.15;
        case ELEM_P:  return 0.16;
        case ELEM_H:  return 0.12;
        default:      return 0.0;
    }
}

/* Molar refractivity contributions */
static double get_mr_contrib(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    switch (atom->element) {
        case ELEM_C:  return 2.50;
        case ELEM_N:  return 2.74;
        case ELEM_O:  return 1.55;
        case ELEM_S:  return 7.69;
        case ELEM_F:  return 0.92;
        case ELEM_Cl: return 5.85;
        case ELEM_Br: return 8.88;
        case ELEM_I:  return 14.0;
        case ELEM_P:  return 6.92;
        case ELEM_H:  return 1.06;
        default:      return 2.0;
    }
}

/* ============================================================================
 * SlogP_VSA Descriptors (12 bins)
 * ============================================================================ */

/* LogP contribution bins */
static const double SLOGP_BINS[] = {
    -0.40, -0.20, 0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60
};
#define NUM_SLOGP_BINS 12

static void compute_slogp_vsa(const molecule_t* mol, double* values) {
    for (int i = 0; i < NUM_SLOGP_BINS; i++) values[i] = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;

        double logp = get_logp_contrib(mol, i);
        double vsa = get_vsa(mol->atoms[i].element);

        int bin = 0;
        for (int b = 0; b < NUM_SLOGP_BINS - 1; b++) {
            if (logp >= SLOGP_BINS[b]) bin = b + 1;
        }
        values[bin] += vsa;
    }
}

/* Individual SlogP_VSA descriptors */
#define SLOGP_VSA_FUNC(n) \
static cchem_status_t slogp_vsa##n(const molecule_t* mol, descriptor_value_t* value) { \
    double vals[NUM_SLOGP_BINS]; \
    compute_slogp_vsa(mol, vals); \
    value->d = vals[n-1]; \
    return CCHEM_OK; \
}

SLOGP_VSA_FUNC(1)
SLOGP_VSA_FUNC(2)
SLOGP_VSA_FUNC(3)
SLOGP_VSA_FUNC(4)
SLOGP_VSA_FUNC(5)
SLOGP_VSA_FUNC(6)
SLOGP_VSA_FUNC(7)
SLOGP_VSA_FUNC(8)
SLOGP_VSA_FUNC(9)
SLOGP_VSA_FUNC(10)
SLOGP_VSA_FUNC(11)
SLOGP_VSA_FUNC(12)

/* ============================================================================
 * SMR_VSA Descriptors (10 bins)
 * ============================================================================ */

static const double SMR_BINS[] = {
    1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 4.25, 5.00
};
#define NUM_SMR_BINS 10

static void compute_smr_vsa(const molecule_t* mol, double* values) {
    for (int i = 0; i < NUM_SMR_BINS; i++) values[i] = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;

        double mr = get_mr_contrib(mol, i);
        double vsa = get_vsa(mol->atoms[i].element);

        int bin = 0;
        for (int b = 0; b < NUM_SMR_BINS - 1; b++) {
            if (mr >= SMR_BINS[b]) bin = b + 1;
        }
        values[bin] += vsa;
    }
}

#define SMR_VSA_FUNC(n) \
static cchem_status_t smr_vsa##n(const molecule_t* mol, descriptor_value_t* value) { \
    double vals[NUM_SMR_BINS]; \
    compute_smr_vsa(mol, vals); \
    value->d = vals[n-1]; \
    return CCHEM_OK; \
}

SMR_VSA_FUNC(1)
SMR_VSA_FUNC(2)
SMR_VSA_FUNC(3)
SMR_VSA_FUNC(4)
SMR_VSA_FUNC(5)
SMR_VSA_FUNC(6)
SMR_VSA_FUNC(7)
SMR_VSA_FUNC(8)
SMR_VSA_FUNC(9)
SMR_VSA_FUNC(10)

/* ============================================================================
 * PEOE_VSA Descriptors (14 bins by partial charge)
 * ============================================================================ */

static const double PEOE_BINS[] = {
    -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30
};
#define NUM_PEOE_BINS 14

/* Simplified Gasteiger charge estimate */
static double get_charge_estimate(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Simple electronegativity-based estimate */
    switch (atom->element) {
        case ELEM_O:  return -0.25;
        case ELEM_N:  return atom->aromatic ? -0.10 : -0.15;
        case ELEM_S:  return -0.05;
        case ELEM_F:  return -0.20;
        case ELEM_Cl: return -0.10;
        case ELEM_Br: return -0.08;
        case ELEM_I:  return -0.05;
        case ELEM_C:
            /* Check for attached electronegative atoms */
            for (int i = 0; i < atom->num_neighbors; i++) {
                element_t ne = mol->atoms[atom->neighbors[i]].element;
                if (ne == ELEM_O || ne == ELEM_N || ne == ELEM_F) {
                    return 0.15;
                }
            }
            return 0.02;
        case ELEM_H:  return 0.08;
        default:      return 0.0;
    }
}

static void compute_peoe_vsa(const molecule_t* mol, double* values) {
    for (int i = 0; i < NUM_PEOE_BINS; i++) values[i] = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;

        double charge = get_charge_estimate(mol, i);
        double vsa = get_vsa(mol->atoms[i].element);

        int bin = 0;
        for (int b = 0; b < NUM_PEOE_BINS - 1; b++) {
            if (charge >= PEOE_BINS[b]) bin = b + 1;
        }
        values[bin] += vsa;
    }
}

#define PEOE_VSA_FUNC(n) \
static cchem_status_t peoe_vsa##n(const molecule_t* mol, descriptor_value_t* value) { \
    double vals[NUM_PEOE_BINS]; \
    compute_peoe_vsa(mol, vals); \
    value->d = vals[n-1]; \
    return CCHEM_OK; \
}

PEOE_VSA_FUNC(1)
PEOE_VSA_FUNC(2)
PEOE_VSA_FUNC(3)
PEOE_VSA_FUNC(4)
PEOE_VSA_FUNC(5)
PEOE_VSA_FUNC(6)
PEOE_VSA_FUNC(7)
PEOE_VSA_FUNC(8)
PEOE_VSA_FUNC(9)
PEOE_VSA_FUNC(10)
PEOE_VSA_FUNC(11)
PEOE_VSA_FUNC(12)
PEOE_VSA_FUNC(13)
PEOE_VSA_FUNC(14)

/* ============================================================================
 * EState_VSA Descriptors (11 bins)
 * ============================================================================ */

/* Intrinsic state value */
static double get_intrinsic_state(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    int heavy = 0, h_count = atom->implicit_h_count;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element == ELEM_H) h_count++;
        else heavy++;
    }

    /* Simplified intrinsic state: (2/N^2) * (delta_v + 1) / delta */
    int delta = heavy;
    if (delta == 0) delta = 1;

    double delta_v = 0.0;
    switch (atom->element) {
        case ELEM_C:  delta_v = 4 - h_count; break;
        case ELEM_N:  delta_v = 5 - h_count; break;
        case ELEM_O:  delta_v = 6 - h_count; break;
        case ELEM_S:  delta_v = 6 - h_count; break;
        case ELEM_F:  delta_v = 7; break;
        case ELEM_Cl: delta_v = 7; break;
        case ELEM_Br: delta_v = 7; break;
        case ELEM_I:  delta_v = 7; break;
        case ELEM_P:  delta_v = 5 - h_count; break;
        default:      delta_v = 4 - h_count; break;
    }

    int n = 2;  /* Principal quantum number estimate */
    if (atom->element >= ELEM_Na) n = 3;
    if (atom->element >= ELEM_K) n = 4;

    return (2.0 / (n * n)) * (delta_v + 1) / delta;
}

static const double ESTATE_BINS[] = {
    -0.39, 0.29, 0.72, 1.17, 1.54, 1.81, 2.05, 2.50, 3.05, 4.69
};
#define NUM_ESTATE_BINS 11

static void compute_estate_vsa(const molecule_t* mol, double* values) {
    for (int i = 0; i < NUM_ESTATE_BINS; i++) values[i] = 0.0;

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;

        double estate = get_intrinsic_state(mol, i);
        double vsa = get_vsa(mol->atoms[i].element);

        int bin = 0;
        for (int b = 0; b < NUM_ESTATE_BINS - 1; b++) {
            if (estate >= ESTATE_BINS[b]) bin = b + 1;
        }
        values[bin] += vsa;
    }
}

#define ESTATE_VSA_FUNC(n) \
static cchem_status_t estate_vsa##n(const molecule_t* mol, descriptor_value_t* value) { \
    double vals[NUM_ESTATE_BINS]; \
    compute_estate_vsa(mol, vals); \
    value->d = vals[n-1]; \
    return CCHEM_OK; \
}

ESTATE_VSA_FUNC(1)
ESTATE_VSA_FUNC(2)
ESTATE_VSA_FUNC(3)
ESTATE_VSA_FUNC(4)
ESTATE_VSA_FUNC(5)
ESTATE_VSA_FUNC(6)
ESTATE_VSA_FUNC(7)
ESTATE_VSA_FUNC(8)
ESTATE_VSA_FUNC(9)
ESTATE_VSA_FUNC(10)
ESTATE_VSA_FUNC(11)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_VSA(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_STERIC; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_vsa(void) {
    /* SlogP_VSA (12 bins) */
    REGISTER_VSA("SlogP_VSA1", "VSA LogP<-0.4", slogp_vsa1);
    REGISTER_VSA("SlogP_VSA2", "VSA LogP[-0.4,-0.2)", slogp_vsa2);
    REGISTER_VSA("SlogP_VSA3", "VSA LogP[-0.2,0.0)", slogp_vsa3);
    REGISTER_VSA("SlogP_VSA4", "VSA LogP[0.0,0.1)", slogp_vsa4);
    REGISTER_VSA("SlogP_VSA5", "VSA LogP[0.1,0.15)", slogp_vsa5);
    REGISTER_VSA("SlogP_VSA6", "VSA LogP[0.15,0.2)", slogp_vsa6);
    REGISTER_VSA("SlogP_VSA7", "VSA LogP[0.2,0.25)", slogp_vsa7);
    REGISTER_VSA("SlogP_VSA8", "VSA LogP[0.25,0.3)", slogp_vsa8);
    REGISTER_VSA("SlogP_VSA9", "VSA LogP[0.3,0.4)", slogp_vsa9);
    REGISTER_VSA("SlogP_VSA10", "VSA LogP[0.4,0.5)", slogp_vsa10);
    REGISTER_VSA("SlogP_VSA11", "VSA LogP[0.5,0.6)", slogp_vsa11);
    REGISTER_VSA("SlogP_VSA12", "VSA LogP>=0.6", slogp_vsa12);

    /* SMR_VSA (10 bins) */
    REGISTER_VSA("SMR_VSA1", "VSA MR<1.29", smr_vsa1);
    REGISTER_VSA("SMR_VSA2", "VSA MR[1.29,1.82)", smr_vsa2);
    REGISTER_VSA("SMR_VSA3", "VSA MR[1.82,2.24)", smr_vsa3);
    REGISTER_VSA("SMR_VSA4", "VSA MR[2.24,2.45)", smr_vsa4);
    REGISTER_VSA("SMR_VSA5", "VSA MR[2.45,2.75)", smr_vsa5);
    REGISTER_VSA("SMR_VSA6", "VSA MR[2.75,3.05)", smr_vsa6);
    REGISTER_VSA("SMR_VSA7", "VSA MR[3.05,3.63)", smr_vsa7);
    REGISTER_VSA("SMR_VSA8", "VSA MR[3.63,4.25)", smr_vsa8);
    REGISTER_VSA("SMR_VSA9", "VSA MR[4.25,5.0)", smr_vsa9);
    REGISTER_VSA("SMR_VSA10", "VSA MR>=5.0", smr_vsa10);

    /* PEOE_VSA (14 bins) */
    REGISTER_VSA("PEOE_VSA1", "VSA charge<-0.30", peoe_vsa1);
    REGISTER_VSA("PEOE_VSA2", "VSA charge[-0.30,-0.25)", peoe_vsa2);
    REGISTER_VSA("PEOE_VSA3", "VSA charge[-0.25,-0.20)", peoe_vsa3);
    REGISTER_VSA("PEOE_VSA4", "VSA charge[-0.20,-0.15)", peoe_vsa4);
    REGISTER_VSA("PEOE_VSA5", "VSA charge[-0.15,-0.10)", peoe_vsa5);
    REGISTER_VSA("PEOE_VSA6", "VSA charge[-0.10,-0.05)", peoe_vsa6);
    REGISTER_VSA("PEOE_VSA7", "VSA charge[-0.05,0.00)", peoe_vsa7);
    REGISTER_VSA("PEOE_VSA8", "VSA charge[0.00,0.05)", peoe_vsa8);
    REGISTER_VSA("PEOE_VSA9", "VSA charge[0.05,0.10)", peoe_vsa9);
    REGISTER_VSA("PEOE_VSA10", "VSA charge[0.10,0.15)", peoe_vsa10);
    REGISTER_VSA("PEOE_VSA11", "VSA charge[0.15,0.20)", peoe_vsa11);
    REGISTER_VSA("PEOE_VSA12", "VSA charge[0.20,0.25)", peoe_vsa12);
    REGISTER_VSA("PEOE_VSA13", "VSA charge[0.25,0.30)", peoe_vsa13);
    REGISTER_VSA("PEOE_VSA14", "VSA charge>=0.30", peoe_vsa14);

    /* EState_VSA (11 bins) */
    REGISTER_VSA("EState_VSA1", "VSA EState<-0.39", estate_vsa1);
    REGISTER_VSA("EState_VSA2", "VSA EState[-0.39,0.29)", estate_vsa2);
    REGISTER_VSA("EState_VSA3", "VSA EState[0.29,0.72)", estate_vsa3);
    REGISTER_VSA("EState_VSA4", "VSA EState[0.72,1.17)", estate_vsa4);
    REGISTER_VSA("EState_VSA5", "VSA EState[1.17,1.54)", estate_vsa5);
    REGISTER_VSA("EState_VSA6", "VSA EState[1.54,1.81)", estate_vsa6);
    REGISTER_VSA("EState_VSA7", "VSA EState[1.81,2.05)", estate_vsa7);
    REGISTER_VSA("EState_VSA8", "VSA EState[2.05,2.50)", estate_vsa8);
    REGISTER_VSA("EState_VSA9", "VSA EState[2.50,3.05)", estate_vsa9);
    REGISTER_VSA("EState_VSA10", "VSA EState[3.05,4.69)", estate_vsa10);
    REGISTER_VSA("EState_VSA11", "VSA EState>=4.69", estate_vsa11);
}
