/**
 * @file infocontent.c
 * @brief Information Content and Molecular Complexity Descriptors
 *
 * Based on Shannon entropy theory applied to molecular structure:
 * - IC: Information content based on atom equivalence classes
 * - SIC: Structural information content
 * - CIC: Complementary information content
 * - BIC: Bonding information content
 * - Neighborhood complexity measures
 *
 * Total: 24 descriptors
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

#define NUM_IC_DESCRIPTORS 24
#define MAX_IC_ATOMS 512

/* ============================================================================
 * Shannon Entropy Calculation
 * H = -Σ(p_i * log2(p_i))
 * ============================================================================ */

static double shannon_entropy(const int* counts, int n_classes, int total) {
    if (total <= 0) return 0.0;

    double h = 0.0;
    for (int i = 0; i < n_classes; i++) {
        if (counts[i] > 0) {
            double p = (double)counts[i] / total;
            h -= p * log2(p);
        }
    }
    return h;
}

/* ============================================================================
 * Atom Equivalence Classes (based on Morgan algorithm iteration)
 * ============================================================================ */

/**
 * Compute atom equivalence classes using extended connectivity
 * Returns number of unique classes
 */
static int compute_atom_classes(const molecule_t* mol, int* classes, int order) {
    int n = mol->num_atoms;
    if (n == 0 || n > MAX_IC_ATOMS) return 0;

    /* Initial class: element type */
    for (int i = 0; i < n; i++) {
        classes[i] = mol->atoms[i].element;
    }

    /* Iteratively refine classes based on neighbor sums */
    int* new_classes = (int*)alloca(n * sizeof(int));

    for (int iter = 0; iter < order; iter++) {
        for (int i = 0; i < n; i++) {
            const atom_t* atom = &mol->atoms[i];
            int sum = classes[i] * 1000;  /* Self contribution */
            for (int j = 0; j < atom->num_neighbors; j++) {
                sum += classes[atom->neighbors[j]];
            }
            new_classes[i] = sum;
        }
        memcpy(classes, new_classes, n * sizeof(int));
    }

    /* Count unique classes */
    int* sorted = (int*)alloca(n * sizeof(int));
    memcpy(sorted, classes, n * sizeof(int));

    /* Simple bubble sort for small n */
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (sorted[j] < sorted[i]) {
                int tmp = sorted[i];
                sorted[i] = sorted[j];
                sorted[j] = tmp;
            }
        }
    }

    int unique = 1;
    for (int i = 1; i < n; i++) {
        if (sorted[i] != sorted[i-1]) unique++;
    }

    return unique;
}

/**
 * Count atoms per equivalence class
 */
static int count_class_distribution(const int* classes, int n, int* counts, int max_classes) {
    /* Map classes to sequential indices */
    int* class_map = (int*)alloca(n * sizeof(int));
    int* unique_classes = (int*)alloca(n * sizeof(int));
    int n_unique = 0;

    for (int i = 0; i < n; i++) {
        int found = -1;
        for (int j = 0; j < n_unique; j++) {
            if (unique_classes[j] == classes[i]) {
                found = j;
                break;
            }
        }
        if (found < 0) {
            if (n_unique < max_classes) {
                unique_classes[n_unique] = classes[i];
                found = n_unique++;
            }
        }
        class_map[i] = found;
    }

    memset(counts, 0, max_classes * sizeof(int));
    for (int i = 0; i < n; i++) {
        if (class_map[i] >= 0 && class_map[i] < max_classes) {
            counts[class_map[i]]++;
        }
    }

    return n_unique;
}

/* ============================================================================
 * Helper Functions from Existing Code
 * ============================================================================ */

static int get_heavy_degree(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    int degree = 0;
    for (int i = 0; i < atom->num_neighbors; i++) {
        if (mol->atoms[atom->neighbors[i]].element != ELEM_H) {
            degree++;
        }
    }
    return degree;
}

static int count_heavy_atoms(const molecule_t* mol) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) count++;
    }
    return count;
}

/* ============================================================================
 * Batch Computation
 * ============================================================================ */

int descriptors_compute_ic_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_IC_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    int n = mol->num_atoms;
    int n_heavy = count_heavy_atoms(mol);

    if (n == 0 || n > MAX_IC_ATOMS) return NUM_IC_DESCRIPTORS;

    int* classes = (int*)alloca(n * sizeof(int));
    int counts[128];

    /* ===== IC0: Information content order 0 (element types) ===== */
    int elem_counts[128] = {0};
    for (int i = 0; i < n; i++) {
        int e = mol->atoms[i].element;
        if (e > 0 && e < 128) elem_counts[e]++;
    }
    double IC0 = shannon_entropy(elem_counts, 128, n);

    /* ===== IC1-IC5: Information content orders 1-5 ===== */
    double IC[6];
    IC[0] = IC0;
    for (int order = 1; order <= 5; order++) {
        (void)compute_atom_classes(mol, classes, order);
        int n_unique = count_class_distribution(classes, n, counts, 128);
        IC[order] = shannon_entropy(counts, n_unique, n);
    }

    /* ===== SIC: Structural Information Content ===== */
    /* SIC = IC / log2(n) - normalized by maximum possible entropy */
    double log2n = (n > 1) ? log2((double)n) : 1.0;
    double SIC0 = IC[0] / log2n;
    double SIC1 = IC[1] / log2n;
    double SIC2 = IC[2] / log2n;

    /* ===== CIC: Complementary Information Content ===== */
    /* CIC = log2(n) - IC */
    double CIC0 = log2n - IC[0];
    double CIC1 = log2n - IC[1];
    double CIC2 = log2n - IC[2];

    /* ===== TIC: Total Information Content ===== */
    /* TIC = n * IC */
    double TIC0 = n * IC[0];
    double TIC1 = n * IC[1];
    double TIC2 = n * IC[2];

    /* ===== BIC: Bonding Information Content ===== */
    /* Based on bond type distribution */
    int bond_counts[8] = {0};  /* single, double, triple, aromatic, etc. */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int idx = 0;
        if (bond->aromatic) idx = 3;
        else if (bond->type == BOND_TRIPLE) idx = 2;
        else if (bond->type == BOND_DOUBLE) idx = 1;
        else idx = 0;
        bond_counts[idx]++;
    }
    double BIC = shannon_entropy(bond_counts, 4, mol->num_bonds);

    /* ===== Neighborhood Complexity ===== */
    /* Based on degree distribution */
    int degree_counts[16] = {0};
    for (int i = 0; i < n; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        int d = get_heavy_degree(mol, i);
        if (d >= 0 && d < 16) degree_counts[d]++;
    }
    double NC = shannon_entropy(degree_counts, 16, n_heavy);

    /* ===== Ring Information Content ===== */
    /* Based on ring size distribution */
    int ring_counts[16] = {0};
    for (int r = 0; r < mol->num_rings && r < 64; r++) {
        int size = mol->rings[r].size;
        if (size >= 3 && size < 16) ring_counts[size]++;
    }
    double RIC = shannon_entropy(ring_counts, 16, mol->num_rings);

    /* ===== Heteroatom Information Content ===== */
    /* Based on heteroatom type distribution */
    int hetero_counts[16] = {0};
    int n_hetero = 0;
    for (int i = 0; i < n; i++) {
        element_t e = mol->atoms[i].element;
        if (e != ELEM_C && e != ELEM_H) {
            int idx = 0;
            switch (e) {
                case ELEM_N: idx = 0; break;
                case ELEM_O: idx = 1; break;
                case ELEM_S: idx = 2; break;
                case ELEM_P: idx = 3; break;
                case ELEM_F: idx = 4; break;
                case ELEM_Cl: idx = 5; break;
                case ELEM_Br: idx = 6; break;
                case ELEM_I: idx = 7; break;
                default: idx = 8; break;
            }
            hetero_counts[idx]++;
            n_hetero++;
        }
    }
    double HIC = shannon_entropy(hetero_counts, 9, n_hetero);

    /* ===== Molecular Complexity Index ===== */
    /* Bertz complexity: based on symmetry and bonding */
    double Bertz = 2 * n_heavy * IC[1] + mol->num_bonds * BIC;

    /* Store results */
    int idx = 0;
    values[idx++].d = IC[0];    /* 0: IC0 */
    values[idx++].d = IC[1];    /* 1: IC1 */
    values[idx++].d = IC[2];    /* 2: IC2 */
    values[idx++].d = IC[3];    /* 3: IC3 */
    values[idx++].d = IC[4];    /* 4: IC4 */
    values[idx++].d = IC[5];    /* 5: IC5 */
    values[idx++].d = SIC0;     /* 6: SIC0 */
    values[idx++].d = SIC1;     /* 7: SIC1 */
    values[idx++].d = SIC2;     /* 8: SIC2 */
    values[idx++].d = CIC0;     /* 9: CIC0 */
    values[idx++].d = CIC1;     /* 10: CIC1 */
    values[idx++].d = CIC2;     /* 11: CIC2 */
    values[idx++].d = TIC0;     /* 12: TIC0 */
    values[idx++].d = TIC1;     /* 13: TIC1 */
    values[idx++].d = TIC2;     /* 14: TIC2 */
    values[idx++].d = BIC;      /* 15: BIC */
    values[idx++].d = NC;       /* 16: Neighborhood Complexity */
    values[idx++].d = RIC;      /* 17: Ring Info Content */
    values[idx++].d = HIC;      /* 18: Heteroatom Info Content */
    values[idx++].d = Bertz;    /* 19: Bertz Complexity */
    values[idx++].d = (double)n_heavy;              /* 20: n for reference */
    values[idx++].d = log2n;                        /* 21: log2(n) */
    values[idx++].d = (n_heavy > 0) ? IC[1] * n_heavy : 0.0;  /* 22: Scaled IC1 */
    values[idx++].d = (n_heavy > 0) ? NC * n_heavy : 0.0;     /* 23: Scaled NC */

    return NUM_IC_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* ic_cached_mol = NULL;
static _Thread_local descriptor_value_t ic_cached_values[NUM_IC_DESCRIPTORS];

static inline void ensure_ic_computed(const molecule_t* mol) {
    if (ic_cached_mol != mol) {
        descriptors_compute_ic_all(mol, ic_cached_values);
        ic_cached_mol = mol;
    }
}

#define DEFINE_IC_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_ic_computed(mol); \
    value->d = ic_cached_values[idx].d; \
    return CCHEM_OK; \
}

DEFINE_IC_FUNC(ic0, 0)
DEFINE_IC_FUNC(ic1, 1)
DEFINE_IC_FUNC(ic2, 2)
DEFINE_IC_FUNC(ic3, 3)
DEFINE_IC_FUNC(ic4, 4)
DEFINE_IC_FUNC(ic5, 5)
DEFINE_IC_FUNC(sic0, 6)
DEFINE_IC_FUNC(sic1, 7)
DEFINE_IC_FUNC(sic2, 8)
DEFINE_IC_FUNC(cic0, 9)
DEFINE_IC_FUNC(cic1, 10)
DEFINE_IC_FUNC(cic2, 11)
DEFINE_IC_FUNC(tic0, 12)
DEFINE_IC_FUNC(tic1, 13)
DEFINE_IC_FUNC(tic2, 14)
DEFINE_IC_FUNC(bic, 15)
DEFINE_IC_FUNC(nc, 16)
DEFINE_IC_FUNC(ric, 17)
DEFINE_IC_FUNC(hic, 18)
DEFINE_IC_FUNC(bertz, 19)
DEFINE_IC_FUNC(ic_n, 20)
DEFINE_IC_FUNC(ic_log2n, 21)
DEFINE_IC_FUNC(ic1_scaled, 22)
DEFINE_IC_FUNC(nc_scaled, 23)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_IC(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_infocontent(void) {
    /* Information Content orders 0-5 */
    REGISTER_IC("IC0", "Information content order 0", desc_ic0);
    REGISTER_IC("IC1", "Information content order 1", desc_ic1);
    REGISTER_IC("IC2", "Information content order 2", desc_ic2);
    REGISTER_IC("IC3", "Information content order 3", desc_ic3);
    REGISTER_IC("IC4", "Information content order 4", desc_ic4);
    REGISTER_IC("IC5", "Information content order 5", desc_ic5);

    /* Structural Information Content */
    REGISTER_IC("SIC0", "Structural IC order 0", desc_sic0);
    REGISTER_IC("SIC1", "Structural IC order 1", desc_sic1);
    REGISTER_IC("SIC2", "Structural IC order 2", desc_sic2);

    /* Complementary Information Content */
    REGISTER_IC("CIC0", "Complementary IC order 0", desc_cic0);
    REGISTER_IC("CIC1", "Complementary IC order 1", desc_cic1);
    REGISTER_IC("CIC2", "Complementary IC order 2", desc_cic2);

    /* Total Information Content */
    REGISTER_IC("TIC0", "Total IC order 0", desc_tic0);
    REGISTER_IC("TIC1", "Total IC order 1", desc_tic1);
    REGISTER_IC("TIC2", "Total IC order 2", desc_tic2);

    /* Bonding and complexity */
    REGISTER_IC("BIC", "Bonding information content", desc_bic);
    REGISTER_IC("NeighborComplexity", "Neighborhood complexity", desc_nc);
    REGISTER_IC("RingIC", "Ring information content", desc_ric);
    REGISTER_IC("HeteroIC", "Heteroatom information content", desc_hic);
    REGISTER_IC("BertzCT", "Bertz complexity index", desc_bertz);

    /* Reference values */
    REGISTER_IC("IC_N", "Number of atoms for IC", desc_ic_n);
    REGISTER_IC("IC_Log2N", "Log2(n) for IC normalization", desc_ic_log2n);
    REGISTER_IC("IC1_Scaled", "Scaled IC1 (IC1 * n)", desc_ic1_scaled);
    REGISTER_IC("NC_Scaled", "Scaled neighborhood complexity", desc_nc_scaled);
}
