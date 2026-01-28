/**
 * @file ringcomplexity.c
 * @brief Ring System Complexity Descriptors
 *
 * Descriptors characterizing ring systems:
 * - Ring counts by size and type
 * - Fusion and spiro patterns
 * - Ring complexity indices
 * - Aromatic system descriptors
 *
 * Total: 18 descriptors
 */

#include "cchem/compat.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

#define NUM_RC_DESCRIPTORS 18
#define MAX_RC_ATOMS 512

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

static int count_heavy_atoms(const molecule_t* mol) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) count++;
    }
    return count;
}

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

/* ============================================================================
 * Ring Size Distribution
 * ============================================================================ */

/* Count rings of specific size */
static int count_rings_of_size(const molecule_t* mol, int size) {
    int count = 0;
    for (int r = 0; r < mol->num_rings; r++) {
        if (mol->rings[r].size == size) count++;
    }
    return count;
}

/* Count aromatic rings */
static int count_aromatic_rings(const molecule_t* mol) {
    int count = 0;
    for (int r = 0; r < mol->num_rings; r++) {
        bool aromatic = true;
        for (int i = 0; i < mol->rings[r].size; i++) {
            if (!mol->atoms[mol->rings[r].atoms[i]].aromatic) {
                aromatic = false;
                break;
            }
        }
        if (aromatic) count++;
    }
    return count;
}

/* Count aliphatic rings */
static int count_aliphatic_rings(const molecule_t* mol) {
    int count = 0;
    for (int r = 0; r < mol->num_rings; r++) {
        bool has_non_aromatic = false;
        for (int i = 0; i < mol->rings[r].size; i++) {
            if (!mol->atoms[mol->rings[r].atoms[i]].aromatic) {
                has_non_aromatic = true;
                break;
            }
        }
        if (has_non_aromatic) count++;
    }
    return count;
}

/* Count heterocyclic rings (rings containing N, O, S) */
static int count_hetero_rings(const molecule_t* mol) {
    int count = 0;
    for (int r = 0; r < mol->num_rings; r++) {
        bool has_hetero = false;
        for (int i = 0; i < mol->rings[r].size; i++) {
            element_t e = mol->atoms[mol->rings[r].atoms[i]].element;
            if (e == ELEM_N || e == ELEM_O || e == ELEM_S) {
                has_hetero = true;
                break;
            }
        }
        if (has_hetero) count++;
    }
    return count;
}

/* ============================================================================
 * Ring Fusion Analysis
 * ============================================================================ */

/* Count fused ring atoms (atoms in multiple rings) */
static int count_fused_atoms(const molecule_t* mol) {
    if (mol->num_atoms > MAX_RC_ATOMS) return 0;

    int ring_count[MAX_RC_ATOMS] = {0};
    for (int r = 0; r < mol->num_rings; r++) {
        for (int i = 0; i < mol->rings[r].size; i++) {
            ring_count[mol->rings[r].atoms[i]]++;
        }
    }

    int fused = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (ring_count[i] > 1) fused++;
    }
    return fused;
}

/* Count fused ring bonds (bonds shared by multiple rings) */
static int count_fused_bonds(const molecule_t* mol) {
    int count = 0;
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        if (!bond->in_ring) continue;

        /* Check if both atoms are in multiple rings */
        if (mol->atoms[bond->atom1].ring_count > 1 &&
            mol->atoms[bond->atom2].ring_count > 1) {
            count++;
        }
    }
    return count;
}

/* Count spiro atoms (atoms connecting rings at single point) */
static int count_spiro_atoms(const molecule_t* mol) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        const atom_t* atom = &mol->atoms[i];
        if (atom->element == ELEM_H) continue;

        /* Spiro: degree 4, in 2+ rings, but rings share only this atom */
        if (atom->ring_count >= 2 && get_heavy_degree(mol, i) == 4) {
            /* Check if any neighbor is also in multiple rings */
            bool neighbor_shared = false;
            for (int j = 0; j < atom->num_neighbors; j++) {
                int nb = atom->neighbors[j];
                if (mol->atoms[nb].element != ELEM_H && mol->atoms[nb].ring_count >= 2) {
                    neighbor_shared = true;
                    break;
                }
            }
            if (!neighbor_shared) count++;
        }
    }
    return count;
}

/* Count bridgehead atoms (in 3+ rings) */
static int count_bridgehead_atoms(const molecule_t* mol) {
    int count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H && mol->atoms[i].ring_count >= 3) {
            count++;
        }
    }
    return count;
}

/* ============================================================================
 * Ring Complexity Indices
 * ============================================================================ */

/* Cyclomatic complexity (number of independent cycles) */
static int compute_cyclomatic(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    int n_heavy_bonds = 0;

    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->atoms[mol->bonds[b].atom1].element != ELEM_H &&
            mol->atoms[mol->bonds[b].atom2].element != ELEM_H) {
            n_heavy_bonds++;
        }
    }

    /* Cyclomatic number = edges - vertices + 1 (for connected graph) */
    return n_heavy_bonds - n_heavy + 1;
}

/* Ring complexity index (weighted sum) */
static double compute_ring_complexity(const molecule_t* mol) {
    double complexity = 0.0;

    for (int r = 0; r < mol->num_rings; r++) {
        int size = mol->rings[r].size;
        bool aromatic = true;
        bool has_hetero = false;

        for (int i = 0; i < size; i++) {
            int atom_idx = mol->rings[r].atoms[i];
            if (!mol->atoms[atom_idx].aromatic) aromatic = false;
            element_t e = mol->atoms[atom_idx].element;
            if (e == ELEM_N || e == ELEM_O || e == ELEM_S) has_hetero = true;
        }

        /* Weight: base by size, multiply by type */
        double weight = (double)size / 6.0;  /* Normalize to benzene */
        if (aromatic) weight *= 1.5;
        if (has_hetero) weight *= 1.2;

        complexity += weight;
    }

    return complexity;
}

/* Ring system complexity (accounts for fusion) */
static double compute_ring_system_complexity(const molecule_t* mol) {
    if (mol->num_rings == 0) return 0.0;

    double base = compute_ring_complexity(mol);
    int fused_atoms = count_fused_atoms(mol);
    int spiro_atoms = count_spiro_atoms(mol);

    /* Add fusion penalty */
    return base + 0.5 * fused_atoms + 0.3 * spiro_atoms;
}

/* Ring count ratio */
static double compute_ring_ratio(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy == 0) return 0.0;

    int ring_atoms = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H && mol->atoms[i].ring_count > 0) {
            ring_atoms++;
        }
    }

    return (double)ring_atoms / n_heavy;
}

/* Aromatic ratio */
static double compute_aromatic_ratio(const molecule_t* mol) {
    int n_heavy = count_heavy_atoms(mol);
    if (n_heavy == 0) return 0.0;

    int aromatic_atoms = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H && mol->atoms[i].aromatic) {
            aromatic_atoms++;
        }
    }

    return (double)aromatic_atoms / n_heavy;
}

/* ============================================================================
 * Batch Computation
 * ============================================================================ */

int descriptors_compute_ringcomplexity_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_RC_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    int idx = 0;

    /* Ring size counts */
    values[idx++].i = count_rings_of_size(mol, 3);   /* 0: Ring3Count */
    values[idx++].i = count_rings_of_size(mol, 4);   /* 1: Ring4Count */
    values[idx++].i = count_rings_of_size(mol, 5);   /* 2: Ring5Count */
    values[idx++].i = count_rings_of_size(mol, 6);   /* 3: Ring6Count */
    values[idx++].i = count_rings_of_size(mol, 7);   /* 4: Ring7Count */

    /* Ring type counts */
    values[idx++].i = count_aromatic_rings(mol);     /* 5: AromaticRingCount */
    values[idx++].i = count_aliphatic_rings(mol);    /* 6: AliphaticRingCount */
    values[idx++].i = count_hetero_rings(mol);       /* 7: HeterocyclicRingCount */

    /* Fusion counts */
    values[idx++].i = count_fused_atoms(mol);        /* 8: FusedRingAtoms */
    values[idx++].i = count_fused_bonds(mol);        /* 9: FusedRingBonds */
    values[idx++].i = count_spiro_atoms(mol);        /* 10: SpiroAtomCount */
    values[idx++].i = count_bridgehead_atoms(mol);   /* 11: BridgeheadAtomCount */

    /* Complexity indices */
    values[idx++].i = compute_cyclomatic(mol);       /* 12: CyclomaticNumber */
    values[idx++].d = compute_ring_complexity(mol);  /* 13: RingComplexityIndex */
    values[idx++].d = compute_ring_system_complexity(mol); /* 14: RingSystemComplexity */
    values[idx++].d = compute_ring_ratio(mol);       /* 15: RingAtomRatio */
    values[idx++].d = compute_aromatic_ratio(mol);   /* 16: AromaticAtomRatio */

    /* Total ring count */
    values[idx++].i = mol->num_rings;                /* 17: TotalRingCount */

    return NUM_RC_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* rc_cached_mol = NULL;
static _Thread_local descriptor_value_t rc_cached_values[NUM_RC_DESCRIPTORS];

static inline void ensure_rc_computed(const molecule_t* mol) {
    if (rc_cached_mol != mol) {
        descriptors_compute_ringcomplexity_all(mol, rc_cached_values);
        rc_cached_mol = mol;
    }
}

#define DEFINE_RC_FUNC_INT(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_rc_computed(mol); \
    value->i = rc_cached_values[idx].i; \
    return CCHEM_OK; \
}

#define DEFINE_RC_FUNC_DBL(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_rc_computed(mol); \
    value->d = rc_cached_values[idx].d; \
    return CCHEM_OK; \
}

DEFINE_RC_FUNC_INT(ring3_count, 0)
DEFINE_RC_FUNC_INT(ring4_count, 1)
DEFINE_RC_FUNC_INT(ring5_count, 2)
DEFINE_RC_FUNC_INT(ring6_count, 3)
DEFINE_RC_FUNC_INT(ring7_count, 4)
DEFINE_RC_FUNC_INT(aromatic_ring_count, 5)
DEFINE_RC_FUNC_INT(aliphatic_ring_count, 6)
DEFINE_RC_FUNC_INT(hetero_ring_count, 7)
DEFINE_RC_FUNC_INT(fused_ring_atoms, 8)
DEFINE_RC_FUNC_INT(fused_ring_bonds, 9)
DEFINE_RC_FUNC_INT(spiro_atom_count, 10)
DEFINE_RC_FUNC_INT(bridgehead_atom_count, 11)
DEFINE_RC_FUNC_INT(cyclomatic_number, 12)
DEFINE_RC_FUNC_DBL(ring_complexity_idx, 13)
DEFINE_RC_FUNC_DBL(ring_system_complexity, 14)
DEFINE_RC_FUNC_DBL(ring_atom_ratio, 15)
DEFINE_RC_FUNC_DBL(aromatic_atom_ratio, 16)
DEFINE_RC_FUNC_INT(total_ring_count, 17)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_RC_INT(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_COUNTS; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

#define REGISTER_RC_DBL(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_ringcomplexity(void) {
    /* Ring size counts */
    REGISTER_RC_INT("Ring3Count", "Number of 3-membered rings", desc_ring3_count);
    REGISTER_RC_INT("Ring4Count", "Number of 4-membered rings", desc_ring4_count);
    REGISTER_RC_INT("Ring5Count", "Number of 5-membered rings", desc_ring5_count);
    REGISTER_RC_INT("Ring6Count", "Number of 6-membered rings", desc_ring6_count);
    REGISTER_RC_INT("Ring7Count", "Number of 7-membered rings", desc_ring7_count);

    /* Ring type counts */
    REGISTER_RC_INT("AromaticRingCount", "Number of aromatic rings", desc_aromatic_ring_count);
    REGISTER_RC_INT("AliphaticRingCount", "Number of aliphatic rings", desc_aliphatic_ring_count);
    REGISTER_RC_INT("HeterocyclicRingCount", "Number of heterocyclic rings", desc_hetero_ring_count);

    /* Fusion counts */
    REGISTER_RC_INT("FusedRingAtoms", "Atoms in fused ring systems", desc_fused_ring_atoms);
    REGISTER_RC_INT("FusedRingBonds", "Bonds shared by rings", desc_fused_ring_bonds);
    REGISTER_RC_INT("SpiroAtomCount", "Spiro junction atoms", desc_spiro_atom_count);
    REGISTER_RC_INT("BridgeheadAtomCount", "Bridgehead atoms", desc_bridgehead_atom_count);

    /* Complexity indices */
    REGISTER_RC_INT("CyclomaticNumber", "Number of independent cycles", desc_cyclomatic_number);
    REGISTER_RC_DBL("RingComplexityIndex", "Weighted ring complexity", desc_ring_complexity_idx);
    REGISTER_RC_DBL("RingSystemComplexity", "Ring system complexity", desc_ring_system_complexity);
    REGISTER_RC_DBL("RingAtomRatio", "Fraction of ring atoms", desc_ring_atom_ratio);
    REGISTER_RC_DBL("AromaticAtomRatio", "Fraction of aromatic atoms", desc_aromatic_atom_ratio);

    /* Total */
    REGISTER_RC_INT("TotalRingCount", "Total number of rings", desc_total_ring_count);
}
