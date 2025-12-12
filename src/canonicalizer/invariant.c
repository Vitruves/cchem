/**
 * @file invariant.c
 * @brief Atom invariant calculation for canonicalization (Morgan algorithm)
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/invariant.h"

/* Thread-local scratch for invariant operations (avoid malloc per call) */
#define INV_SCRATCH_MAX_ATOMS 1024
static __thread uint64_t tl_inv_scratch[INV_SCRATCH_MAX_ATOMS];
static __thread uint64_t tl_new_invariants[INV_SCRATCH_MAX_ATOMS];

/*
 * Initial invariant encoding (packed into 64-bit integer):
 * Bits 56-63: Number of connections (8 bits)
 * Bits 48-55: Number of attached hydrogens (8 bits)
 * Bits 40-47: Sign + formal charge (8 bits)
 * Bits 32-39: Atomic number (8 bits)
 * Bits 24-31: Mass difference from natural (8 bits)
 * Bits 16-23: Number of rings (8 bits)
 * Bits 8-15:  Aromatic flag (8 bits)
 * Bits 0-7:   Reserved
 */

uint64_t invariant_pack(const invariant_components_t* comp) {
    uint64_t inv = 0;

    inv |= ((uint64_t)(comp->num_connections & 0xFF)) << 56;
    inv |= ((uint64_t)(comp->num_h & 0xFF)) << 48;
    inv |= ((uint64_t)((comp->charge + 128) & 0xFF)) << 40;
    inv |= ((uint64_t)(comp->atomic_num & 0xFF)) << 32;
    inv |= ((uint64_t)(comp->mass & 0xFF)) << 24;
    inv |= ((uint64_t)(comp->num_rings & 0xFF)) << 16;
    inv |= ((uint64_t)(comp->aromatic ? 1 : 0)) << 8;

    return inv;
}

void invariant_unpack(uint64_t inv, invariant_components_t* comp) {
    comp->num_connections = (inv >> 56) & 0xFF;
    comp->num_h = (inv >> 48) & 0xFF;
    comp->charge = ((inv >> 40) & 0xFF) - 128;
    comp->atomic_num = (inv >> 32) & 0xFF;
    comp->mass = (inv >> 24) & 0xFF;
    comp->num_rings = (inv >> 16) & 0xFF;
    comp->aromatic = ((inv >> 8) & 0xFF) != 0;
}

uint64_t invariant_calc_atom(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = molecule_get_atom_const(mol, atom_idx);
    if (!atom) return 0;

    invariant_components_t comp;

    comp.num_connections = (uint8_t)atom->num_neighbors;
    comp.num_h = (uint8_t)(atom->implicit_h_count >= 0 ? atom->implicit_h_count : 0);
    comp.charge = (int8_t)atom->charge;
    comp.atomic_num = (uint8_t)atom->element;
    comp.mass = (uint8_t)(atom->isotope > 0 ? atom->isotope : 0);
    comp.num_rings = (uint8_t)atom->ring_count;
    comp.aromatic = atom->aromatic;

    return invariant_pack(&comp);
}

cchem_status_t invariant_calc_initial(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Ensure rings are computed for ring count */
    if (!mol->rings_computed) {
        cchem_status_t status = molecule_find_rings(mol);
        if (status != CCHEM_OK) return status;
    }

    /* Calculate initial invariant for each atom */
    for (int i = 0; i < mol->num_atoms; i++) {
        mol->atoms[i].invariant = invariant_calc_atom(mol, i);
    }

    return CCHEM_OK;
}

/*
 * Hash combine function for refinement.
 * Using a variant of boost::hash_combine.
 */
static uint64_t hash_combine(uint64_t seed, uint64_t value) {
    return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

/* Single refinement step using Morgan algorithm */
int invariant_refine_step(molecule_t* mol) {
    if (!mol || mol->num_atoms == 0) return 0;

    int n = mol->num_atoms;

    /* Use thread-local scratch for small molecules, heap for large */
    uint64_t* new_invariants;
    bool use_heap = (n > INV_SCRATCH_MAX_ATOMS);

    if (use_heap) {
        new_invariants = (uint64_t*)calloc(n, sizeof(uint64_t));
        if (!new_invariants) return -1;
    } else {
        new_invariants = tl_new_invariants;
    }

    /* Count unique values before */
    int unique_before = invariant_count_distinct(mol);

    /* For each atom, compute new invariant using hash combining */
    for (int i = 0; i < n; i++) {
        atom_t* atom = &mol->atoms[i];

        /* Start with self invariant */
        uint64_t hash = atom->invariant;

        /* Collect and sort neighbor invariants for deterministic ordering */
        uint64_t neighbor_hashes[16];  /* Max 16 neighbors typically */
        int num_neighbors = atom->num_neighbors < 16 ? atom->num_neighbors : 16;

        for (int j = 0; j < num_neighbors; j++) {
            int neighbor = atom->neighbors[j];
            int bond_idx = atom->neighbor_bonds[j];

            uint64_t neighbor_inv = mol->atoms[neighbor].invariant;

            /* Include bond order in the hash */
            int bond_order = 1;
            if (bond_idx >= 0 && bond_idx < mol->num_bonds) {
                bond_order = bond_get_int_order(&mol->bonds[bond_idx]);
            }

            /* Hash the neighbor invariant with bond order */
            neighbor_hashes[j] = hash_combine(neighbor_inv, (uint64_t)bond_order);
        }

        /* Sort neighbor hashes for consistent ordering (insertion sort for small n) */
        for (int j = 1; j < num_neighbors; j++) {
            uint64_t key = neighbor_hashes[j];
            int k = j - 1;
            while (k >= 0 && neighbor_hashes[k] > key) {
                neighbor_hashes[k + 1] = neighbor_hashes[k];
                k--;
            }
            neighbor_hashes[k + 1] = key;
        }

        /* Combine all neighbor hashes into final invariant */
        for (int j = 0; j < num_neighbors; j++) {
            hash = hash_combine(hash, neighbor_hashes[j]);
        }

        new_invariants[i] = hash;
    }

    /* Update invariants */
    for (int i = 0; i < n; i++) {
        mol->atoms[i].invariant = new_invariants[i];
    }

    if (use_heap) {
        free(new_invariants);
    }

    int unique_after = invariant_count_distinct(mol);

    /* Return change in number of unique values */
    return unique_after - unique_before;
}

cchem_status_t invariant_refine(molecule_t* mol, int max_iterations) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    if (max_iterations <= 0) {
        max_iterations = mol->num_atoms;  /* Default: at most n iterations */
    }

    /* Iterate until no change or max iterations */
    for (int iter = 0; iter < max_iterations; iter++) {
        int change = invariant_refine_step(mol);

        if (change < 0) {
            return CCHEM_ERROR_MEMORY;
        }

        /* Stop if no improvement in distinguishing power */
        if (change == 0) {
            break;
        }

        /* Stop if all atoms are unique */
        if (invariant_all_unique(mol)) {
            break;
        }
    }

    return CCHEM_OK;
}

/* Comparison function for sorting */
typedef struct {
    int index;
    uint64_t invariant;
} inv_pair_t;

/* Comparison for descending order (higher invariant = lower rank) */
static int inv_pair_compare_desc(const void* a, const void* b) {
    const inv_pair_t* pa = (const inv_pair_t*)a;
    const inv_pair_t* pb = (const inv_pair_t*)b;

    if (pa->invariant > pb->invariant) return -1;
    if (pa->invariant < pb->invariant) return 1;
    return 0;
}

cchem_status_t invariant_to_ranks(molecule_t* mol) {
    if (!mol || mol->num_atoms == 0) return CCHEM_ERROR_INVALID_INPUT;

    /* Create pairs of (index, invariant) */
    inv_pair_t* pairs = (inv_pair_t*)calloc(mol->num_atoms, sizeof(inv_pair_t));
    if (!pairs) return CCHEM_ERROR_MEMORY;

    for (int i = 0; i < mol->num_atoms; i++) {
        pairs[i].index = i;
        pairs[i].invariant = mol->atoms[i].invariant;
    }

    /* Sort by invariant (descending - higher invariant = lower rank) */
    qsort(pairs, mol->num_atoms, sizeof(inv_pair_t), inv_pair_compare_desc);

    /* Assign ranks (equal invariants get same rank) */
    int current_rank = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (i > 0 && pairs[i].invariant != pairs[i - 1].invariant) {
            current_rank = i;
        }
        mol->atoms[pairs[i].index].canon_rank = current_rank;
    }

    free(pairs);
    return CCHEM_OK;
}

/* Comparison for qsort */
static int uint64_compare(const void* a, const void* b) {
    uint64_t va = *(const uint64_t*)a;
    uint64_t vb = *(const uint64_t*)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

int invariant_count_distinct(const molecule_t* mol) {
    if (!mol || mol->num_atoms == 0) return 0;

    int n = mol->num_atoms;

    /* Use thread-local scratch for small molecules, heap for large */
    uint64_t* sorted;
    bool use_heap = (n > INV_SCRATCH_MAX_ATOMS);

    if (use_heap) {
        sorted = (uint64_t*)malloc(n * sizeof(uint64_t));
        if (!sorted) return 0;
    } else {
        sorted = tl_inv_scratch;
    }

    /* Copy invariants */
    for (int i = 0; i < n; i++) {
        sorted[i] = mol->atoms[i].invariant;
    }

    /* Sort O(n log n) */
    qsort(sorted, n, sizeof(uint64_t), uint64_compare);

    /* Count distinct in O(n) */
    int count = 1;
    for (int i = 1; i < n; i++) {
        if (sorted[i] != sorted[i - 1]) {
            count++;
        }
    }

    if (use_heap) {
        free(sorted);
    }

    return count;
}

bool invariant_all_unique(const molecule_t* mol) {
    if (!mol) return false;
    return invariant_count_distinct(mol) == mol->num_atoms;
}

/* Structure for tie-breaking sort - needs molecule pointer for degree lookup */
typedef struct {
    int index;
    uint64_t invariant;
    int degree;
} tie_pair_t;

/* Comparison for tie-breaking: descending by invariant, then by degree, then by index */
static int tie_pair_compare(const void* a, const void* b) {
    const tie_pair_t* pa = (const tie_pair_t*)a;
    const tie_pair_t* pb = (const tie_pair_t*)b;

    /* Higher invariant = lower rank = comes first (descending) */
    if (pa->invariant > pb->invariant) return -1;
    if (pa->invariant < pb->invariant) return 1;

    /* Same invariant: lower degree comes first */
    if (pa->degree < pb->degree) return -1;
    if (pa->degree > pb->degree) return 1;

    /* Same degree: lower index comes first */
    if (pa->index < pb->index) return -1;
    if (pa->index > pb->index) return 1;

    return 0;
}

cchem_status_t invariant_break_ties(molecule_t* mol) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Convert to ranks first */
    cchem_status_t status = invariant_to_ranks(mol);
    if (status != CCHEM_OK) return status;

    /* If all unique already, done */
    if (invariant_all_unique(mol)) {
        return CCHEM_OK;
    }

    int n = mol->num_atoms;

    /* Use stack for small molecules, heap for large */
    tie_pair_t stack_pairs[256];
    tie_pair_t* pairs;
    bool use_heap = (n > 256);

    if (use_heap) {
        pairs = (tie_pair_t*)calloc(n, sizeof(tie_pair_t));
        if (!pairs) return CCHEM_ERROR_MEMORY;
    } else {
        pairs = stack_pairs;
    }

    /* Populate pairs with invariant and degree */
    for (int i = 0; i < n; i++) {
        pairs[i].index = i;
        pairs[i].invariant = mol->atoms[i].invariant;
        pairs[i].degree = mol->atoms[i].num_neighbors;
    }

    /* Sort using O(n log n) qsort */
    qsort(pairs, n, sizeof(tie_pair_t), tie_pair_compare);

    /* Assign unique ranks based on sorted position */
    for (int i = 0; i < n; i++) {
        mol->atoms[pairs[i].index].canon_rank = i;
    }

    if (use_heap) {
        free(pairs);
    }

    return CCHEM_OK;
}
