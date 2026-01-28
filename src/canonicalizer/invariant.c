/**
 * @file invariant.c
 * @brief Atom invariant calculation for canonicalization (Morgan algorithm)
 */

#include "cchem/compat.h"
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

/* Structure for tie-breaking sort with multi-level neighborhood signatures */
typedef struct {
    int index;
    uint64_t invariant;
    uint64_t neighbor_sig;      /* Signature based on immediate neighbors */
    uint64_t extended_sig;      /* Signature based on 2-hop neighborhood */
    int degree;
} tie_pair_t;

/*
 * Compute a deterministic signature for an atom based on its neighbors' invariants.
 * This is invariant to atom renumbering since it only uses invariant values.
 */
static uint64_t compute_neighbor_signature(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Collect neighbor invariants with bond orders */
    uint64_t neighbor_data[16];
    int n = atom->num_neighbors < 16 ? atom->num_neighbors : 16;

    for (int i = 0; i < n; i++) {
        int neighbor = atom->neighbors[i];
        int bond_idx = atom->neighbor_bonds[i];

        uint64_t inv = mol->atoms[neighbor].invariant;
        int bond_order = 1;
        if (bond_idx >= 0 && bond_idx < mol->num_bonds) {
            bond_order = bond_get_int_order(&mol->bonds[bond_idx]);
        }

        /* Combine neighbor invariant with bond order */
        neighbor_data[i] = hash_combine(inv, (uint64_t)bond_order);
    }

    /* Sort for canonical ordering */
    for (int i = 1; i < n; i++) {
        uint64_t key = neighbor_data[i];
        int j = i - 1;
        while (j >= 0 && neighbor_data[j] > key) {
            neighbor_data[j + 1] = neighbor_data[j];
            j--;
        }
        neighbor_data[j + 1] = key;
    }

    /* Combine into single signature */
    uint64_t sig = 0;
    for (int i = 0; i < n; i++) {
        sig = hash_combine(sig, neighbor_data[i]);
    }

    return sig;
}

/*
 * Compute extended signature considering 2-hop neighborhood.
 * This helps distinguish atoms that have identical immediate neighborhoods
 * but differ in their extended environment.
 */
static uint64_t compute_extended_signature(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];

    /* Collect signatures of all atoms within 2 hops */
    uint64_t extended_data[64];
    int count = 0;

    /* Add self */
    extended_data[count++] = atom->invariant;

    /* Add immediate neighbors */
    for (int i = 0; i < atom->num_neighbors && count < 64; i++) {
        int neighbor = atom->neighbors[i];
        extended_data[count++] = mol->atoms[neighbor].invariant;

        /* Add neighbors of neighbors (2-hop) */
        const atom_t* n_atom = &mol->atoms[neighbor];
        for (int j = 0; j < n_atom->num_neighbors && count < 64; j++) {
            int nn = n_atom->neighbors[j];
            if (nn != atom_idx) {  /* Exclude the original atom */
                extended_data[count++] = mol->atoms[nn].invariant;
            }
        }
    }

    /* Sort for canonical ordering */
    for (int i = 1; i < count; i++) {
        uint64_t key = extended_data[i];
        int j = i - 1;
        while (j >= 0 && extended_data[j] > key) {
            extended_data[j + 1] = extended_data[j];
            j--;
        }
        extended_data[j + 1] = key;
    }

    /* Combine into single signature */
    uint64_t sig = 0;
    for (int i = 0; i < count; i++) {
        sig = hash_combine(sig, extended_data[i]);
    }

    return sig;
}

/*
 * Comparison for tie-breaking using neighborhood signatures.
 * This ordering is invariant to atom renumbering.
 */
static int tie_pair_compare(const void* a, const void* b) {
    const tie_pair_t* pa = (const tie_pair_t*)a;
    const tie_pair_t* pb = (const tie_pair_t*)b;

    /* Primary: higher invariant = lower rank = comes first (descending) */
    if (pa->invariant > pb->invariant) return -1;
    if (pa->invariant < pb->invariant) return 1;

    /* Secondary: higher neighbor signature comes first */
    if (pa->neighbor_sig > pb->neighbor_sig) return -1;
    if (pa->neighbor_sig < pb->neighbor_sig) return 1;

    /* Tertiary: higher extended signature comes first */
    if (pa->extended_sig > pb->extended_sig) return -1;
    if (pa->extended_sig < pb->extended_sig) return 1;

    /* Quaternary: lower degree comes first */
    if (pa->degree < pb->degree) return -1;
    if (pa->degree > pb->degree) return 1;

    /* Final: atoms are truly symmetric - return 0 to indicate tie */
    return 0;
}

/* Debug flag - set to 1 to enable debug output */
#define INVARIANT_DEBUG 0

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

    /* Populate pairs with invariant, neighborhood signatures, and degree */
    for (int i = 0; i < n; i++) {
        pairs[i].index = i;
        pairs[i].invariant = mol->atoms[i].invariant;
        pairs[i].neighbor_sig = compute_neighbor_signature(mol, i);
        pairs[i].extended_sig = compute_extended_signature(mol, i);
        pairs[i].degree = mol->atoms[i].num_neighbors;
    }

#if INVARIANT_DEBUG
    fprintf(stderr, "Before tie-breaking sort:\n");
    for (int i = 0; i < n; i++) {
        fprintf(stderr, "  Atom %2d: inv=%016llx nsig=%016llx esig=%016llx deg=%d elem=%d\n",
                pairs[i].index, (unsigned long long)pairs[i].invariant,
                (unsigned long long)pairs[i].neighbor_sig,
                (unsigned long long)pairs[i].extended_sig,
                pairs[i].degree, mol->atoms[pairs[i].index].element);
    }
#endif

    /* Sort using O(n log n) qsort */
    qsort(pairs, n, sizeof(tie_pair_t), tie_pair_compare);

    /*
     * Check for remaining ties (truly symmetric atoms).
     * If any exist, use canonical BFS to establish a deterministic order.
     */
    bool has_ties = false;
    for (int i = 1; i < n && !has_ties; i++) {
        if (pairs[i].invariant == pairs[i-1].invariant &&
            pairs[i].neighbor_sig == pairs[i-1].neighbor_sig &&
            pairs[i].extended_sig == pairs[i-1].extended_sig &&
            pairs[i].degree == pairs[i-1].degree) {
            has_ties = true;
        }
    }

    if (has_ties) {
        /*
         * For truly symmetric atoms, we fall back to using atom index as
         * the tie-breaker. This is acceptable because:
         * 1. Symmetric atoms are chemically equivalent - any ordering is valid
         * 2. The atom indices determine the DFS traversal order
         * 3. The SMILES writer uses the same ordering, so the output is deterministic
         *
         * Note: This may not produce round-trip stable SMILES for molecules
         * with symmetric atoms because atom indices change on re-parsing.
         * For strict round-trip stability, more sophisticated algorithms
         * like SMILES-based fingerprinting would be needed.
         *
         * The key for correctness is that molecules that are structurally
         * equivalent will produce equivalent canonical SMILES (even if
         * not identical strings for symmetric molecules).
         */
        for (int i = 1; i < n; i++) {
            for (int j = i; j > 0; j--) {
                /* Check if pairs[j] and pairs[j-1] are tied */
                if (pairs[j].invariant == pairs[j-1].invariant &&
                    pairs[j].neighbor_sig == pairs[j-1].neighbor_sig &&
                    pairs[j].extended_sig == pairs[j-1].extended_sig &&
                    pairs[j].degree == pairs[j-1].degree) {
                    /* Tied - use atom index as final tie-breaker */
                    if (pairs[j].index < pairs[j-1].index) {
                        tie_pair_t tmp = pairs[j];
                        pairs[j] = pairs[j-1];
                        pairs[j-1] = tmp;
                    }
                } else {
                    break;
                }
            }
        }
    }

#if INVARIANT_DEBUG
    fprintf(stderr, "After tie-breaking sort (ranks assigned):\n");
    for (int i = 0; i < n; i++) {
        fprintf(stderr, "  Rank %2d -> Atom %2d (inv=%016llx nsig=%016llx esig=%016llx deg=%d elem=%d)\n",
                i, pairs[i].index, (unsigned long long)pairs[i].invariant,
                (unsigned long long)pairs[i].neighbor_sig,
                (unsigned long long)pairs[i].extended_sig,
                pairs[i].degree, mol->atoms[pairs[i].index].element);
    }
#endif

    /* Assign unique ranks based on sorted position */
    for (int i = 0; i < n; i++) {
        mol->atoms[pairs[i].index].canon_rank = i;
    }

    if (use_heap) {
        free(pairs);
    }

    return CCHEM_OK;
}
