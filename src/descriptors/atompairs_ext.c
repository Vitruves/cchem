/**
 * @file atompairs_ext.c
 * @brief Extended Atom Pair Descriptors - Detailed topological pair analysis
 *
 * Categories:
 * 1. Aromatic-Aliphatic Pairs (14 descriptors)
 *    - Aromatic-aromatic pairs at distance 1-7
 *    - Aromatic-aliphatic pairs at distance 1-7
 *
 * 2. Degree-Based Pairs (14 descriptors)
 *    - Terminal-terminal, terminal-branch, branch-branch pairs
 *
 * 3. Charge-Based Pairs (14 descriptors)
 *    - Positive-positive, negative-negative, positive-negative pairs
 *
 * 4. Ring-Based Pairs (14 descriptors)
 *    - Ring-ring, ring-chain, chain-chain pairs at various distances
 *
 * Total: 56 descriptors
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cchem/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* Total atom pair extension descriptors */
#define NUM_APEXT_DESCRIPTORS 56

/* Maximum atoms for stack allocation */
#define MAX_ATOMS_STACK 256

/* Maximum path length for BFS */
#define MAX_PATH_LEN 8

/* ============================================================================
 * BFS Distance Matrix Computation
 * ============================================================================ */

static void compute_distance_matrix(const molecule_t* mol,
                                     int n_heavy, int* reverse_map, int* dist) {
    /* Initialize distances to infinity */
    for (int i = 0; i < n_heavy * n_heavy; i++) {
        dist[i] = 9999;
    }
    for (int i = 0; i < n_heavy; i++) {
        dist[i * n_heavy + i] = 0;
    }

    /* Set distance 1 for bonded atoms */
    for (int b = 0; b < mol->num_bonds; b++) {
        const bond_t* bond = &mol->bonds[b];
        int i = reverse_map[bond->atom1];
        int j = reverse_map[bond->atom2];
        if (i >= 0 && j >= 0) {
            dist[i * n_heavy + j] = 1;
            dist[j * n_heavy + i] = 1;
        }
    }

    /* Floyd-Warshall for shortest paths */
    for (int k = 0; k < n_heavy; k++) {
        for (int i = 0; i < n_heavy; i++) {
            for (int j = 0; j < n_heavy; j++) {
                if (dist[i * n_heavy + k] + dist[k * n_heavy + j] < dist[i * n_heavy + j]) {
                    dist[i * n_heavy + j] = dist[i * n_heavy + k] + dist[k * n_heavy + j];
                }
            }
        }
    }
}

/* ============================================================================
 * Atom Classification Helpers
 * ============================================================================ */

/* Check if atom is aromatic (using atom's aromatic flag) */
static bool is_aromatic(const molecule_t* mol, int atom_idx) {
    return mol->atoms[atom_idx].aromatic;
}

/* Check if atom is in any ring */
static bool is_ring_atom(const molecule_t* mol, int atom_idx) {
    for (int r = 0; r < mol->num_rings && r < MAX_RINGS; r++) {
        const ring_t* ring = &mol->rings[r];
        for (int i = 0; i < ring->size; i++) {
            if (ring->atoms[i] == atom_idx) return true;
        }
    }
    return false;
}

/* Get atom degree (non-H neighbors) */
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

/* Check if terminal (degree 1) */
static bool is_terminal(const molecule_t* mol, int atom_idx) {
    return get_heavy_degree(mol, atom_idx) == 1;
}

/* Check if branch point (degree >= 3) */
static bool is_branch(const molecule_t* mol, int atom_idx) {
    return get_heavy_degree(mol, atom_idx) >= 3;
}

/* Estimate partial charge (simplified) */
static double get_partial_charge(const molecule_t* mol, int atom_idx) {
    const atom_t* atom = &mol->atoms[atom_idx];
    /* Use formal charge if set */
    if (atom->charge != 0) return atom->charge;

    /* Simple electronegativity-based estimate */
    element_t elem = atom->element;
    switch (elem) {
        case ELEM_N: return -0.3;
        case ELEM_O: return -0.4;
        case ELEM_F: return -0.3;
        case ELEM_Cl: return -0.2;
        case ELEM_S: return -0.1;
        case ELEM_C: return 0.0;
        default: return 0.0;
    }
}

/* ============================================================================
 * Main Computation
 * ============================================================================ */

int descriptors_compute_atompairs_ext_all(const molecule_t* mol, descriptor_value_t* values) {
    if (!mol || !values) return -1;

    for (int i = 0; i < NUM_APEXT_DESCRIPTORS; i++) {
        values[i].d = 0.0;
    }

    /* Count and index heavy atoms */
    int n_heavy = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) n_heavy++;
    }

    if (n_heavy < 2) return NUM_APEXT_DESCRIPTORS;

    /* Allocate arrays */
    bool heap_alloc = (n_heavy > MAX_ATOMS_STACK);
    int* heavy_idx;
    int* reverse_map;
    int* dist;
    bool* is_arom;
    bool* is_ring;
    bool* is_term;
    bool* is_br;
    double* charges;

    size_t arr_size = n_heavy * sizeof(int);
    size_t map_size = mol->num_atoms * sizeof(int);
    size_t dist_size = n_heavy * n_heavy * sizeof(int);
    size_t bool_size = n_heavy * sizeof(bool);
    size_t dbl_size = n_heavy * sizeof(double);

    if (heap_alloc) {
        heavy_idx = (int*)malloc(arr_size);
        reverse_map = (int*)malloc(map_size);
        dist = (int*)malloc(dist_size);
        is_arom = (bool*)malloc(bool_size);
        is_ring = (bool*)malloc(bool_size);
        is_term = (bool*)malloc(bool_size);
        is_br = (bool*)malloc(bool_size);
        charges = (double*)malloc(dbl_size);
        if (!heavy_idx || !reverse_map || !dist || !is_arom || !is_ring ||
            !is_term || !is_br || !charges) {
            free(heavy_idx); free(reverse_map); free(dist);
            free(is_arom); free(is_ring); free(is_term); free(is_br); free(charges);
            return -1;
        }
    } else {
        heavy_idx = (int*)alloca(arr_size);
        reverse_map = (int*)alloca(map_size);
        dist = (int*)alloca(dist_size);
        is_arom = (bool*)alloca(bool_size);
        is_ring = (bool*)alloca(bool_size);
        is_term = (bool*)alloca(bool_size);
        is_br = (bool*)alloca(bool_size);
        charges = (double*)alloca(dbl_size);
    }

    /* Build heavy atom index and reverse map */
    memset(reverse_map, -1, map_size);
    int idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element != ELEM_H) {
            heavy_idx[idx] = i;
            reverse_map[i] = idx;
            idx++;
        }
    }

    /* Compute distance matrix */
    compute_distance_matrix(mol, n_heavy, reverse_map, dist);

    /* Classify all heavy atoms */
    for (int i = 0; i < n_heavy; i++) {
        int atom_idx = heavy_idx[i];
        is_arom[i] = is_aromatic(mol, atom_idx);
        is_ring[i] = is_ring_atom(mol, atom_idx);
        is_term[i] = is_terminal(mol, atom_idx);
        is_br[i] = is_branch(mol, atom_idx);
        charges[i] = get_partial_charge(mol, atom_idx);
    }

    /* ====================================================================
     * Section 1: Aromatic-Aliphatic Pairs (14 descriptors, idx 0-13)
     * ==================================================================== */

    /* AP_AromArom at distance 1-7 */
    int arom_arom[7] = {0};
    /* AP_AromAlip at distance 1-7 */
    int arom_alip[7] = {0};

    for (int i = 0; i < n_heavy; i++) {
        for (int j = i + 1; j < n_heavy; j++) {
            int d = dist[i * n_heavy + j];
            if (d >= 1 && d <= 7) {
                if (is_arom[i] && is_arom[j]) {
                    arom_arom[d - 1]++;
                } else if ((is_arom[i] && !is_arom[j]) || (!is_arom[i] && is_arom[j])) {
                    arom_alip[d - 1]++;
                }
            }
        }
    }

    int out_idx = 0;
    for (int d = 0; d < 7; d++) values[out_idx++].d = arom_arom[d];
    for (int d = 0; d < 7; d++) values[out_idx++].d = arom_alip[d];

    /* ====================================================================
     * Section 2: Degree-Based Pairs (14 descriptors, idx 14-27)
     * ==================================================================== */

    /* Terminal-terminal at distance 1-4 */
    int term_term[4] = {0};
    /* Terminal-branch at distance 1-4 */
    int term_branch[4] = {0};
    /* Branch-branch at distance 1-3 */
    int branch_branch[3] = {0};
    /* Branch-core at distance 1-3 (core = degree 2) */
    int branch_core[3] = {0};

    for (int i = 0; i < n_heavy; i++) {
        for (int j = i + 1; j < n_heavy; j++) {
            int d = dist[i * n_heavy + j];
            if (d >= 1 && d <= 4) {
                if (is_term[i] && is_term[j]) {
                    term_term[d - 1]++;
                }
                if ((is_term[i] && is_br[j]) || (is_br[i] && is_term[j])) {
                    term_branch[d - 1]++;
                }
            }
            if (d >= 1 && d <= 3) {
                if (is_br[i] && is_br[j]) {
                    branch_branch[d - 1]++;
                }
                /* Core = not terminal and not branch */
                bool core_i = !is_term[i] && !is_br[i];
                bool core_j = !is_term[j] && !is_br[j];
                if ((is_br[i] && core_j) || (core_i && is_br[j])) {
                    branch_core[d - 1]++;
                }
            }
        }
    }

    for (int d = 0; d < 4; d++) values[out_idx++].d = term_term[d];
    for (int d = 0; d < 4; d++) values[out_idx++].d = term_branch[d];
    for (int d = 0; d < 3; d++) values[out_idx++].d = branch_branch[d];
    for (int d = 0; d < 3; d++) values[out_idx++].d = branch_core[d];

    /* ====================================================================
     * Section 3: Charge-Based Pairs (14 descriptors, idx 28-41)
     * ==================================================================== */

    /* Positive-positive at distance 1-4 */
    int pos_pos[4] = {0};
    /* Negative-negative at distance 1-4 */
    int neg_neg[4] = {0};
    /* Positive-negative at distance 1-6 */
    int pos_neg[6] = {0};

    for (int i = 0; i < n_heavy; i++) {
        for (int j = i + 1; j < n_heavy; j++) {
            int d = dist[i * n_heavy + j];
            bool pos_i = charges[i] > 0.05;
            bool neg_i = charges[i] < -0.05;
            bool pos_j = charges[j] > 0.05;
            bool neg_j = charges[j] < -0.05;

            if (d >= 1 && d <= 4) {
                if (pos_i && pos_j) pos_pos[d - 1]++;
                if (neg_i && neg_j) neg_neg[d - 1]++;
            }
            if (d >= 1 && d <= 6) {
                if ((pos_i && neg_j) || (neg_i && pos_j)) {
                    pos_neg[d - 1]++;
                }
            }
        }
    }

    for (int d = 0; d < 4; d++) values[out_idx++].d = pos_pos[d];
    for (int d = 0; d < 4; d++) values[out_idx++].d = neg_neg[d];
    for (int d = 0; d < 6; d++) values[out_idx++].d = pos_neg[d];

    /* ====================================================================
     * Section 4: Ring-Based Pairs (14 descriptors, idx 42-55)
     * ==================================================================== */

    /* Ring-ring at distance 1-5 */
    int ring_ring[5] = {0};
    /* Ring-chain at distance 1-5 */
    int ring_chain[5] = {0};
    /* Chain-chain at distance 1-4 */
    int chain_chain[4] = {0};

    for (int i = 0; i < n_heavy; i++) {
        for (int j = i + 1; j < n_heavy; j++) {
            int d = dist[i * n_heavy + j];

            if (d >= 1 && d <= 5) {
                if (is_ring[i] && is_ring[j]) {
                    ring_ring[d - 1]++;
                }
                if ((is_ring[i] && !is_ring[j]) || (!is_ring[i] && is_ring[j])) {
                    ring_chain[d - 1]++;
                }
            }
            if (d >= 1 && d <= 4) {
                if (!is_ring[i] && !is_ring[j]) {
                    chain_chain[d - 1]++;
                }
            }
        }
    }

    for (int d = 0; d < 5; d++) values[out_idx++].d = ring_ring[d];
    for (int d = 0; d < 5; d++) values[out_idx++].d = ring_chain[d];
    for (int d = 0; d < 4; d++) values[out_idx++].d = chain_chain[d];

    if (heap_alloc) {
        free(heavy_idx); free(reverse_map); free(dist);
        free(is_arom); free(is_ring); free(is_term); free(is_br); free(charges);
    }

    return NUM_APEXT_DESCRIPTORS;
}

/* ============================================================================
 * Individual Descriptor Functions with Caching
 * ============================================================================ */

static _Thread_local const molecule_t* apext_cached_mol = NULL;
static _Thread_local uint64_t apext_cached_gen = 0;
static _Thread_local descriptor_value_t apext_cached_values[NUM_APEXT_DESCRIPTORS];

static inline void ensure_apext_computed(const molecule_t* mol) {
    uint64_t current_gen = descriptor_cache_generation();
    if (apext_cached_mol != mol || apext_cached_gen != current_gen) {
        descriptors_compute_atompairs_ext_all(mol, apext_cached_values);
        apext_cached_mol = mol;
        apext_cached_gen = current_gen;
    }
}

#define DEFINE_APEXT_FUNC(name, idx) \
static cchem_status_t desc_##name(const molecule_t* mol, descriptor_value_t* value) { \
    if (!mol || !value) return CCHEM_ERROR_INVALID_INPUT; \
    ensure_apext_computed(mol); \
    value->d = apext_cached_values[idx].d; \
    return CCHEM_OK; \
}

/* Aromatic-Aromatic pairs (7) */
DEFINE_APEXT_FUNC(ap_aromarom1, 0)
DEFINE_APEXT_FUNC(ap_aromarom2, 1)
DEFINE_APEXT_FUNC(ap_aromarom3, 2)
DEFINE_APEXT_FUNC(ap_aromarom4, 3)
DEFINE_APEXT_FUNC(ap_aromarom5, 4)
DEFINE_APEXT_FUNC(ap_aromarom6, 5)
DEFINE_APEXT_FUNC(ap_aromarom7, 6)

/* Aromatic-Aliphatic pairs (7) */
DEFINE_APEXT_FUNC(ap_aromalip1, 7)
DEFINE_APEXT_FUNC(ap_aromalip2, 8)
DEFINE_APEXT_FUNC(ap_aromalip3, 9)
DEFINE_APEXT_FUNC(ap_aromalip4, 10)
DEFINE_APEXT_FUNC(ap_aromalip5, 11)
DEFINE_APEXT_FUNC(ap_aromalip6, 12)
DEFINE_APEXT_FUNC(ap_aromalip7, 13)

/* Terminal-Terminal pairs (4) */
DEFINE_APEXT_FUNC(ap_term1, 14)
DEFINE_APEXT_FUNC(ap_term2, 15)
DEFINE_APEXT_FUNC(ap_term3, 16)
DEFINE_APEXT_FUNC(ap_term4, 17)

/* Terminal-Branch pairs (4) */
DEFINE_APEXT_FUNC(ap_termbranch1, 18)
DEFINE_APEXT_FUNC(ap_termbranch2, 19)
DEFINE_APEXT_FUNC(ap_termbranch3, 20)
DEFINE_APEXT_FUNC(ap_termbranch4, 21)

/* Branch-Branch pairs (3) */
DEFINE_APEXT_FUNC(ap_branch1, 22)
DEFINE_APEXT_FUNC(ap_branch2, 23)
DEFINE_APEXT_FUNC(ap_branch3, 24)

/* Branch-Core pairs (3) */
DEFINE_APEXT_FUNC(ap_branchcore1, 25)
DEFINE_APEXT_FUNC(ap_branchcore2, 26)
DEFINE_APEXT_FUNC(ap_branchcore3, 27)

/* Positive-Positive pairs (4) */
DEFINE_APEXT_FUNC(ap_pospos1, 28)
DEFINE_APEXT_FUNC(ap_pospos2, 29)
DEFINE_APEXT_FUNC(ap_pospos3, 30)
DEFINE_APEXT_FUNC(ap_pospos4, 31)

/* Negative-Negative pairs (4) */
DEFINE_APEXT_FUNC(ap_negneg1, 32)
DEFINE_APEXT_FUNC(ap_negneg2, 33)
DEFINE_APEXT_FUNC(ap_negneg3, 34)
DEFINE_APEXT_FUNC(ap_negneg4, 35)

/* Positive-Negative pairs (6) */
DEFINE_APEXT_FUNC(ap_posneg1, 36)
DEFINE_APEXT_FUNC(ap_posneg2, 37)
DEFINE_APEXT_FUNC(ap_posneg3, 38)
DEFINE_APEXT_FUNC(ap_posneg4, 39)
DEFINE_APEXT_FUNC(ap_posneg5, 40)
DEFINE_APEXT_FUNC(ap_posneg6, 41)

/* Ring-Ring pairs (5) */
DEFINE_APEXT_FUNC(ap_ringring1, 42)
DEFINE_APEXT_FUNC(ap_ringring2, 43)
DEFINE_APEXT_FUNC(ap_ringring3, 44)
DEFINE_APEXT_FUNC(ap_ringring4, 45)
DEFINE_APEXT_FUNC(ap_ringring5, 46)

/* Ring-Chain pairs (5) */
DEFINE_APEXT_FUNC(ap_ringchain1, 47)
DEFINE_APEXT_FUNC(ap_ringchain2, 48)
DEFINE_APEXT_FUNC(ap_ringchain3, 49)
DEFINE_APEXT_FUNC(ap_ringchain4, 50)
DEFINE_APEXT_FUNC(ap_ringchain5, 51)

/* Chain-Chain pairs (4) */
DEFINE_APEXT_FUNC(ap_chainchain1, 52)
DEFINE_APEXT_FUNC(ap_chainchain2, 53)
DEFINE_APEXT_FUNC(ap_chainchain3, 54)
DEFINE_APEXT_FUNC(ap_chainchain4, 55)

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_APEXT(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_PROPERTIES; \
    def.value_type = DESC_VALUE_DOUBLE; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_atompairs_ext(void) {
    /* Aromatic-Aromatic pairs */
    REGISTER_APEXT("AP_AromArom1", "Aromatic-aromatic pairs at distance 1", desc_ap_aromarom1);
    REGISTER_APEXT("AP_AromArom2", "Aromatic-aromatic pairs at distance 2", desc_ap_aromarom2);
    REGISTER_APEXT("AP_AromArom3", "Aromatic-aromatic pairs at distance 3", desc_ap_aromarom3);
    REGISTER_APEXT("AP_AromArom4", "Aromatic-aromatic pairs at distance 4", desc_ap_aromarom4);
    REGISTER_APEXT("AP_AromArom5", "Aromatic-aromatic pairs at distance 5", desc_ap_aromarom5);
    REGISTER_APEXT("AP_AromArom6", "Aromatic-aromatic pairs at distance 6", desc_ap_aromarom6);
    REGISTER_APEXT("AP_AromArom7", "Aromatic-aromatic pairs at distance 7", desc_ap_aromarom7);

    /* Aromatic-Aliphatic pairs */
    REGISTER_APEXT("AP_AromAlip1", "Aromatic-aliphatic pairs at distance 1", desc_ap_aromalip1);
    REGISTER_APEXT("AP_AromAlip2", "Aromatic-aliphatic pairs at distance 2", desc_ap_aromalip2);
    REGISTER_APEXT("AP_AromAlip3", "Aromatic-aliphatic pairs at distance 3", desc_ap_aromalip3);
    REGISTER_APEXT("AP_AromAlip4", "Aromatic-aliphatic pairs at distance 4", desc_ap_aromalip4);
    REGISTER_APEXT("AP_AromAlip5", "Aromatic-aliphatic pairs at distance 5", desc_ap_aromalip5);
    REGISTER_APEXT("AP_AromAlip6", "Aromatic-aliphatic pairs at distance 6", desc_ap_aromalip6);
    REGISTER_APEXT("AP_AromAlip7", "Aromatic-aliphatic pairs at distance 7", desc_ap_aromalip7);

    /* Terminal-Terminal pairs */
    REGISTER_APEXT("AP_Term1", "Terminal-terminal pairs at distance 1", desc_ap_term1);
    REGISTER_APEXT("AP_Term2", "Terminal-terminal pairs at distance 2", desc_ap_term2);
    REGISTER_APEXT("AP_Term3", "Terminal-terminal pairs at distance 3", desc_ap_term3);
    REGISTER_APEXT("AP_Term4", "Terminal-terminal pairs at distance 4", desc_ap_term4);

    /* Terminal-Branch pairs */
    REGISTER_APEXT("AP_TermBranch1", "Terminal-branch pairs at distance 1", desc_ap_termbranch1);
    REGISTER_APEXT("AP_TermBranch2", "Terminal-branch pairs at distance 2", desc_ap_termbranch2);
    REGISTER_APEXT("AP_TermBranch3", "Terminal-branch pairs at distance 3", desc_ap_termbranch3);
    REGISTER_APEXT("AP_TermBranch4", "Terminal-branch pairs at distance 4", desc_ap_termbranch4);

    /* Branch-Branch pairs */
    REGISTER_APEXT("AP_Branch1", "Branch-branch pairs at distance 1", desc_ap_branch1);
    REGISTER_APEXT("AP_Branch2", "Branch-branch pairs at distance 2", desc_ap_branch2);
    REGISTER_APEXT("AP_Branch3", "Branch-branch pairs at distance 3", desc_ap_branch3);

    /* Branch-Core pairs */
    REGISTER_APEXT("AP_BranchCore1", "Branch-core pairs at distance 1", desc_ap_branchcore1);
    REGISTER_APEXT("AP_BranchCore2", "Branch-core pairs at distance 2", desc_ap_branchcore2);
    REGISTER_APEXT("AP_BranchCore3", "Branch-core pairs at distance 3", desc_ap_branchcore3);

    /* Positive-Positive pairs */
    REGISTER_APEXT("AP_PosPos1", "Positive-positive pairs at distance 1", desc_ap_pospos1);
    REGISTER_APEXT("AP_PosPos2", "Positive-positive pairs at distance 2", desc_ap_pospos2);
    REGISTER_APEXT("AP_PosPos3", "Positive-positive pairs at distance 3", desc_ap_pospos3);
    REGISTER_APEXT("AP_PosPos4", "Positive-positive pairs at distance 4", desc_ap_pospos4);

    /* Negative-Negative pairs */
    REGISTER_APEXT("AP_NegNeg1", "Negative-negative pairs at distance 1", desc_ap_negneg1);
    REGISTER_APEXT("AP_NegNeg2", "Negative-negative pairs at distance 2", desc_ap_negneg2);
    REGISTER_APEXT("AP_NegNeg3", "Negative-negative pairs at distance 3", desc_ap_negneg3);
    REGISTER_APEXT("AP_NegNeg4", "Negative-negative pairs at distance 4", desc_ap_negneg4);

    /* Positive-Negative pairs */
    REGISTER_APEXT("AP_PosNeg1", "Positive-negative pairs at distance 1", desc_ap_posneg1);
    REGISTER_APEXT("AP_PosNeg2", "Positive-negative pairs at distance 2", desc_ap_posneg2);
    REGISTER_APEXT("AP_PosNeg3", "Positive-negative pairs at distance 3", desc_ap_posneg3);
    REGISTER_APEXT("AP_PosNeg4", "Positive-negative pairs at distance 4", desc_ap_posneg4);
    REGISTER_APEXT("AP_PosNeg5", "Positive-negative pairs at distance 5", desc_ap_posneg5);
    REGISTER_APEXT("AP_PosNeg6", "Positive-negative pairs at distance 6", desc_ap_posneg6);

    /* Ring-Ring pairs */
    REGISTER_APEXT("AP_RingRing1", "Ring-ring pairs at distance 1", desc_ap_ringring1);
    REGISTER_APEXT("AP_RingRing2", "Ring-ring pairs at distance 2", desc_ap_ringring2);
    REGISTER_APEXT("AP_RingRing3", "Ring-ring pairs at distance 3", desc_ap_ringring3);
    REGISTER_APEXT("AP_RingRing4", "Ring-ring pairs at distance 4", desc_ap_ringring4);
    REGISTER_APEXT("AP_RingRing5", "Ring-ring pairs at distance 5", desc_ap_ringring5);

    /* Ring-Chain pairs */
    REGISTER_APEXT("AP_RingChain1", "Ring-chain pairs at distance 1", desc_ap_ringchain1);
    REGISTER_APEXT("AP_RingChain2", "Ring-chain pairs at distance 2", desc_ap_ringchain2);
    REGISTER_APEXT("AP_RingChain3", "Ring-chain pairs at distance 3", desc_ap_ringchain3);
    REGISTER_APEXT("AP_RingChain4", "Ring-chain pairs at distance 4", desc_ap_ringchain4);
    REGISTER_APEXT("AP_RingChain5", "Ring-chain pairs at distance 5", desc_ap_ringchain5);

    /* Chain-Chain pairs */
    REGISTER_APEXT("AP_ChainChain1", "Chain-chain pairs at distance 1", desc_ap_chainchain1);
    REGISTER_APEXT("AP_ChainChain2", "Chain-chain pairs at distance 2", desc_ap_chainchain2);
    REGISTER_APEXT("AP_ChainChain3", "Chain-chain pairs at distance 3", desc_ap_chainchain3);
    REGISTER_APEXT("AP_ChainChain4", "Chain-chain pairs at distance 4", desc_ap_chainchain4);
}
