/**
 * @file atompairs.c
 * @brief Topological atom pair descriptors
 *
 * Counts pairs of atom types at various topological distances.
 * Fast O(n^2) BFS-based implementation.
 *
 * Atom types: C, N, O, S, Hal (F/Cl/Br/I), Other
 * Distances: 1-7 bonds
 */

#include <string.h>
#include <stdbool.h>
#include "cchem/utils/descriptors.h"
#include "cchem/canonicalizer/molecule.h"
#include "cchem/canonicalizer/atom.h"
#include "cchem/canonicalizer/bond.h"

/* ============================================================================
 * Distance Matrix Computation (BFS)
 * ============================================================================ */

#define MAX_DIST 8
#define MAX_AP_ATOMS 256

/* Compute shortest path distances using BFS */
static void compute_distances(const molecule_t* mol, int dist[MAX_AP_ATOMS][MAX_AP_ATOMS]) {
    int n = mol->num_atoms < MAX_AP_ATOMS ? mol->num_atoms : MAX_AP_ATOMS;

    /* Initialize with -1 (unreachable) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dist[i][j] = (i == j) ? 0 : -1;
        }
    }

    /* BFS from each atom */
    for (int start = 0; start < n; start++) {
        if (mol->atoms[start].element == ELEM_H) continue;

        int queue[MAX_AP_ATOMS];
        int front = 0, back = 0;
        queue[back++] = start;

        while (front < back) {
            int curr = queue[front++];
            int curr_dist = dist[start][curr];

            if (curr_dist >= MAX_DIST - 1) continue;

            const atom_t* atom = &mol->atoms[curr];
            for (int i = 0; i < atom->num_neighbors; i++) {
                int next = atom->neighbors[i];
                if (next >= n) continue;
                if (mol->atoms[next].element == ELEM_H) continue;

                if (dist[start][next] < 0) {
                    dist[start][next] = curr_dist + 1;
                    dist[next][start] = curr_dist + 1;
                    queue[back++] = next;
                }
            }
        }
    }
}

/* Atom type classification */
typedef enum {
    ATYPE_C = 0,
    ATYPE_N = 1,
    ATYPE_O = 2,
    ATYPE_S = 3,
    ATYPE_HAL = 4,
    ATYPE_OTHER = 5,
    NUM_ATYPES = 6
} atom_type_t;

static atom_type_t get_atom_type(element_t elem) {
    switch (elem) {
        case ELEM_C:  return ATYPE_C;
        case ELEM_N:  return ATYPE_N;
        case ELEM_O:  return ATYPE_O;
        case ELEM_S:  return ATYPE_S;
        case ELEM_F:
        case ELEM_Cl:
        case ELEM_Br:
        case ELEM_I:  return ATYPE_HAL;
        default:      return ATYPE_OTHER;
    }
}

/* ============================================================================
 * Atom Pair Count Structure
 * ============================================================================ */

/* We'll compute all atom pairs in one pass and cache results */
typedef struct {
    int pairs[NUM_ATYPES][NUM_ATYPES][MAX_DIST];  /* [type1][type2][dist] */
    bool computed;
} atom_pair_cache_t;

static void compute_atom_pairs(const molecule_t* mol, atom_pair_cache_t* cache) {
    if (cache->computed) return;

    memset(cache->pairs, 0, sizeof(cache->pairs));

    /* Compute distance matrix */
    int dist[MAX_AP_ATOMS][MAX_AP_ATOMS];
    compute_distances(mol, dist);

    int n = mol->num_atoms < MAX_AP_ATOMS ? mol->num_atoms : MAX_AP_ATOMS;

    /* Count pairs */
    for (int i = 0; i < n; i++) {
        if (mol->atoms[i].element == ELEM_H) continue;
        atom_type_t ti = get_atom_type(mol->atoms[i].element);

        for (int j = i + 1; j < n; j++) {
            if (mol->atoms[j].element == ELEM_H) continue;
            if (dist[i][j] < 1 || dist[i][j] >= MAX_DIST) continue;

            atom_type_t tj = get_atom_type(mol->atoms[j].element);
            int d = dist[i][j];

            /* Count both orderings for symmetric access */
            cache->pairs[ti][tj][d]++;
            if (ti != tj) {
                cache->pairs[tj][ti][d]++;
            }
        }
    }

    cache->computed = true;
}

/* ============================================================================
 * Individual Atom Pair Descriptors
 * ============================================================================ */

/* Macro to generate pair count functions */
#define AP_FUNC(name, t1, t2, d) \
static cchem_status_t name(const molecule_t* mol, descriptor_value_t* value) { \
    atom_pair_cache_t cache = {0}; \
    compute_atom_pairs(mol, &cache); \
    value->i = cache.pairs[t1][t2][d]; \
    return CCHEM_OK; \
}

/* C-C pairs at distances 1-7 */
AP_FUNC(ap_CC1, ATYPE_C, ATYPE_C, 1)
AP_FUNC(ap_CC2, ATYPE_C, ATYPE_C, 2)
AP_FUNC(ap_CC3, ATYPE_C, ATYPE_C, 3)
AP_FUNC(ap_CC4, ATYPE_C, ATYPE_C, 4)
AP_FUNC(ap_CC5, ATYPE_C, ATYPE_C, 5)
AP_FUNC(ap_CC6, ATYPE_C, ATYPE_C, 6)
AP_FUNC(ap_CC7, ATYPE_C, ATYPE_C, 7)

/* C-N pairs */
AP_FUNC(ap_CN1, ATYPE_C, ATYPE_N, 1)
AP_FUNC(ap_CN2, ATYPE_C, ATYPE_N, 2)
AP_FUNC(ap_CN3, ATYPE_C, ATYPE_N, 3)
AP_FUNC(ap_CN4, ATYPE_C, ATYPE_N, 4)
AP_FUNC(ap_CN5, ATYPE_C, ATYPE_N, 5)
AP_FUNC(ap_CN6, ATYPE_C, ATYPE_N, 6)
AP_FUNC(ap_CN7, ATYPE_C, ATYPE_N, 7)

/* C-O pairs */
AP_FUNC(ap_CO1, ATYPE_C, ATYPE_O, 1)
AP_FUNC(ap_CO2, ATYPE_C, ATYPE_O, 2)
AP_FUNC(ap_CO3, ATYPE_C, ATYPE_O, 3)
AP_FUNC(ap_CO4, ATYPE_C, ATYPE_O, 4)
AP_FUNC(ap_CO5, ATYPE_C, ATYPE_O, 5)

/* C-S pairs */
AP_FUNC(ap_CS1, ATYPE_C, ATYPE_S, 1)
AP_FUNC(ap_CS2, ATYPE_C, ATYPE_S, 2)
AP_FUNC(ap_CS3, ATYPE_C, ATYPE_S, 3)

/* C-Hal pairs */
AP_FUNC(ap_CHal1, ATYPE_C, ATYPE_HAL, 1)
AP_FUNC(ap_CHal2, ATYPE_C, ATYPE_HAL, 2)
AP_FUNC(ap_CHal3, ATYPE_C, ATYPE_HAL, 3)

/* N-N pairs */
AP_FUNC(ap_NN1, ATYPE_N, ATYPE_N, 1)
AP_FUNC(ap_NN2, ATYPE_N, ATYPE_N, 2)
AP_FUNC(ap_NN3, ATYPE_N, ATYPE_N, 3)
AP_FUNC(ap_NN4, ATYPE_N, ATYPE_N, 4)
AP_FUNC(ap_NN5, ATYPE_N, ATYPE_N, 5)

/* N-O pairs */
AP_FUNC(ap_NO1, ATYPE_N, ATYPE_O, 1)
AP_FUNC(ap_NO2, ATYPE_N, ATYPE_O, 2)
AP_FUNC(ap_NO3, ATYPE_N, ATYPE_O, 3)
AP_FUNC(ap_NO4, ATYPE_N, ATYPE_O, 4)

/* O-O pairs */
AP_FUNC(ap_OO1, ATYPE_O, ATYPE_O, 1)
AP_FUNC(ap_OO2, ATYPE_O, ATYPE_O, 2)
AP_FUNC(ap_OO3, ATYPE_O, ATYPE_O, 3)
AP_FUNC(ap_OO4, ATYPE_O, ATYPE_O, 4)

/* S-S pairs */
AP_FUNC(ap_SS1, ATYPE_S, ATYPE_S, 1)
AP_FUNC(ap_SS2, ATYPE_S, ATYPE_S, 2)

/* N-S pairs */
AP_FUNC(ap_NS1, ATYPE_N, ATYPE_S, 1)
AP_FUNC(ap_NS2, ATYPE_N, ATYPE_S, 2)
AP_FUNC(ap_NS3, ATYPE_N, ATYPE_S, 3)

/* O-S pairs */
AP_FUNC(ap_OS1, ATYPE_O, ATYPE_S, 1)
AP_FUNC(ap_OS2, ATYPE_O, ATYPE_S, 2)

/* Hal-Hal pairs */
AP_FUNC(ap_HalHal1, ATYPE_HAL, ATYPE_HAL, 1)
AP_FUNC(ap_HalHal2, ATYPE_HAL, ATYPE_HAL, 2)
AP_FUNC(ap_HalHal3, ATYPE_HAL, ATYPE_HAL, 3)

/* ============================================================================
 * Summary Descriptors
 * ============================================================================ */

/* Total atom pairs at each distance */
static cchem_status_t ap_total_d1(const molecule_t* mol, descriptor_value_t* value) {
    atom_pair_cache_t cache = {0};
    compute_atom_pairs(mol, &cache);
    int sum = 0;
    for (int i = 0; i < NUM_ATYPES; i++) {
        for (int j = i; j < NUM_ATYPES; j++) {
            sum += cache.pairs[i][j][1];
        }
    }
    value->i = sum;
    return CCHEM_OK;
}

static cchem_status_t ap_total_d2(const molecule_t* mol, descriptor_value_t* value) {
    atom_pair_cache_t cache = {0};
    compute_atom_pairs(mol, &cache);
    int sum = 0;
    for (int i = 0; i < NUM_ATYPES; i++) {
        for (int j = i; j < NUM_ATYPES; j++) {
            sum += cache.pairs[i][j][2];
        }
    }
    value->i = sum;
    return CCHEM_OK;
}

static cchem_status_t ap_total_d3(const molecule_t* mol, descriptor_value_t* value) {
    atom_pair_cache_t cache = {0};
    compute_atom_pairs(mol, &cache);
    int sum = 0;
    for (int i = 0; i < NUM_ATYPES; i++) {
        for (int j = i; j < NUM_ATYPES; j++) {
            sum += cache.pairs[i][j][3];
        }
    }
    value->i = sum;
    return CCHEM_OK;
}

static cchem_status_t ap_total_d4(const molecule_t* mol, descriptor_value_t* value) {
    atom_pair_cache_t cache = {0};
    compute_atom_pairs(mol, &cache);
    int sum = 0;
    for (int i = 0; i < NUM_ATYPES; i++) {
        for (int j = i; j < NUM_ATYPES; j++) {
            sum += cache.pairs[i][j][4];
        }
    }
    value->i = sum;
    return CCHEM_OK;
}

static cchem_status_t ap_total_d5(const molecule_t* mol, descriptor_value_t* value) {
    atom_pair_cache_t cache = {0};
    compute_atom_pairs(mol, &cache);
    int sum = 0;
    for (int i = 0; i < NUM_ATYPES; i++) {
        for (int j = i; j < NUM_ATYPES; j++) {
            sum += cache.pairs[i][j][5];
        }
    }
    value->i = sum;
    return CCHEM_OK;
}

/* Heteroatom pairs */
static cchem_status_t ap_hetero_total(const molecule_t* mol, descriptor_value_t* value) {
    atom_pair_cache_t cache = {0};
    compute_atom_pairs(mol, &cache);
    int sum = 0;
    /* Sum N-N, N-O, N-S, O-O, O-S, S-S at all distances */
    for (int d = 1; d < MAX_DIST; d++) {
        sum += cache.pairs[ATYPE_N][ATYPE_N][d];
        sum += cache.pairs[ATYPE_N][ATYPE_O][d];
        sum += cache.pairs[ATYPE_N][ATYPE_S][d];
        sum += cache.pairs[ATYPE_O][ATYPE_O][d];
        sum += cache.pairs[ATYPE_O][ATYPE_S][d];
        sum += cache.pairs[ATYPE_S][ATYPE_S][d];
    }
    value->i = sum;
    return CCHEM_OK;
}

/* ============================================================================
 * Registration
 * ============================================================================ */

#define REGISTER_AP(name_str, desc_str, func) do { \
    descriptor_def_t def = {0}; \
    strncpy(def.name, name_str, MAX_DESCRIPTOR_NAME - 1); \
    strncpy(def.description, desc_str, sizeof(def.description) - 1); \
    def.category = DESC_CATEGORY_COUNTS; \
    def.value_type = DESC_VALUE_INT; \
    def.compute = func; \
    descriptor_register(&def); \
} while(0)

void descriptors_register_atompairs(void) {
    /* C-C pairs */
    REGISTER_AP("AP_CC1", "C-C pairs at distance 1", ap_CC1);
    REGISTER_AP("AP_CC2", "C-C pairs at distance 2", ap_CC2);
    REGISTER_AP("AP_CC3", "C-C pairs at distance 3", ap_CC3);
    REGISTER_AP("AP_CC4", "C-C pairs at distance 4", ap_CC4);
    REGISTER_AP("AP_CC5", "C-C pairs at distance 5", ap_CC5);
    REGISTER_AP("AP_CC6", "C-C pairs at distance 6", ap_CC6);
    REGISTER_AP("AP_CC7", "C-C pairs at distance 7", ap_CC7);

    /* C-N pairs */
    REGISTER_AP("AP_CN1", "C-N pairs at distance 1", ap_CN1);
    REGISTER_AP("AP_CN2", "C-N pairs at distance 2", ap_CN2);
    REGISTER_AP("AP_CN3", "C-N pairs at distance 3", ap_CN3);
    REGISTER_AP("AP_CN4", "C-N pairs at distance 4", ap_CN4);
    REGISTER_AP("AP_CN5", "C-N pairs at distance 5", ap_CN5);
    REGISTER_AP("AP_CN6", "C-N pairs at distance 6", ap_CN6);
    REGISTER_AP("AP_CN7", "C-N pairs at distance 7", ap_CN7);

    /* C-O pairs */
    REGISTER_AP("AP_CO1", "C-O pairs at distance 1", ap_CO1);
    REGISTER_AP("AP_CO2", "C-O pairs at distance 2", ap_CO2);
    REGISTER_AP("AP_CO3", "C-O pairs at distance 3", ap_CO3);
    REGISTER_AP("AP_CO4", "C-O pairs at distance 4", ap_CO4);
    REGISTER_AP("AP_CO5", "C-O pairs at distance 5", ap_CO5);

    /* C-S pairs */
    REGISTER_AP("AP_CS1", "C-S pairs at distance 1", ap_CS1);
    REGISTER_AP("AP_CS2", "C-S pairs at distance 2", ap_CS2);
    REGISTER_AP("AP_CS3", "C-S pairs at distance 3", ap_CS3);

    /* C-Hal pairs */
    REGISTER_AP("AP_CHal1", "C-Hal pairs at distance 1", ap_CHal1);
    REGISTER_AP("AP_CHal2", "C-Hal pairs at distance 2", ap_CHal2);
    REGISTER_AP("AP_CHal3", "C-Hal pairs at distance 3", ap_CHal3);

    /* N-N pairs */
    REGISTER_AP("AP_NN1", "N-N pairs at distance 1", ap_NN1);
    REGISTER_AP("AP_NN2", "N-N pairs at distance 2", ap_NN2);
    REGISTER_AP("AP_NN3", "N-N pairs at distance 3", ap_NN3);
    REGISTER_AP("AP_NN4", "N-N pairs at distance 4", ap_NN4);
    REGISTER_AP("AP_NN5", "N-N pairs at distance 5", ap_NN5);

    /* N-O pairs */
    REGISTER_AP("AP_NO1", "N-O pairs at distance 1", ap_NO1);
    REGISTER_AP("AP_NO2", "N-O pairs at distance 2", ap_NO2);
    REGISTER_AP("AP_NO3", "N-O pairs at distance 3", ap_NO3);
    REGISTER_AP("AP_NO4", "N-O pairs at distance 4", ap_NO4);

    /* O-O pairs */
    REGISTER_AP("AP_OO1", "O-O pairs at distance 1", ap_OO1);
    REGISTER_AP("AP_OO2", "O-O pairs at distance 2", ap_OO2);
    REGISTER_AP("AP_OO3", "O-O pairs at distance 3", ap_OO3);
    REGISTER_AP("AP_OO4", "O-O pairs at distance 4", ap_OO4);

    /* S-S pairs */
    REGISTER_AP("AP_SS1", "S-S pairs at distance 1", ap_SS1);
    REGISTER_AP("AP_SS2", "S-S pairs at distance 2", ap_SS2);

    /* N-S pairs */
    REGISTER_AP("AP_NS1", "N-S pairs at distance 1", ap_NS1);
    REGISTER_AP("AP_NS2", "N-S pairs at distance 2", ap_NS2);
    REGISTER_AP("AP_NS3", "N-S pairs at distance 3", ap_NS3);

    /* O-S pairs */
    REGISTER_AP("AP_OS1", "O-S pairs at distance 1", ap_OS1);
    REGISTER_AP("AP_OS2", "O-S pairs at distance 2", ap_OS2);

    /* Hal-Hal pairs */
    REGISTER_AP("AP_HalHal1", "Hal-Hal pairs at distance 1", ap_HalHal1);
    REGISTER_AP("AP_HalHal2", "Hal-Hal pairs at distance 2", ap_HalHal2);
    REGISTER_AP("AP_HalHal3", "Hal-Hal pairs at distance 3", ap_HalHal3);

    /* Summary descriptors */
    REGISTER_AP("AP_TotalD1", "Total pairs at distance 1", ap_total_d1);
    REGISTER_AP("AP_TotalD2", "Total pairs at distance 2", ap_total_d2);
    REGISTER_AP("AP_TotalD3", "Total pairs at distance 3", ap_total_d3);
    REGISTER_AP("AP_TotalD4", "Total pairs at distance 4", ap_total_d4);
    REGISTER_AP("AP_TotalD5", "Total pairs at distance 5", ap_total_d5);
    REGISTER_AP("AP_HeteroTotal", "Total heteroatom pairs", ap_hetero_total);
}
