/**
 * @file mol_soa.h
 * @brief Structure of Arrays (SoA) molecule view for SIMD-friendly descriptors
 *
 * The mol_soa_t provides a cache-friendly layout where each property is
 * stored in a contiguous array. This enables:
 * - Better cache utilization (only load what you need)
 * - SIMD vectorization (contiguous data for sum/max/min operations)
 * - 4-8x speedup for descriptor calculations
 *
 * Usage:
 *   mol_soa_t soa;
 *   mol_soa_init(&soa, mol, arena);
 *   // Use soa.elements[], soa.charges[], etc.
 *   // Arena reset clears the SoA automatically
 */

#ifndef CCHEM_MOL_SOA_H
#define CCHEM_MOL_SOA_H

#include <stdbool.h>
#include <stdint.h>
#include "cchem/arena.h"
#include "cchem/canonicalizer/types.h"

/* Forward declaration */
struct molecule;

/**
 * Structure of Arrays view of a molecule.
 * All arrays are num_atoms in length and allocated from arena.
 */
typedef struct mol_soa {
    int num_atoms;
    int num_heavy_atoms;
    int num_bonds;
    int num_rings;

    /* Per-atom arrays (contiguous, cache-friendly) */
    element_t* elements;         /* Atomic numbers */
    int8_t* charges;             /* Formal charges */
    uint8_t* num_neighbors;      /* Degree (number of bonds) */
    uint8_t* implicit_h;         /* Implicit hydrogen count */
    uint8_t* ring_count;         /* Number of rings containing atom */
    uint8_t* aromatic;           /* Is aromatic (0/1) */
    uint8_t* is_heavy;           /* Is heavy atom (0/1) */

    /* Pre-computed counts for fast lookup */
    int count_C;
    int count_N;
    int count_O;
    int count_S;
    int count_F;
    int count_Cl;
    int count_Br;
    int count_I;
    int count_P;
    int count_H;
    int count_aromatic;
    int count_in_ring;

    /* Bond info for batch descriptors */
    int count_single_bonds;
    int count_double_bonds;
    int count_triple_bonds;
    int count_aromatic_bonds;
    int count_rotatable_bonds;

    /* Total charges and sums */
    int total_charge;
    int total_implicit_h;
} mol_soa_t;

/**
 * Initialize SoA view from molecule using arena allocation.
 * Returns true on success, false on allocation failure.
 *
 * Note: The SoA is valid until arena_reset() is called.
 */
bool mol_soa_init(mol_soa_t* soa, const struct molecule* mol, arena_t* arena);

/**
 * Initialize SoA using thread-local arena.
 * Convenience wrapper around mol_soa_init().
 */
bool mol_soa_init_tl(mol_soa_t* soa, const struct molecule* mol);

#endif /* CCHEM_MOL_SOA_H */
