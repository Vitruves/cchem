/**
 * @file memory.h
 * @brief Memory management utilities for cchem
 *
 * This module provides:
 * - Linear memory arena for zero-allocation batch processing
 * - Structure of Arrays (SoA) molecule view for SIMD-friendly descriptors
 */

#ifndef CCHEM_MEMORY_H
#define CCHEM_MEMORY_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "cchem/canonicalizer/types.h"

/* ============================================================================
 * Arena Allocator
 *
 * The arena allocator provides fast, bump-pointer allocation with O(1) reset.
 * Perfect for batch processing where all allocations can be freed together.
 * ============================================================================ */

/* Default arena size: 1MB per thread */
#define ARENA_DEFAULT_SIZE (1024 * 1024)

/* Arena structure */
typedef struct arena {
    uint8_t* data;       /* Base pointer */
    size_t size;         /* Total size */
    size_t offset;       /* Current allocation offset */
    size_t peak;         /* Peak usage (for monitoring) */
} arena_t;

/* Create arena with specified size */
arena_t* arena_create(size_t size);

/* Create arena with default size */
arena_t* arena_create_default(void);

/* Free arena and all memory */
void arena_free(arena_t* arena);

/* Reset arena to empty (O(1) - just resets offset) */
static inline void arena_reset(arena_t* arena) {
    if (arena) {
        if (arena->offset > arena->peak) {
            arena->peak = arena->offset;
        }
        arena->offset = 0;
    }
}

/* Allocate memory from arena (returns NULL if out of space) */
static inline void* arena_alloc(arena_t* arena, size_t size) {
    if (!arena) return NULL;

    /* Align to 8 bytes */
    size_t aligned_size = (size + 7) & ~(size_t)7;

    if (arena->offset + aligned_size > arena->size) {
        return NULL;  /* Out of arena space */
    }

    void* ptr = arena->data + arena->offset;
    arena->offset += aligned_size;
    return ptr;
}

/* Allocate and zero memory */
static inline void* arena_calloc(arena_t* arena, size_t count, size_t elem_size) {
    size_t total = count * elem_size;
    void* ptr = arena_alloc(arena, total);
    if (ptr) {
        __builtin_memset(ptr, 0, total);
    }
    return ptr;
}

/* Get remaining space in arena */
static inline size_t arena_remaining(const arena_t* arena) {
    if (!arena) return 0;
    return arena->size - arena->offset;
}

/* Get peak usage (for tuning arena size) */
static inline size_t arena_peak_usage(const arena_t* arena) {
    if (!arena) return 0;
    return arena->peak > arena->offset ? arena->peak : arena->offset;
}

/* Thread-local arena access */
arena_t* arena_get_thread_local(void);

/* ============================================================================
 * Structure of Arrays (SoA) Molecule View
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
 * ============================================================================ */

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

#endif /* CCHEM_MEMORY_H */
