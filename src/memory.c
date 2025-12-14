/**
 * @file memory.c
 * @brief Memory management utilities implementation
 */

#include <stdlib.h>
#include <string.h>
#include "cchem/memory.h"
#include "cchem/canonicalizer/molecule.h"

/* ============================================================================
 * Arena Allocator Implementation
 * ============================================================================ */

arena_t* arena_create(size_t size) {
    arena_t* arena = (arena_t*)malloc(sizeof(arena_t));
    if (!arena) return NULL;

    arena->data = (uint8_t*)malloc(size);
    if (!arena->data) {
        free(arena);
        return NULL;
    }

    arena->size = size;
    arena->offset = 0;
    arena->peak = 0;

    return arena;
}

arena_t* arena_create_default(void) {
    return arena_create(ARENA_DEFAULT_SIZE);
}

void arena_free(arena_t* arena) {
    if (arena) {
        if (arena->data) free(arena->data);
        free(arena);
    }
}

/* Thread-local arena - one per thread */
static __thread arena_t* tl_arena = NULL;

arena_t* arena_get_thread_local(void) {
    if (!tl_arena) {
        tl_arena = arena_create_default();
    }
    return tl_arena;
}

/* ============================================================================
 * Structure of Arrays (SoA) Molecule View Implementation
 * ============================================================================ */

bool mol_soa_init(mol_soa_t* soa, const molecule_t* mol, arena_t* arena) {
    if (!soa || !mol || !arena) return false;

    int n = mol->num_atoms;
    soa->num_atoms = n;
    soa->num_bonds = mol->num_bonds;
    soa->num_rings = mol->num_rings;

    if (n == 0) {
        /* Empty molecule */
        soa->elements = NULL;
        soa->charges = NULL;
        soa->num_neighbors = NULL;
        soa->implicit_h = NULL;
        soa->ring_count = NULL;
        soa->aromatic = NULL;
        soa->is_heavy = NULL;
        soa->num_heavy_atoms = 0;
        memset(&soa->count_C, 0, sizeof(int) * 20);  /* Zero all counts */
        return true;
    }

    /* Allocate all arrays from arena (single allocation block, aligned) */
    size_t elem_size = n * sizeof(element_t);
    size_t byte_size = n * sizeof(uint8_t);
    size_t i8_size = n * sizeof(int8_t);

    soa->elements = (element_t*)arena_alloc(arena, elem_size);
    soa->charges = (int8_t*)arena_alloc(arena, i8_size);
    soa->num_neighbors = (uint8_t*)arena_alloc(arena, byte_size);
    soa->implicit_h = (uint8_t*)arena_alloc(arena, byte_size);
    soa->ring_count = (uint8_t*)arena_alloc(arena, byte_size);
    soa->aromatic = (uint8_t*)arena_alloc(arena, byte_size);
    soa->is_heavy = (uint8_t*)arena_alloc(arena, byte_size);

    if (!soa->elements || !soa->charges || !soa->num_neighbors ||
        !soa->implicit_h || !soa->ring_count || !soa->aromatic || !soa->is_heavy) {
        return false;  /* Arena out of memory */
    }

    /* Initialize counters */
    soa->count_C = 0;
    soa->count_N = 0;
    soa->count_O = 0;
    soa->count_S = 0;
    soa->count_F = 0;
    soa->count_Cl = 0;
    soa->count_Br = 0;
    soa->count_I = 0;
    soa->count_P = 0;
    soa->count_H = 0;
    soa->count_aromatic = 0;
    soa->count_in_ring = 0;
    soa->num_heavy_atoms = 0;
    soa->total_charge = 0;
    soa->total_implicit_h = 0;

    /* Copy data from AoS to SoA with counting in single pass */
    for (int i = 0; i < n; i++) {
        const atom_t* atom = &mol->atoms[i];

        element_t elem = atom->element;
        soa->elements[i] = elem;
        soa->charges[i] = (int8_t)atom->charge;
        soa->num_neighbors[i] = (uint8_t)atom->num_neighbors;
        soa->implicit_h[i] = (uint8_t)(atom->implicit_h_count > 0 ? atom->implicit_h_count : 0);
        soa->ring_count[i] = (uint8_t)atom->ring_count;
        soa->aromatic[i] = atom->aromatic ? 1 : 0;

        bool is_heavy = (elem != ELEM_H);
        soa->is_heavy[i] = is_heavy ? 1 : 0;

        /* Update counters */
        soa->total_charge += atom->charge;
        soa->total_implicit_h += soa->implicit_h[i];

        if (is_heavy) {
            soa->num_heavy_atoms++;
        }

        if (atom->aromatic) {
            soa->count_aromatic++;
        }

        if (atom->ring_count > 0) {
            soa->count_in_ring++;
        }

        /* Element counting */
        switch (elem) {
            case ELEM_C:  soa->count_C++; break;
            case ELEM_N:  soa->count_N++; break;
            case ELEM_O:  soa->count_O++; break;
            case ELEM_S:  soa->count_S++; break;
            case ELEM_F:  soa->count_F++; break;
            case ELEM_Cl: soa->count_Cl++; break;
            case ELEM_Br: soa->count_Br++; break;
            case ELEM_I:  soa->count_I++; break;
            case ELEM_P:  soa->count_P++; break;
            case ELEM_H:  soa->count_H++; break;
            default: break;
        }
    }

    /* Count bond types */
    soa->count_single_bonds = 0;
    soa->count_double_bonds = 0;
    soa->count_triple_bonds = 0;
    soa->count_aromatic_bonds = 0;
    soa->count_rotatable_bonds = 0;

    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = &mol->bonds[i];

        switch (bond->type) {
            case BOND_SINGLE:   soa->count_single_bonds++; break;
            case BOND_DOUBLE:   soa->count_double_bonds++; break;
            case BOND_TRIPLE:   soa->count_triple_bonds++; break;
            case BOND_AROMATIC: soa->count_aromatic_bonds++; break;
            default: break;
        }

        /* Count rotatable bonds (single, non-ring, between heavy atoms) */
        if (bond->type == BOND_SINGLE && !bond->in_ring) {
            int a1 = bond->atom1;
            int a2 = bond->atom2;
            if (mol->atoms[a1].element != ELEM_H &&
                mol->atoms[a2].element != ELEM_H &&
                mol->atoms[a1].num_neighbors > 1 &&
                mol->atoms[a2].num_neighbors > 1) {
                soa->count_rotatable_bonds++;
            }
        }
    }

    return true;
}

bool mol_soa_init_tl(mol_soa_t* soa, const molecule_t* mol) {
    return mol_soa_init(soa, mol, arena_get_thread_local());
}
