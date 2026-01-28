/**
 * @file atom.h
 * @brief Atom data structure and operations
 */

#ifndef CCHEM_CANONICALIZER_ATOM_H
#define CCHEM_CANONICALIZER_ATOM_H

#include "types.h"
#include "element.h"

/* Forward declaration */
struct molecule;

/* Atom structure */
typedef struct atom {
    /* Basic properties */
    element_t element;           /* Atomic number */
    int isotope;                 /* Isotope mass (0 = natural) */
    int charge;                  /* Formal charge */
    int h_count;                 /* Explicit hydrogen count (-1 = implicit) */
    bool aromatic;               /* Is aromatic */
    chirality_t chirality;       /* Stereochemistry */
    int chirality_class;         /* For extended chirality (OH, TB) */

    /* Graph connectivity */
    int index;                   /* Index in molecule */
    int neighbors[MAX_NEIGHBORS]; /* Indices of neighboring atoms */
    int neighbor_bonds[MAX_NEIGHBORS]; /* Indices of bonds to neighbors */
    int num_neighbors;           /* Number of neighbors */

    /* Computed properties */
    int ring_count;              /* Number of rings containing this atom */
    int implicit_h_count;        /* Computed implicit hydrogens */
    int total_bond_order;        /* Sum of bond orders */
    int pi_electrons;            /* Pi electrons contributed to aromaticity (1=pyridine-type, 2=pyrrole-type) */

    /* Canonicalization */
    uint64_t invariant;          /* Canonical invariant */
    int canon_rank;              /* Canonical rank */
    bool visited;                /* For traversal algorithms */
    int atom_class;              /* Atom class [atom:n] */

    /* Original input order for stereo */
    int input_order;             /* Position in original SMILES */
    int stereo_neighbors[MAX_NEIGHBORS]; /* Neighbors in original order */
    int num_stereo_neighbors;

    /* Ring opening tracking for correct stereo neighbor order.
     * When a ring opens at a chiral atom, we record where it should go
     * in stereo_neighbors, then insert the closing atom there later. */
    int ring_open_slots[10];     /* Position in stereo_neighbors for each ring opening */
    int ring_open_nums[10];      /* Ring number for each slot */
    int num_ring_opens;          /* Number of pending ring openings */
} atom_t;

/* Create a new atom */
atom_t* atom_create(element_t element);

/* Free atom memory */
void atom_free(atom_t* atom);

/* Initialize atom with defaults */
void atom_init(atom_t* atom, element_t element);

/* Reset atom to default state */
void atom_reset(atom_t* atom);

/* Add neighbor to atom */
cchem_status_t atom_add_neighbor(atom_t* atom, int neighbor_idx, int bond_idx);

/* Remove neighbor from atom */
cchem_status_t atom_remove_neighbor(atom_t* atom, int neighbor_idx);

/* Get bond index to specific neighbor */
int atom_get_bond_to(const atom_t* atom, int neighbor_idx);

/* Check if atoms are bonded */
bool atom_is_bonded_to(const atom_t* atom, int other_idx);

/* Calculate implicit hydrogen count */
int atom_calc_implicit_h(const atom_t* atom, const struct molecule* mol);

/* Get total degree (including implicit H) */
int atom_get_degree(const atom_t* atom);

/* Get heavy atom degree (excluding H) */
int atom_get_heavy_degree(const atom_t* atom);

/* Check if atom is in ring */
bool atom_in_ring(const atom_t* atom);

/* Copy atom data */
void atom_copy(atom_t* dest, const atom_t* src);

#endif /* CCHEM_CANONICALIZER_ATOM_H */
