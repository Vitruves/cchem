/**
 * @file molecule.h
 * @brief Molecular graph data structure and operations
 */

#ifndef CCHEM_CANONICALIZER_MOLECULE_H
#define CCHEM_CANONICALIZER_MOLECULE_H

#include "types.h"
#include "atom.h"
#include "bond.h"

/* Ring structure */
typedef struct ring {
    int atoms[MAX_RING_SIZE];    /* Atom indices in ring */
    int bonds[MAX_RING_SIZE];    /* Bond indices in ring */
    int size;                    /* Ring size */
    bool aromatic;               /* Is aromatic ring */
    int index;                   /* Ring index */
} ring_t;

/* Molecule structure */
typedef struct molecule {
    /* Graph components */
    atom_t* atoms;               /* Array of atoms */
    int num_atoms;               /* Number of atoms */
    int atoms_capacity;          /* Allocated capacity for atoms */

    bond_t* bonds;               /* Array of bonds */
    int num_bonds;               /* Number of bonds */
    int bonds_capacity;          /* Allocated capacity for bonds */

    /* Ring information */
    ring_t* rings;               /* Array of rings (SSSR) */
    int num_rings;               /* Number of rings */
    int rings_capacity;          /* Allocated capacity for rings */
    bool rings_computed;         /* Have rings been computed */

    /* Connectivity */
    int num_fragments;           /* Number of disconnected fragments */
    int* fragment_ids;           /* Fragment ID for each atom */

    /* Properties */
    int total_charge;            /* Sum of formal charges */
    double molecular_weight;     /* Molecular weight */

    /* Canonicalization state */
    int* canon_order;            /* Canonical atom ordering */
    bool is_canonical;           /* Has molecule been canonicalized */

    /* Original SMILES (for reference) */
    char* original_smiles;       /* Input SMILES string */

    /* Error handling */
    char error_msg[256];         /* Last error message */
} molecule_t;

/* Create a new empty molecule */
molecule_t* molecule_create(void);

/* Create molecule with initial capacity */
molecule_t* molecule_create_with_capacity(int atom_capacity, int bond_capacity);

/* Free molecule memory */
void molecule_free(molecule_t* mol);

/* Reset molecule to empty state */
void molecule_reset(molecule_t* mol);

/* Add an atom to the molecule, returns atom index or -1 on error */
int molecule_add_atom(molecule_t* mol, element_t element);

/* Add a pre-configured atom, returns atom index or -1 on error */
int molecule_add_atom_full(molecule_t* mol, const atom_t* atom);

/* Add a bond between atoms, returns bond index or -1 on error */
int molecule_add_bond(molecule_t* mol, int atom1, int atom2, bond_type_t type);

/* Get atom by index */
atom_t* molecule_get_atom(molecule_t* mol, int idx);
const atom_t* molecule_get_atom_const(const molecule_t* mol, int idx);

/* Get bond by index */
bond_t* molecule_get_bond(molecule_t* mol, int idx);
const bond_t* molecule_get_bond_const(const molecule_t* mol, int idx);

/* Get bond between two atoms, returns NULL if not found */
bond_t* molecule_get_bond_between(molecule_t* mol, int atom1, int atom2);
const bond_t* molecule_get_bond_between_const(const molecule_t* mol, int atom1, int atom2);

/* Get bond index between two atoms, returns -1 if not found */
int molecule_get_bond_index(const molecule_t* mol, int atom1, int atom2);

/* Remove atom (and its bonds), returns status */
cchem_status_t molecule_remove_atom(molecule_t* mol, int idx);

/* Remove bond, returns status */
cchem_status_t molecule_remove_bond(molecule_t* mol, int idx);

/* Calculate implicit hydrogens for all atoms */
void molecule_calc_implicit_h(molecule_t* mol);

/* Find rings using SSSR algorithm */
cchem_status_t molecule_find_rings(molecule_t* mol);

/* Detect aromaticity */
cchem_status_t molecule_perceive_aromaticity(molecule_t* mol);

/* Identify disconnected fragments */
cchem_status_t molecule_find_fragments(molecule_t* mol);

/* Calculate molecular weight */
double molecule_calc_weight(const molecule_t* mol);

/* Calculate total charge */
int molecule_calc_charge(const molecule_t* mol);

/* Validate molecule structure */
cchem_status_t molecule_validate(const molecule_t* mol);

/* Check if molecule is connected (single fragment) */
bool molecule_is_connected(const molecule_t* mol);

/* Get number of heavy atoms */
int molecule_num_heavy_atoms(const molecule_t* mol);

/* Get number of rotatable bonds */
int molecule_num_rotatable_bonds(const molecule_t* mol);

/* Clone a molecule */
molecule_t* molecule_clone(const molecule_t* mol);

/* Get last error message */
const char* molecule_get_error(const molecule_t* mol);

#endif /* CCHEM_CANONICALIZER_MOLECULE_H */
