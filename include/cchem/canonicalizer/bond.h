/**
 * @file bond.h
 * @brief Bond data structure and operations
 */

#ifndef CCHEM_CANONICALIZER_BOND_H
#define CCHEM_CANONICALIZER_BOND_H

#include "types.h"

/* Forward declaration */
struct molecule;

/* Bond structure */
typedef struct bond {
    int atom1;              /* First atom index */
    int atom2;              /* Second atom index */
    bond_type_t type;       /* Bond type */
    bool aromatic;          /* Is aromatic bond */
    bool in_ring;           /* Is part of a ring */
    stereo_ez_t stereo;     /* E/Z stereochemistry */

    /* For stereo representation */
    bond_type_t stereo_type; /* Original stereo bond type (UP/DOWN) */
    int stereo_atom;         /* Reference atom for stereo direction */

    /* Canonicalization */
    int index;              /* Index in molecule */
    int canon_idx;          /* Canonical index */
    bool visited;           /* For traversal */

    /* Ring membership */
    int ring_ids[MAX_RINGS]; /* IDs of rings containing this bond */
    int num_rings;           /* Number of rings */
} bond_t;

/* Create a new bond */
bond_t* bond_create(int atom1, int atom2, bond_type_t type);

/* Free bond memory */
void bond_free(bond_t* bond);

/* Initialize bond with defaults */
void bond_init(bond_t* bond, int atom1, int atom2, bond_type_t type);

/* Reset bond to default state */
void bond_reset(bond_t* bond);

/* Get bond order (1, 2, 3, or 1.5 for aromatic) */
double bond_get_order(const bond_t* bond);

/* Get integer bond order */
int bond_get_int_order(const bond_t* bond);

/* Get the other atom in the bond */
int bond_get_other_atom(const bond_t* bond, int atom_idx);

/* Check if bond contains atom */
bool bond_contains_atom(const bond_t* bond, int atom_idx);

/* Get bond type character for SMILES output */
char bond_to_smiles_char(bond_type_t type);

/* Parse bond type from SMILES character */
bond_type_t bond_from_smiles_char(char c);

/* Check if bond type is stereo (UP/DOWN) */
bool bond_is_stereo(bond_type_t type);

/* Copy bond data */
void bond_copy(bond_t* dest, const bond_t* src);

#endif /* CCHEM_CANONICALIZER_BOND_H */
