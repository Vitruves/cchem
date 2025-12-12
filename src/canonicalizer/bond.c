/**
 * @file bond.c
 * @brief Bond data structure and operations implementation
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/bond.h"

bond_t* bond_create(int atom1, int atom2, bond_type_t type) {
    bond_t* bond = (bond_t*)calloc(1, sizeof(bond_t));
    if (!bond) return NULL;

    bond_init(bond, atom1, atom2, type);
    return bond;
}

void bond_free(bond_t* bond) {
    if (bond) {
        free(bond);
    }
}

void bond_init(bond_t* bond, int atom1, int atom2, bond_type_t type) {
    if (!bond) return;

    memset(bond, 0, sizeof(bond_t));
    bond->atom1 = atom1;
    bond->atom2 = atom2;
    bond->type = type;
    bond->aromatic = (type == BOND_AROMATIC || type == BOND_RING_AROMATIC);
    bond->in_ring = false;
    bond->stereo = STEREO_NONE;
    bond->stereo_type = BOND_NONE;
    bond->stereo_atom = -1;
    bond->index = -1;
    bond->canon_idx = -1;
    bond->visited = false;
    bond->num_rings = 0;
}

void bond_reset(bond_t* bond) {
    if (bond) {
        int a1 = bond->atom1;
        int a2 = bond->atom2;
        bond_type_t t = bond->type;
        bond_init(bond, a1, a2, t);
    }
}

double bond_get_order(const bond_t* bond) {
    if (!bond) return 0.0;

    switch (bond->type) {
        case BOND_NONE:         return 0.0;
        case BOND_SINGLE:       return 1.0;
        case BOND_DOUBLE:       return 2.0;
        case BOND_TRIPLE:       return 3.0;
        case BOND_AROMATIC:     return 1.5;
        case BOND_UP:           return 1.0;
        case BOND_DOWN:         return 1.0;
        case BOND_RING_SINGLE:  return 1.0;
        case BOND_RING_AROMATIC: return 1.5;
        default:                return 1.0;
    }
}

int bond_get_int_order(const bond_t* bond) {
    if (!bond) return 0;

    switch (bond->type) {
        case BOND_NONE:         return 0;
        case BOND_SINGLE:       return 1;
        case BOND_DOUBLE:       return 2;
        case BOND_TRIPLE:       return 3;
        case BOND_AROMATIC:     return 1;  /* Treat as single for valence */
        case BOND_UP:           return 1;
        case BOND_DOWN:         return 1;
        case BOND_RING_SINGLE:  return 1;
        case BOND_RING_AROMATIC: return 1;
        default:                return 1;
    }
}

int bond_get_other_atom(const bond_t* bond, int atom_idx) {
    if (!bond) return -1;

    if (bond->atom1 == atom_idx) return bond->atom2;
    if (bond->atom2 == atom_idx) return bond->atom1;
    return -1;
}

bool bond_contains_atom(const bond_t* bond, int atom_idx) {
    if (!bond) return false;
    return (bond->atom1 == atom_idx || bond->atom2 == atom_idx);
}

char bond_to_smiles_char(bond_type_t type) {
    switch (type) {
        case BOND_NONE:         return '\0';
        case BOND_SINGLE:       return '-';
        case BOND_DOUBLE:       return '=';
        case BOND_TRIPLE:       return '#';
        case BOND_AROMATIC:     return ':';
        case BOND_UP:           return '/';
        case BOND_DOWN:         return '\\';
        case BOND_RING_SINGLE:  return '-';
        case BOND_RING_AROMATIC: return ':';
        default:                return '\0';
    }
}

bond_type_t bond_from_smiles_char(char c) {
    switch (c) {
        case '-':  return BOND_SINGLE;
        case '=':  return BOND_DOUBLE;
        case '#':  return BOND_TRIPLE;
        case '$':  return BOND_TRIPLE;  /* Quadruple bond, treat as triple */
        case ':':  return BOND_AROMATIC;
        case '/':  return BOND_UP;
        case '\\': return BOND_DOWN;
        default:   return BOND_NONE;
    }
}

bool bond_is_stereo(bond_type_t type) {
    return (type == BOND_UP || type == BOND_DOWN);
}

void bond_copy(bond_t* dest, const bond_t* src) {
    if (!dest || !src) return;
    memcpy(dest, src, sizeof(bond_t));
}
