/**
 * @file ring_finder.h
 * @brief Ring detection algorithms (SSSR - Smallest Set of Smallest Rings)
 */

#ifndef CCHEM_CANONICALIZER_RING_FINDER_H
#define CCHEM_CANONICALIZER_RING_FINDER_H

#include "types.h"
#include "molecule.h"

/* Ring detection using BFS-based SSSR algorithm */
cchem_status_t ring_finder_sssr(molecule_t* mol);

/* Find all rings up to a maximum size */
cchem_status_t ring_finder_all(molecule_t* mol, int max_ring_size);

/* Get ring membership for atom */
int ring_finder_atom_ring_count(const molecule_t* mol, int atom_idx);

/* Get ring membership for bond */
int ring_finder_bond_ring_count(const molecule_t* mol, int bond_idx);

/* Classify ring as aromatic based on Huckel's rule */
bool ring_finder_is_aromatic_ring(const molecule_t* mol, const ring_t* ring);

/* Detect aromatic systems and mark atoms/bonds */
cchem_status_t ring_finder_perceive_aromaticity(molecule_t* mol);

/* Count pi electrons in ring for aromaticity */
int ring_finder_count_pi_electrons(const molecule_t* mol, const ring_t* ring);

#endif /* CCHEM_CANONICALIZER_RING_FINDER_H */
