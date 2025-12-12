/**
 * @file invariant.h
 * @brief Atom invariant calculation for canonicalization
 *
 * Implements the initial atom invariants and refinement steps
 * based on the Morgan/Weininger algorithm.
 */

#ifndef CCHEM_CANONICALIZER_INVARIANT_H
#define CCHEM_CANONICALIZER_INVARIANT_H

#include "types.h"
#include "molecule.h"

/* Initial invariant components (bit-packed) */
typedef struct {
    uint8_t num_connections;    /* Number of heavy atom connections */
    uint8_t num_h;              /* Number of attached hydrogens */
    int8_t charge;              /* Formal charge */
    uint8_t atomic_num;         /* Atomic number */
    uint8_t mass;               /* Isotope mass (0 = natural) */
    uint8_t num_rings;          /* Number of rings atom is in */
    bool aromatic;              /* Is aromatic */
} invariant_components_t;

/* Calculate initial invariants for all atoms */
cchem_status_t invariant_calc_initial(molecule_t* mol);

/* Calculate initial invariant for single atom */
uint64_t invariant_calc_atom(const molecule_t* mol, int atom_idx);

/* Refine invariants iteratively (Morgan algorithm) */
cchem_status_t invariant_refine(molecule_t* mol, int max_iterations);

/* Single refinement step */
int invariant_refine_step(molecule_t* mol);

/* Convert invariants to ranks */
cchem_status_t invariant_to_ranks(molecule_t* mol);

/* Get number of distinct ranks */
int invariant_count_distinct(const molecule_t* mol);

/* Break ties in ranks using various strategies */
cchem_status_t invariant_break_ties(molecule_t* mol);

/* Check if all atoms have unique ranks */
bool invariant_all_unique(const molecule_t* mol);

/* Pack invariant components into 64-bit integer */
uint64_t invariant_pack(const invariant_components_t* comp);

/* Unpack 64-bit invariant to components */
void invariant_unpack(uint64_t inv, invariant_components_t* comp);

#endif /* CCHEM_CANONICALIZER_INVARIANT_H */
