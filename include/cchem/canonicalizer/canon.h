/**
 * @file canon.h
 * @brief Main canonicalization algorithm
 *
 * Implements a complete canonicalization algorithm based on:
 * - Weininger's original SMILES canonicalization
 * - Morgan algorithm for topological equivalence
 * - Tie-breaking strategies for unique ordering
 */

#ifndef CCHEM_CANONICALIZER_CANON_H
#define CCHEM_CANONICALIZER_CANON_H

#include "types.h"
#include "molecule.h"
#include "stereo.h"
#include "parser.h"

/* Canonicalization options */
typedef struct {
    bool preserve_stereo;        /* Keep stereochemistry */
    bool kekulize;               /* Convert aromatic to Kekule form */
    bool remove_explicit_h;      /* Remove explicit H atoms */
    bool remove_atom_maps;       /* Remove atom class numbers */
    bool canonical_stereo;       /* Use canonical stereo representation */
    bool use_isomeric;           /* Include stereochemistry in output */
    bool use_isotopes;           /* Include isotopes in output */
    bool use_charges;            /* Include charges in output */
} canon_options_t;

/* Default canonicalization options */
extern const canon_options_t CANON_OPTIONS_DEFAULT;

/* Strict canonical options (all features) */
extern const canon_options_t CANON_OPTIONS_STRICT;

/* Minimal options (no stereo, no isotopes) */
extern const canon_options_t CANON_OPTIONS_MINIMAL;

/* Canonicalization context */
typedef struct {
    molecule_t* mol;             /* Molecule being canonicalized */
    stereo_info_t* stereo;       /* Stereochemistry info */
    canon_options_t options;     /* Canonicalization options */

    /* Working data */
    int* ranks;                  /* Current atom ranks */
    uint64_t* invariants;        /* Atom invariants */
    int* order;                  /* Canonical ordering */
    int num_atoms;               /* Number of atoms */

    /* State tracking */
    bool computed;               /* Has canonicalization been done */
    int iterations;              /* Number of refinement iterations */

    /* Error handling */
    char error_msg[256];
} canon_context_t;

/* Create canonicalization context */
canon_context_t* canon_context_create(molecule_t* mol, const canon_options_t* options);

/* Free canonicalization context */
void canon_context_free(canon_context_t* ctx);

/* Main canonicalization function */
cchem_status_t canon_canonicalize(canon_context_t* ctx);

/* Step 1: Calculate initial atom invariants */
cchem_status_t canon_calc_invariants(canon_context_t* ctx);

/* Step 2: Iteratively refine invariants */
cchem_status_t canon_refine_invariants(canon_context_t* ctx);

/* Step 3: Break ties to get unique ordering */
cchem_status_t canon_break_ties(canon_context_t* ctx);

/* Step 4: Generate canonical ordering */
cchem_status_t canon_generate_order(canon_context_t* ctx);

/* Step 5: Canonicalize stereochemistry */
cchem_status_t canon_stereo(canon_context_t* ctx);

/* Get canonical order array */
const int* canon_get_order(const canon_context_t* ctx);

/* Get canonical rank for atom */
int canon_get_rank(const canon_context_t* ctx, int atom_idx);

/* High-level: canonicalize molecule in place */
cchem_status_t molecule_canonicalize(molecule_t* mol, const canon_options_t* options);

/* High-level: get canonical SMILES from molecule */
char* molecule_to_canonical_smiles(const molecule_t* mol, const canon_options_t* options);

/* Highest level: canonical SMILES from SMILES */
char* smiles_canonicalize(const char* smiles, const canon_options_t* options,
                          char* error_buf, size_t error_buf_size);

/* Check if two SMILES represent the same molecule */
bool smiles_are_equivalent(const char* smiles1, const char* smiles2);

#endif /* CCHEM_CANONICALIZER_CANON_H */
