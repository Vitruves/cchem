/**
 * @file stereo.h
 * @brief Stereochemistry detection and canonicalization
 *
 * Handles:
 * - Tetrahedral chirality (R/S, @/@@)
 * - Double bond stereochemistry (E/Z, cis/trans)
 * - Extended tetrahedral (allene-like)
 */

#ifndef CCHEM_CANONICALIZER_STEREO_H
#define CCHEM_CANONICALIZER_STEREO_H

#include "types.h"
#include "molecule.h"

/* Stereocenter types */
typedef enum {
    STEREO_CENTER_NONE = 0,
    STEREO_CENTER_TETRAHEDRAL,     /* sp3 carbon, etc */
    STEREO_CENTER_DOUBLE_BOND,     /* C=C, C=N, etc */
    STEREO_CENTER_ALLENE,          /* Extended tetrahedral */
    STEREO_CENTER_SQUARE_PLANAR,   /* sp, tb stereochemistry */
    STEREO_CENTER_TRIGONAL_BIPYR,  /* tb stereochemistry */
    STEREO_CENTER_OCTAHEDRAL       /* oh stereochemistry */
} stereo_center_type_t;

/* Stereocenter information */
typedef struct {
    stereo_center_type_t type;
    int center_atom;              /* Central atom for tetrahedral */
    int atoms[6];                 /* Surrounding atoms (up to 6 for octahedral) */
    int num_atoms;                /* Number of surrounding atoms */
    chirality_t chirality;        /* Current chirality assignment */
    bool is_specified;            /* Is stereochemistry specified */
    bool can_be_stereo;           /* Could this be a stereocenter */
} stereo_center_t;

/* Double bond stereo information */
typedef struct {
    int bond_idx;                 /* Index of double bond */
    int atom1, atom2;             /* Atoms of double bond */
    int ref_atom1, ref_atom2;     /* Reference atoms for E/Z */
    stereo_ez_t stereo;           /* E/Z designation */
    bool is_specified;            /* Is stereochemistry specified */
} stereo_double_bond_t;

/* Stereochemistry context for molecule */
typedef struct {
    stereo_center_t* centers;
    int num_centers;
    int centers_capacity;

    stereo_double_bond_t* double_bonds;
    int num_double_bonds;
    int double_bonds_capacity;
} stereo_info_t;

/* Create stereochemistry info structure */
stereo_info_t* stereo_info_create(void);

/* Free stereochemistry info */
void stereo_info_free(stereo_info_t* info);

/* Detect all stereocenters in molecule */
cchem_status_t stereo_detect_centers(molecule_t* mol, stereo_info_t* info);

/* Detect tetrahedral stereocenters */
cchem_status_t stereo_detect_tetrahedral(molecule_t* mol, stereo_info_t* info);

/* Detect double bond stereocenters */
cchem_status_t stereo_detect_double_bonds(molecule_t* mol, stereo_info_t* info);

/* Assign stereochemistry from SMILES parsing info */
cchem_status_t stereo_assign_from_smiles(molecule_t* mol, stereo_info_t* info);

/* Check if atom could be tetrahedral stereocenter */
bool stereo_is_potential_center(const molecule_t* mol, int atom_idx);

/* Check if bond could be E/Z stereocenter */
bool stereo_is_potential_double_bond(const molecule_t* mol, int bond_idx);

/* Determine chirality from neighbor ordering */
chirality_t stereo_calc_chirality(const molecule_t* mol, int center_atom,
                                  const int* neighbors, int num_neighbors);

/* Determine E/Z from neighbor positions */
stereo_ez_t stereo_calc_ez(const molecule_t* mol, int bond_idx,
                           int ref1, int ref2);

/* Canonicalize stereochemistry representation */
cchem_status_t stereo_canonicalize(molecule_t* mol, stereo_info_t* info,
                                   const int* canon_order);

/* Invert chirality at atom */
void stereo_invert_chirality(atom_t* atom);

/* Convert chirality based on new neighbor ordering */
chirality_t stereo_permute_chirality(chirality_t orig, const int* old_order,
                                     const int* new_order, int n);

/* Get canonical chirality symbol for output
 * from_atom is the atom we came from in DFS traversal (-1 if starting atom) */
chirality_t stereo_get_canonical_chirality(const molecule_t* mol, int atom_idx,
                                           int from_atom, const int* output_order);

/* Check stereochemistry equivalence */
bool stereo_centers_equivalent(const stereo_center_t* a, const stereo_center_t* b);

#endif /* CCHEM_CANONICALIZER_STEREO_H */
