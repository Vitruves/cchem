/**
 * @file depictor.h
 * @brief Main molecular depiction API
 */

#ifndef CCHEM_DEPICTOR_DEPICTOR_H
#define CCHEM_DEPICTOR_DEPICTOR_H

#include "cchem/depictor/types.h"
#include "cchem/depictor/coords2d.h"
#include "cchem/depictor/coords3d.h"
#include "cchem/depictor/render.h"
#include "cchem/canonicalizer/molecule.h"

/* Depiction result info (for verbose output) */
typedef struct {
    char canonical_smiles[512];   /* Canonical SMILES used */
    double energy_initial;        /* Initial MMFF94 energy (3D only) */
    double energy_final;          /* Final MMFF94 energy (3D only) */
    int num_atoms;                /* Number of atoms */
    int num_bonds;                /* Number of bonds */
    int num_rings;                /* Number of rings */
} depict_info_t;

cchem_status_t depict_molecule(const molecule_t* mol, const char* filename,
                               const depictor_options_t* options);
cchem_status_t depict_smiles(const char* smiles, const char* filename,
                             const depictor_options_t* options,
                             char* error_buf, size_t error_buf_size);

/* Verbose version that returns depiction info */
cchem_status_t depict_smiles_verbose(const char* smiles, const char* filename,
                                     const depictor_options_t* options,
                                     depict_info_t* info,
                                     char* error_buf, size_t error_buf_size);

#endif /* CCHEM_DEPICTOR_DEPICTOR_H */
