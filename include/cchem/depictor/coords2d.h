/**
 * @file coords2d.h
 * @brief 2D coordinate generation for molecular depiction
 */

#ifndef CCHEM_DEPICTOR_COORDS2D_H
#define CCHEM_DEPICTOR_COORDS2D_H

#include "cchem/depictor/types.h"
#include "cchem/canonicalizer/molecule.h"

typedef struct {
    double bond_length;
    bool use_templates;
    int max_iterations;
    double clash_threshold;
    bool debug;
} coords2d_options_t;

#define COORDS2D_OPTIONS_DEFAULT { \
    .bond_length = 1.5, \
    .use_templates = true, \
    .max_iterations = 500, \
    .clash_threshold = 0.7, \
    .debug = false \
}

mol_coords_t* coords2d_generate(const molecule_t* mol, const coords2d_options_t* options);
int coords2d_refine(mol_coords_t* coords, const molecule_t* mol, const coords2d_options_t* options);
void coords2d_scale_to_fit(mol_coords_t* coords, double width, double height, double margin);
void coords2d_center(mol_coords_t* coords);

#endif /* CCHEM_DEPICTOR_COORDS2D_H */
