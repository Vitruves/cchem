/**
 * @file coords3d.h
 * @brief 3D coordinate generation for molecular depiction
 */

#ifndef CCHEM_DEPICTOR_COORDS3D_H
#define CCHEM_DEPICTOR_COORDS3D_H

#include "cchem/depictor/types.h"
#include "cchem/canonicalizer/molecule.h"

typedef struct {
    int random_seed;
    bool enforce_planarity;
    bool use_experimental_torsions;
    int max_iterations;
    double energy_tolerance;
    double gradient_tolerance;
    bool use_mmff94;
} coords3d_options_t;

#define COORDS3D_OPTIONS_DEFAULT { \
    .random_seed = 42, \
    .enforce_planarity = true, \
    .use_experimental_torsions = true, \
    .max_iterations = 500, \
    .energy_tolerance = 1e-4, \
    .gradient_tolerance = 1e-3, \
    .use_mmff94 = true \
}

mol_coords_t* coords3d_generate(const molecule_t* mol, const coords3d_options_t* options);
int coords3d_count_clashes(const mol_coords_t* coords, const molecule_t* mol, double threshold);
double coords3d_minimize(mol_coords_t* coords, const molecule_t* mol, const coords3d_options_t* options);

/* Generate with energy reporting */
mol_coords_t* coords3d_generate_with_energy(const molecule_t* mol, const coords3d_options_t* options,
                                            double* energy_initial, double* energy_final);

#endif /* CCHEM_DEPICTOR_COORDS3D_H */
