/**
 * @file force_directed.h
 * @brief Force-directed layout refinement for 2D molecular depiction
 *
 * Implements a force-directed algorithm to refine molecular layouts
 * by optimizing bond lengths, angles, and preventing overlaps.
 */

#ifndef CCHEM_DEPICTOR_FORCE_DIRECTED_H
#define CCHEM_DEPICTOR_FORCE_DIRECTED_H

#include "cchem/depictor/types.h"
#include "cchem/depictor/layout2d.h"
#include <stdbool.h>

/* ============== Force Model Parameters ============== */

typedef struct {
    /* Spring force (bond length) */
    double spring_strength;     /* Spring constant for bonds (default: 1.0) */
    double ideal_bond_length;   /* Target bond length */

    /* Repulsion force (atom overlap prevention) */
    double repulsion_strength;  /* Repulsion constant (default: 1.0) */
    double min_atom_distance;   /* Minimum allowed atom distance */

    /* Angular force (ideal bond angles) */
    double angular_strength;    /* Angular spring constant (default: 0.5) */

    /* Planarity force (ring flatness) */
    double planarity_strength;  /* Force to keep rings planar (default: 0.3) */

    /* Integration parameters */
    double damping;             /* Velocity damping (default: 0.85) */
    double max_displacement;    /* Maximum displacement per step */
    double convergence_threshold; /* Energy change threshold for convergence */

    /* Constraints */
    bool fix_ring_atoms;        /* Don't move atoms in aromatic rings */
    bool fix_placed_anchors;    /* Don't move atoms used as anchors */
} force_params_t;

/* Default force parameters */
#define FORCE_PARAMS_DEFAULT { \
    .spring_strength = 1.0, \
    .ideal_bond_length = 1.5, \
    .repulsion_strength = 1.0, \
    .min_atom_distance = 0.5, \
    .angular_strength = 0.5, \
    .planarity_strength = 0.3, \
    .damping = 0.85, \
    .max_displacement = 0.5, \
    .convergence_threshold = 0.001, \
    .fix_ring_atoms = false, \
    .fix_placed_anchors = false \
}

/* ============== Force Calculation ============== */

/**
 * Calculate total energy of the current layout
 * @param ctx Layout context
 * @param params Force parameters
 * @return Total energy (lower is better)
 */
double force_calculate_energy(const layout_context_t* ctx,
                              const force_params_t* params);

/**
 * Calculate forces on all atoms
 * @param ctx Layout context
 * @param params Force parameters
 * @param forces Output: force vectors for each atom
 */
void force_calculate_forces(const layout_context_t* ctx,
                            const force_params_t* params,
                            point2d_t* forces);

/**
 * Calculate spring force between bonded atoms
 * @param p1 Position of first atom
 * @param p2 Position of second atom
 * @param ideal_length Target bond length
 * @param strength Spring constant
 * @return Force vector on first atom
 */
point2d_t force_spring(point2d_t p1, point2d_t p2,
                       double ideal_length, double strength);

/**
 * Calculate repulsion force between atoms
 * @param p1 Position of first atom
 * @param p2 Position of second atom
 * @param min_distance Minimum allowed distance
 * @param strength Repulsion constant
 * @return Force vector on first atom
 */
point2d_t force_repulsion(point2d_t p1, point2d_t p2,
                          double min_distance, double strength);

/**
 * Calculate angular force at a central atom
 * @param center Position of central atom
 * @param neighbor1 Position of first neighbor
 * @param neighbor2 Position of second neighbor
 * @param ideal_angle Target angle in radians
 * @param strength Angular spring constant
 * @param out_force1 Output: force on neighbor1
 * @param out_force2 Output: force on neighbor2
 */
void force_angular(point2d_t center, point2d_t neighbor1, point2d_t neighbor2,
                   double ideal_angle, double strength,
                   point2d_t* out_force1, point2d_t* out_force2);

/* ============== Integration ============== */

/**
 * Perform one step of force-directed refinement
 * @param ctx Layout context
 * @param params Force parameters
 * @param velocities Current velocities (modified in place)
 * @return Energy after step
 */
double force_step(layout_context_t* ctx, const force_params_t* params,
                  point2d_t* velocities);

/**
 * Run force-directed refinement until convergence or max iterations
 * @param ctx Layout context
 * @param params Force parameters (or NULL for defaults)
 * @param max_iterations Maximum iterations
 * @return Number of iterations performed
 */
int force_refine(layout_context_t* ctx, const force_params_t* params,
                 int max_iterations);

/* ============== Ideal Angle Calculation ============== */

/**
 * Get ideal bond angle for an atom based on hybridization
 * @param mol Molecule
 * @param atom_idx Atom index
 * @return Ideal angle in radians (e.g., 2π/3 for sp2)
 */
double force_ideal_angle(const molecule_t* mol, int atom_idx);

/**
 * Get ideal angles around an atom based on neighbor count
 * @param num_neighbors Number of neighbors
 * @param angles Output: ideal angles between consecutive neighbors
 */
void force_ideal_angles_for_neighbors(int num_neighbors, double* angles);

#endif /* CCHEM_DEPICTOR_FORCE_DIRECTED_H */
