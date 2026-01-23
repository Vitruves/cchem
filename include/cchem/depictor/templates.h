/**
 * @file templates.h
 * @brief Ring template definitions for 2D molecular layout
 *
 * Pre-defined templates for common ring systems with optimal 2D coordinates.
 * Templates are matched against ring systems and transformed to fit.
 */

#ifndef CCHEM_DEPICTOR_TEMPLATES_H
#define CCHEM_DEPICTOR_TEMPLATES_H

#include "cchem/depictor/types.h"
#include "cchem/depictor/layout2d.h"
#include "cchem/canonicalizer/molecule.h"
#include <stdbool.h>

/* Maximum atoms in a template */
#define MAX_TEMPLATE_ATOMS 32

/* Maximum rings in a template */
#define MAX_TEMPLATE_RINGS 8

/* ============== Template Structure ============== */

typedef struct {
    const char* name;               /* Template name (e.g., "benzene", "naphthalene") */
    int num_atoms;                  /* Number of atoms in template */
    int ring_sizes[MAX_TEMPLATE_RINGS]; /* Sizes of rings in template */
    int num_rings;                  /* Number of rings */
    point2d_t coords[MAX_TEMPLATE_ATOMS]; /* Template coordinates (normalized) */
    int adjacency[MAX_TEMPLATE_ATOMS][4]; /* Adjacency list per atom (max 4 neighbors, -1 terminated) */
} ring_template_t;

/* ============== Template Matching ============== */

/**
 * Find a template matching the given ring system
 * @param ctx Layout context
 * @param system Ring system to match
 * @return Matching template or NULL if no match
 */
const ring_template_t* template_find_match(const layout_context_t* ctx,
                                           const ring_system_t* system);

/**
 * Apply a template to a ring system
 * Places atoms according to template coordinates with proper transformation
 * @param ctx Layout context
 * @param system Ring system to place
 * @param templ Template to apply
 * @param center Center point for placement
 * @param scale Scale factor (typically bond_length)
 * @param rotation Rotation angle in radians
 * @return true on success
 */
bool template_apply(layout_context_t* ctx, const ring_system_t* system,
                    const ring_template_t* templ, point2d_t center,
                    double scale, double rotation);

/**
 * Calculate the best transformation to apply a template
 * Finds optimal rotation to align template with any already-placed atoms
 * @param ctx Layout context
 * @param system Ring system
 * @param templ Template to apply
 * @param out_center Output: optimal center
 * @param out_rotation Output: optimal rotation
 * @return true if transformation found
 */
bool template_find_transformation(const layout_context_t* ctx,
                                  const ring_system_t* system,
                                  const ring_template_t* templ,
                                  point2d_t* out_center, double* out_rotation);

/* ============== Built-in Templates ============== */

/**
 * Get the array of built-in ring templates
 * @param count Output: number of templates
 * @return Array of templates
 */
const ring_template_t* template_get_builtin(int* count);

/* ============== Template Helper Functions ============== */

/**
 * Check if a ring system matches a template's ring signature
 * @param system Ring system
 * @param templ Template to check
 * @param mol Molecule (for ring access)
 * @return true if ring sizes match
 */
bool template_ring_signature_matches(const ring_system_t* system,
                                     const ring_template_t* templ,
                                     const molecule_t* mol);

/**
 * Build atom mapping between ring system and template
 * @param ctx Layout context
 * @param system Ring system
 * @param templ Template
 * @param mapping Output: mapping[system_atom] = template_atom (-1 if unmapped)
 * @return true if mapping found
 */
bool template_build_atom_mapping(const layout_context_t* ctx,
                                 const ring_system_t* system,
                                 const ring_template_t* templ,
                                 int* mapping);

#endif /* CCHEM_DEPICTOR_TEMPLATES_H */
