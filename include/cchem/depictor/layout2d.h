/**
 * @file layout2d.h
 * @brief Internal interfaces for 2D molecular layout system
 *
 * This header defines the internal data structures and function interfaces
 * used by the modular 2D coordinate generation system.
 */

#ifndef CCHEM_DEPICTOR_LAYOUT2D_H
#define CCHEM_DEPICTOR_LAYOUT2D_H

#include "cchem/depictor/types.h"
#include "cchem/depictor/coords2d.h"
#include "cchem/canonicalizer/molecule.h"
#include <stdbool.h>

/* Maximum number of rings in a ring system */
#define MAX_RINGS_PER_SYSTEM 32

/* Maximum shared atoms between rings */
#define MAX_SHARED_ATOMS 16

/* ============== Ring System Classification ============== */

typedef enum {
    RING_SYSTEM_ISOLATED,   /* Single ring */
    RING_SYSTEM_FUSED,      /* Edge-fused (naphthalene-like) */
    RING_SYSTEM_SPIRO,      /* Single atom shared */
    RING_SYSTEM_BRIDGED     /* Multiple paths between atoms */
} ring_system_type_t;

typedef struct {
    int* ring_indices;          /* Indices into molecule's ring array */
    int num_rings;              /* Number of rings in this system */
    ring_system_type_t type;    /* Classification of the system */
    int* shared_atoms;          /* Atoms shared between rings */
    int num_shared_atoms;       /* Number of shared atoms */
    int* all_atoms;             /* All atoms in the ring system */
    int num_atoms;              /* Total atoms in system */
    bool placed;                /* Has this system been placed */
} ring_system_t;

/* ============== Layout Context ============== */

typedef struct {
    const molecule_t* mol;      /* Source molecule (not owned) */
    mol_coords_t* coords;       /* Coordinates being generated (not owned) */
    bool* placed;               /* Placement flags per atom */
    double bond_length;         /* Target bond length */
    int* atom_to_ring_system;   /* Map from atom index to ring system index (-1 if not in ring system) */
    ring_system_t* ring_systems;/* Array of ring systems */
    int num_ring_systems;       /* Number of ring systems */
    int* chain_dir;             /* Direction for zigzag chain placement */
    const coords2d_options_t* options; /* Layout options */
} layout_context_t;

/* ============== Layout Context Management ============== */

/**
 * Create a new layout context
 * @param mol Source molecule
 * @param coords Coordinate storage (must be pre-allocated)
 * @param options Layout options
 * @return New context or NULL on error
 */
layout_context_t* layout_context_create(const molecule_t* mol,
                                        mol_coords_t* coords,
                                        const coords2d_options_t* options);

/**
 * Free a layout context and its resources
 * @param ctx Context to free
 */
void layout_context_free(layout_context_t* ctx);

/* ============== Geometry Helpers ============== */

/**
 * Calculate circumradius for a regular polygon
 * @param num_sides Number of sides
 * @param edge_length Length of each edge
 * @return Circumradius
 */
double layout_ring_circumradius(int num_sides, double edge_length);

/**
 * Calculate the centroid of a set of points
 * @param points Array of points
 * @param num_points Number of points
 * @return Centroid point
 */
point2d_t layout_centroid(const point2d_t* points, int num_points);

/**
 * Calculate centroid of specified atoms
 * @param ctx Layout context
 * @param atom_indices Array of atom indices
 * @param num_atoms Number of atoms
 * @return Centroid point
 */
point2d_t layout_atom_centroid(const layout_context_t* ctx,
                               const int* atom_indices, int num_atoms);

/* ============== Ring System Detection (ring_systems.c) ============== */

/**
 * Detect and classify ring systems in the molecule
 * Groups rings by connectivity and classifies as fused/spiro/bridged
 * @param ctx Layout context (ring_systems will be populated)
 * @return 0 on success, -1 on error
 */
int layout_detect_ring_systems(layout_context_t* ctx);

/**
 * Count shared atoms between two rings
 * @param r1 First ring
 * @param r2 Second ring
 * @param shared1_idx Optional output: indices in r1 of shared atoms
 * @param shared2_idx Optional output: indices in r2 of shared atoms
 * @param max_shared Maximum shared atoms to find
 * @return Number of shared atoms
 */
int layout_count_shared_atoms(const ring_t* r1, const ring_t* r2,
                              int* shared1_idx, int* shared2_idx, int max_shared);

/**
 * Check if two rings are fused (share at least 2 atoms)
 * @param r1 First ring
 * @param r2 Second ring
 * @return true if fused
 */
bool layout_rings_are_fused(const ring_t* r1, const ring_t* r2);

/**
 * Free ring systems in a context
 * @param ctx Layout context
 */
void layout_free_ring_systems(layout_context_t* ctx);

/* ============== Ring Ordering (ring_systems.c) ============== */

/**
 * Get ring atoms in bond-connected order starting from a given atom
 * @param ring Ring to traverse
 * @param mol Molecule containing the ring
 * @param start_atom Atom to start from
 * @param ordered Output array (must be ring->size capacity)
 * @param num_ordered Output: number of atoms in order
 */
void layout_get_ring_order(const ring_t* ring, const molecule_t* mol,
                           int start_atom, int* ordered, int* num_ordered);

/* ============== Ring Placement (ring_placement.c) ============== */

/**
 * Place all ring systems in the molecule
 * Uses templates when available, falls back to algorithmic placement
 * @param ctx Layout context
 * @return 0 on success, -1 on error
 */
int layout_place_ring_systems(layout_context_t* ctx);

/**
 * Place a ring as a regular polygon
 * @param ctx Layout context
 * @param ring Ring to place
 * @param center Center point
 * @param start_angle Starting angle for first atom
 * @param start_atom Atom to place first
 */
void layout_place_ring_polygon(layout_context_t* ctx, const ring_t* ring,
                               point2d_t center, double start_angle, int start_atom);

/**
 * Place a fused ring using already-placed atoms as anchors
 * @param ctx Layout context
 * @param ring Ring to place
 * @return true on success, false if not enough anchors
 */
bool layout_place_fused_ring(layout_context_t* ctx, const ring_t* ring);

/**
 * Place a spiro ring (shares exactly one atom with placed atoms)
 * @param ctx Layout context
 * @param ring Ring to place
 * @param anchor_atom The shared (already placed) atom
 * @return true on success
 */
bool layout_place_spiro_ring(layout_context_t* ctx, const ring_t* ring, int anchor_atom);

/**
 * Place a substituent ring (connected via bond to placed atoms)
 * @param ctx Layout context
 * @param ring Ring to place
 * @param anchor_atom Ring atom connected to placed neighbor
 * @param placed_neighbor The placed atom anchor_atom is bonded to
 * @return true on success
 */
bool layout_place_substituent_ring(layout_context_t* ctx, const ring_t* ring,
                                   int anchor_atom, int placed_neighbor);

/* ============== Chain Placement (chain_placement.c) ============== */

/**
 * Place all chain (non-ring) atoms
 * Uses BFS from placed atoms with zigzag angle calculation
 * @param ctx Layout context
 * @return 0 on success, -1 on error
 */
int layout_place_chains(layout_context_t* ctx);

/**
 * Trace a chain from a starting atom to find if it connects to a ring
 * @param mol Molecule
 * @param start_chain_atom Starting chain atom
 * @param from_ring_atom Atom we came from
 * @param chain_length Output: length of chain (optional)
 * @return Ring atom at other end, or -1 if terminal
 */
int layout_trace_chain_to_ring(const molecule_t* mol, int start_chain_atom,
                               int from_ring_atom, int* chain_length);

/* ============== Collision Detection (collision.c) ============== */

/**
 * Check if a proposed ring placement causes collisions
 * @param ctx Layout context
 * @param center Proposed ring center
 * @param radius Ring circumradius
 * @param ring Ring being placed (atoms to exclude)
 * @param exclude_atom Additional atom to exclude (e.g., anchor neighbor)
 * @return true if collision detected
 */
bool layout_check_ring_collision(const layout_context_t* ctx,
                                 point2d_t center, double radius,
                                 const ring_t* ring, int exclude_atom);

/**
 * Check for atom-atom overlaps in the current layout
 * @param ctx Layout context
 * @param threshold Minimum distance threshold
 * @return Number of overlapping pairs
 */
int layout_count_atom_overlaps(const layout_context_t* ctx, double threshold);

/**
 * Check for bond-bond crossings in the current layout
 * @param ctx Layout context
 * @return Number of crossing pairs
 */
int layout_count_bond_crossings(const layout_context_t* ctx);

/**
 * Calculate overall layout quality score (lower is better)
 * @param ctx Layout context
 * @return Quality score
 */
double layout_quality_score(const layout_context_t* ctx);

/**
 * Find collision-free angle for ring placement
 * @param ctx Layout context
 * @param neighbor_pos Position of neighbor atom
 * @param ring Ring being placed
 * @param base_angle Initial angle to try
 * @param radius Ring circumradius
 * @param exclude_atom Atom to exclude from collision checks (e.g., connection point)
 * @param out_angle Output: collision-free angle
 * @return true if found, false if all angles have collisions
 */
bool layout_find_collision_free_angle(const layout_context_t* ctx,
                                      point2d_t neighbor_pos,
                                      const ring_t* ring,
                                      double base_angle, double radius,
                                      int exclude_atom,
                                      double* out_angle);

/* ============== Force-Directed Refinement (force_directed.c) ============== */

/**
 * Refine layout using force-directed algorithm
 * @param ctx Layout context
 * @param iterations Maximum iterations
 * @return Number of iterations performed
 */
int layout_force_directed_refine(layout_context_t* ctx, int iterations);

#endif /* CCHEM_DEPICTOR_LAYOUT2D_H */
