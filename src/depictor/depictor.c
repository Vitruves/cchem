/**
 * @file depictor.c
 * @brief Main molecular depiction implementation
 */

#include "cchem/depictor/depictor.h"
#include "cchem/depictor/colors.h"
#include "cchem/canonicalizer/parser.h"
#include "cchem/cchem.h"
#include <stdlib.h>
#include <string.h>

static const depictor_options_t DEFAULT_OPTIONS = DEPICTOR_OPTIONS_DEFAULT;

/* Kekulize aromatic bonds using backtracking algorithm */

typedef struct {
    int* arom_bond_idx;     /* Indices of aromatic bonds */
    int num_arom_bonds;
    bool* has_double;       /* has_double[atom] = true if atom already has a double bond */
    bool* bond_is_double;   /* Result: which aromatic bonds become double */
    molecule_t* mol;
} kekule_state_t;

/* Check if aromatic atom contributes lone pair instead of needing double bond */
static bool atom_has_lone_pair_contribution(const atom_t* atom, const molecule_t* mol) {
    if (!atom->aromatic) return false;

    /* Pyrrole-type nitrogen: aromatic N with 2 neighbors that has implicit H */
    if (atom->element == ELEM_N) {
        /* Check for explicit H first */
        if (atom->h_count > 0) return true;

        /* Aromatic N with exactly 2 neighbors:
         * - Pyrrole-type: contributes lone pair, needs implicit H -> no double bond
         * - Pyridine-type: participates in C=N double bond, no H
         *
         * Distinguish by ring size: 5-ring = pyrrole, 6-ring = pyridine
         * BUT in complex fused systems, SSSR might not find the 5-ring.
         *
         * Alternative: count aromatic bonds. If N has 2 aromatic bonds,
         * check if it's in a 6-ring. If not found in 6-ring, assume pyrrole-type.
         */
        if (atom->num_neighbors == 2 && atom->ring_count > 0) {
            /* First check if in a 5-membered ring */
            bool in_5_ring = false;
            bool in_6_ring = false;

            for (int r = 0; r < mol->num_rings; r++) {
                const ring_t* ring = &mol->rings[r];
                for (int i = 0; i < ring->size; i++) {
                    if (ring->atoms[i] == atom->index) {
                        if (ring->size == 5) in_5_ring = true;
                        if (ring->size == 6) in_6_ring = true;
                        break;
                    }
                }
            }

            /* Pyrrole-type if in 5-ring, or if NOT in any 6-ring
             * (handles cases where 5-ring isn't in SSSR but we know
             * it's not pyridine because not in a 6-ring) */
            if (in_5_ring || !in_6_ring) {
                return true;
            }
        }
    }

    /* Furan-type oxygen: aromatic O with 2 neighbors contributes lone pair */
    if (atom->element == ELEM_O && atom->num_neighbors == 2) {
        return true;
    }

    /* Thiophene-type sulfur: aromatic S with 2 neighbors contributes lone pair */
    if (atom->element == ELEM_S && atom->num_neighbors == 2) {
        return true;
    }

    return false;
}

/* Check if atom needs a double bond (aromatic atom without one yet) */
static bool atom_needs_double(kekule_state_t* state, int atom, int up_to_bond) {
    if (state->has_double[atom]) return false;

    /* Check if atom is in aromatic system */
    bool in_aromatic = false;
    for (int i = 0; i < state->num_arom_bonds; i++) {
        int b = state->arom_bond_idx[i];
        if (state->mol->bonds[b].atom1 == atom || state->mol->bonds[b].atom2 == atom) {
            in_aromatic = true;
            break;
        }
    }
    if (!in_aromatic) return false;

    /* Check if already has double from assigned bonds */
    for (int i = 0; i < up_to_bond; i++) {
        if (!state->bond_is_double[i]) continue;
        int b = state->arom_bond_idx[i];
        if (state->mol->bonds[b].atom1 == atom || state->mol->bonds[b].atom2 == atom) {
            return false;
        }
    }

    /* Atoms with lone pair contribution don't need a double bond */
    atom_t* a = &state->mol->atoms[atom];
    if (atom_has_lone_pair_contribution(a, state->mol)) {
        return false;
    }

    return true;
}

/* Recursive backtracking to find valid Kekule structure */
static bool kekule_solve(kekule_state_t* state, int bond_idx) {
    if (bond_idx >= state->num_arom_bonds) {
        /* Check all aromatic atoms are satisfied */
        for (int i = 0; i < state->mol->num_atoms; i++) {
            if (atom_needs_double(state, i, state->num_arom_bonds)) {
                return false;
            }
        }
        return true;
    }

    int b = state->arom_bond_idx[bond_idx];
    int a1 = state->mol->bonds[b].atom1;
    int a2 = state->mol->bonds[b].atom2;

    /* If this bond was pre-assigned by constraint propagation, skip it */
    if (state->bond_is_double[bond_idx]) {
        return kekule_solve(state, bond_idx + 1);
    }

    /* Check if either atom already has a double bond */
    bool a1_has = state->has_double[a1];
    bool a2_has = state->has_double[a2];

    /* Atoms with lone pair contribution (pyrrole-type N, furan O, thiophene S)
     * should be treated as already having their "double" - they contribute
     * electrons via lone pair, not via double bond */
    if (atom_has_lone_pair_contribution(&state->mol->atoms[a1], state->mol)) {
        a1_has = true;
    }
    if (atom_has_lone_pair_contribution(&state->mol->atoms[a2], state->mol)) {
        a2_has = true;
    }

    /* Also check bonds we've assigned in this recursion */
    for (int i = 0; i < bond_idx; i++) {
        if (!state->bond_is_double[i]) continue;
        int bi = state->arom_bond_idx[i];
        if (state->mol->bonds[bi].atom1 == a1 || state->mol->bonds[bi].atom2 == a1) a1_has = true;
        if (state->mol->bonds[bi].atom1 == a2 || state->mol->bonds[bi].atom2 == a2) a2_has = true;
    }

    /* If both atoms already have doubles, this must be single */
    if (a1_has && a2_has) {
        state->bond_is_double[bond_idx] = false;
        return kekule_solve(state, bond_idx + 1);
    }

    /* If one atom has a double, this must be single */
    if (a1_has || a2_has) {
        state->bond_is_double[bond_idx] = false;
        return kekule_solve(state, bond_idx + 1);
    }

    /* Check if either atom MUST get their double from this bond */
    bool a1_needs = !atom_has_lone_pair_contribution(&state->mol->atoms[a1], state->mol);
    bool a2_needs = !atom_has_lone_pair_contribution(&state->mol->atoms[a2], state->mol);

    /* Count remaining options for each atom */
    int a1_remaining = 0, a2_remaining = 0;
    for (int i = bond_idx; i < state->num_arom_bonds; i++) {
        if (state->bond_is_double[i]) continue;  /* Already assigned */
        int bi = state->arom_bond_idx[i];
        int bi_a1 = state->mol->bonds[bi].atom1;
        int bi_a2 = state->mol->bonds[bi].atom2;

        if (bi_a1 == a1 || bi_a2 == a1) {
            int other = (bi_a1 == a1) ? bi_a2 : bi_a1;
            if (!state->has_double[other]) a1_remaining++;
        }
        if (bi_a1 == a2 || bi_a2 == a2) {
            int other = (bi_a1 == a2) ? bi_a2 : bi_a1;
            if (!state->has_double[other]) a2_remaining++;
        }
    }

    /* If an atom needs a double and this is their only chance, force double */
    if ((a1_needs && a1_remaining == 1) || (a2_needs && a2_remaining == 1)) {
        state->bond_is_double[bond_idx] = true;
        return kekule_solve(state, bond_idx + 1);
    }

    /* Try double first, then single */
    state->bond_is_double[bond_idx] = true;
    if (kekule_solve(state, bond_idx + 1)) {
        return true;
    }

    state->bond_is_double[bond_idx] = false;
    return kekule_solve(state, bond_idx + 1);
}

static void kekulize_molecule(molecule_t* mol) {
    /* Collect aromatic bonds and atoms */
    kekule_state_t state;
    state.mol = mol;
    state.arom_bond_idx = malloc(mol->num_bonds * sizeof(int));
    state.num_arom_bonds = 0;

    /* Track which atoms are in aromatic system */
    bool* in_aromatic = calloc(mol->num_atoms, sizeof(bool));

    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->bonds[b].type == BOND_AROMATIC) {
            state.arom_bond_idx[state.num_arom_bonds++] = b;
            in_aromatic[mol->bonds[b].atom1] = true;
            in_aromatic[mol->bonds[b].atom2] = true;
        }
    }

    if (state.num_arom_bonds == 0) {
        free(state.arom_bond_idx);
        free(in_aromatic);
        return;
    }

    /* Mark atoms that have double bonds from NON-aromatic bonds
     * Only mark if the double bond is entirely outside the aromatic system */
    state.has_double = calloc(mol->num_atoms, sizeof(bool));
    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->bonds[b].type == BOND_DOUBLE) {
            int a1 = mol->bonds[b].atom1;
            int a2 = mol->bonds[b].atom2;
            /* Only mark as having double if BOTH atoms are outside aromatic system
             * OR if the non-aromatic double is to an atom outside the system */
            if (!in_aromatic[a1]) {
                state.has_double[a1] = true;
            }
            if (!in_aromatic[a2]) {
                state.has_double[a2] = true;
            }
            /* If one atom is aromatic and one isn't, the aromatic one gets
             * its double bond need satisfied by this exocyclic double bond */
            if (in_aromatic[a1] && !in_aromatic[a2]) {
                state.has_double[a1] = true;
            }
            if (in_aromatic[a2] && !in_aromatic[a1]) {
                state.has_double[a2] = true;
            }
        }
    }

    free(in_aromatic);

    /* Initialize bond assignment */
    state.bond_is_double = calloc(state.num_arom_bonds, sizeof(bool));

    /* Solve using backtracking - constraint propagation is built into the solver */
    kekule_solve(&state, 0);

    /* Apply the solution */
    for (int i = 0; i < state.num_arom_bonds; i++) {
        int b = state.arom_bond_idx[i];
        mol->bonds[b].type = state.bond_is_double[i] ? BOND_DOUBLE : BOND_SINGLE;
    }

    free(state.bond_is_double);
    free(state.has_double);
    free(state.arom_bond_idx);
}

cchem_status_t depict_molecule(const molecule_t* mol, const char* filename,
                               const depictor_options_t* options) {
    if (!mol || !filename) return CCHEM_ERROR_INVALID_INPUT;

    depictor_options_t opts = options ? *options : DEFAULT_OPTIONS;

    /* Clone molecule for kekulization */
    molecule_t* work_mol = molecule_clone(mol);
    if (!work_mol) return CCHEM_ERROR_MEMORY;

    /* Ensure rings are computed */
    if (!work_mol->rings_computed) {
        molecule_find_rings(work_mol);
        molecule_perceive_aromaticity(work_mol);
    }

    /* Kekulize if not drawing aromatic circles */
    if (!opts.draw_aromatic_circles) {
        kekulize_molecule(work_mol);
    }

    /* Generate coordinates */
    mol_coords_t* coords = NULL;

    if (opts.mode == DEPICT_MODE_3D) {
        coords3d_options_t c3d_opts = COORDS3D_OPTIONS_DEFAULT;
        c3d_opts.max_iterations = opts.max_iterations;
        coords = coords3d_generate(work_mol, &c3d_opts);

        if (coords && coords->has_3d) {
            /* Project 3D to 2D (simple orthographic projection) */
            coords->has_2d = true;
            for (int i = 0; i < coords->num_atoms; i++) {
                coords->coords_2d[i].x = coords->coords_3d[i].x;
                coords->coords_2d[i].y = coords->coords_3d[i].y;
            }
        }
    } else {
        coords2d_options_t c2d_opts = COORDS2D_OPTIONS_DEFAULT;
        c2d_opts.bond_length = opts.bond_length / 25.0;
        coords = coords2d_generate(work_mol, &c2d_opts);
    }

    if (!coords) {
        molecule_free(work_mol);
        return CCHEM_ERROR_MEMORY;
    }

    /* Apply scale factor for higher resolution output */
    double scale = (opts.scale_factor > 0.0) ? opts.scale_factor : 1.0;
    int scaled_width = (int)(opts.width * scale);
    int scaled_height = (int)(opts.height * scale);
    int scaled_margin = (int)(opts.margin * scale);

    /* Scale coordinates to fit image */
    coords2d_scale_to_fit(coords, scaled_width, scaled_height, scaled_margin);

    /* Create scaled options for rendering */
    depictor_options_t scaled_opts = opts;
    scaled_opts.bond_width *= scale;
    scaled_opts.font_size *= scale;
    scaled_opts.bond_length *= scale;

    /* Render - use appropriate context for format */
    render_context_t* ctx = NULL;
    if (opts.format == IMG_FORMAT_SVG) {
        ctx = render_context_create_ex(scaled_width, scaled_height, opts.background,
                                        IMG_FORMAT_SVG, filename);
    } else {
        ctx = render_context_create_ex(scaled_width, scaled_height, opts.background,
                                        IMG_FORMAT_PNG, NULL);
    }

    if (!ctx) {
        mol_coords_free(coords);
        molecule_free(work_mol);
        return CCHEM_ERROR_MEMORY;
    }

    cchem_status_t status = render_molecule(ctx, work_mol, coords, &scaled_opts);

    if (status == CCHEM_OK) {
        if (opts.format == IMG_FORMAT_SVG) {
            status = render_save_svg(ctx, filename);
        } else if (opts.format == IMG_FORMAT_PNG) {
            status = render_save_png(ctx, filename);
        } else {
            status = render_save_jpeg(ctx, filename, opts.jpeg_quality);
        }
    }

    render_context_free(ctx);
    mol_coords_free(coords);
    molecule_free(work_mol);

    return status;
}

cchem_status_t depict_smiles(const char* smiles, const char* filename,
                             const depictor_options_t* options,
                             char* error_buf, size_t error_buf_size) {
    return depict_smiles_verbose(smiles, filename, options, NULL, error_buf, error_buf_size);
}

cchem_status_t depict_smiles_verbose(const char* smiles, const char* filename,
                                     const depictor_options_t* options,
                                     depict_info_t* info,
                                     char* error_buf, size_t error_buf_size) {
    if (!smiles || !filename) return CCHEM_ERROR_INVALID_INPUT;

    /* Canonicalize input SMILES first for consistent depiction */
    char* canonical = cchem_canonicalize(smiles, error_buf, error_buf_size);
    if (!canonical) {
        return CCHEM_ERROR_PARSE;
    }

    /* Store canonical SMILES in info if requested */
    if (info) {
        strncpy(info->canonical_smiles, canonical, sizeof(info->canonical_smiles) - 1);
        info->canonical_smiles[sizeof(info->canonical_smiles) - 1] = '\0';
        info->energy_initial = 0.0;
        info->energy_final = 0.0;
    }

    molecule_t* mol = smiles_to_molecule(canonical, error_buf, error_buf_size);
    free(canonical);

    if (!mol) {
        return CCHEM_ERROR_PARSE;
    }

    molecule_calc_implicit_h(mol);
    molecule_find_rings(mol);
    molecule_perceive_aromaticity(mol);

    /* Store molecule info */
    if (info) {
        info->num_atoms = mol->num_atoms;
        info->num_bonds = mol->num_bonds;
        info->num_rings = mol->num_rings;
    }

    /* Depict with energy reporting for 3D mode */
    depictor_options_t opts = options ? *options : DEFAULT_OPTIONS;

    /* Clone molecule for kekulization */
    molecule_t* work_mol = molecule_clone(mol);
    if (!work_mol) {
        molecule_free(mol);
        return CCHEM_ERROR_MEMORY;
    }

    /* Ensure rings are computed */
    if (!work_mol->rings_computed) {
        molecule_find_rings(work_mol);
        molecule_perceive_aromaticity(work_mol);
    }

    /* Kekulize if not drawing aromatic circles */
    if (!opts.draw_aromatic_circles) {
        kekulize_molecule(work_mol);
    }

    /* Generate coordinates */
    mol_coords_t* coords = NULL;

    if (opts.mode == DEPICT_MODE_3D) {
        coords3d_options_t c3d_opts = COORDS3D_OPTIONS_DEFAULT;
        c3d_opts.max_iterations = opts.max_iterations;

        if (info) {
            coords = coords3d_generate_with_energy(work_mol, &c3d_opts,
                                                   &info->energy_initial,
                                                   &info->energy_final);
        } else {
            coords = coords3d_generate(work_mol, &c3d_opts);
        }

        if (coords && coords->has_3d) {
            coords->has_2d = true;
            for (int i = 0; i < coords->num_atoms; i++) {
                coords->coords_2d[i].x = coords->coords_3d[i].x;
                coords->coords_2d[i].y = coords->coords_3d[i].y;
            }
        }
    } else {
        coords2d_options_t c2d_opts = COORDS2D_OPTIONS_DEFAULT;
        c2d_opts.bond_length = opts.bond_length / 25.0;
        coords = coords2d_generate(work_mol, &c2d_opts);
    }

    if (!coords) {
        molecule_free(work_mol);
        molecule_free(mol);
        return CCHEM_ERROR_MEMORY;
    }

    /* Calculate effective margin based on render style */
    int effective_margin = opts.margin;
    if (opts.render_style == RENDER_STYLE_SPACEFILL) {
        /* Spacefill needs much larger margin for VDW spheres */
        effective_margin = opts.margin + opts.width / 6;
    } else if (opts.render_style == RENDER_STYLE_BALLS_AND_STICKS) {
        /* Balls-and-sticks needs moderate extra margin */
        effective_margin = opts.margin + opts.width / 12;
    } else if (opts.render_style == RENDER_STYLE_SURFACE) {
        /* Surface also uses VDW radii */
        effective_margin = opts.margin + opts.width / 6;
    }

    /* Apply scale factor for higher resolution output */
    double scale = (opts.scale_factor > 0.0) ? opts.scale_factor : 1.0;
    int scaled_width = (int)(opts.width * scale);
    int scaled_height = (int)(opts.height * scale);
    int scaled_margin = (int)(effective_margin * scale);

    /* Scale coordinates to fit image */
    coords2d_scale_to_fit(coords, scaled_width, scaled_height, scaled_margin);

    /* Create scaled options for rendering */
    depictor_options_t scaled_opts = opts;
    scaled_opts.bond_width *= scale;
    scaled_opts.font_size *= scale;
    scaled_opts.bond_length *= scale;

    /* Render - use appropriate context for format */
    render_context_t* ctx = NULL;
    if (opts.format == IMG_FORMAT_SVG) {
        ctx = render_context_create_ex(scaled_width, scaled_height, opts.background,
                                        IMG_FORMAT_SVG, filename);
    } else {
        ctx = render_context_create_ex(scaled_width, scaled_height, opts.background,
                                        IMG_FORMAT_PNG, NULL);
    }

    if (!ctx) {
        mol_coords_free(coords);
        molecule_free(work_mol);
        molecule_free(mol);
        return CCHEM_ERROR_MEMORY;
    }

    cchem_status_t status = render_molecule(ctx, work_mol, coords, &scaled_opts);

    if (status == CCHEM_OK) {
        if (opts.format == IMG_FORMAT_SVG) {
            status = render_save_svg(ctx, filename);
        } else if (opts.format == IMG_FORMAT_PNG) {
            status = render_save_png(ctx, filename);
        } else {
            status = render_save_jpeg(ctx, filename, opts.jpeg_quality);
        }
    }

    render_context_free(ctx);
    mol_coords_free(coords);
    molecule_free(work_mol);
    molecule_free(mol);

    return status;
}
