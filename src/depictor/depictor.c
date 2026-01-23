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
    if (!atom->aromatic) {
        if (atom->element == ELEM_N) {
            fprintf(stderr, "DEBUG kekule: N atom %d NOT aromatic, skipping lone pair check\n", atom->index);
        }
        return false;
    }

    /* Pyrrole-type nitrogen: aromatic N that contributes lone pair */
    if (atom->element == ELEM_N) {
        /* Check for explicit H first */
        if (atom->h_count > 0) return true;

        /* Check what rings this N is in */
        bool in_5_ring = false;

        for (int r = 0; r < mol->num_rings; r++) {
            const ring_t* ring = &mol->rings[r];
            for (int i = 0; i < ring->size; i++) {
                if (ring->atoms[i] == atom->index) {
                    if (ring->size == 5) in_5_ring = true;
                    break;
                }
            }
        }

        /* Aromatic N with 3 neighbors in a 5-membered ring:
         * This is like pyrrole N with a substituent (e.g., N-phenyl triazole).
         * The N contributes its lone pair to aromaticity and should NOT
         * have a double bond (already has 3 single bonds). */
        if (atom->num_neighbors == 3 && in_5_ring) {
            return true;
        }

        /* For N with 2 neighbors in a 5-membered ring:
         * In rings like imidazole/benzimidazole with 2 N's (both with 2 neighbors),
         * exactly one contributes the lone pair. Use heuristic: if there's another
         * N with 3 neighbors (substituted), this one doesn't contribute. Otherwise,
         * if there are multiple N's with 2 neighbors, the highest-indexed one
         * contributes (arbitrary but consistent choice). */
        if (atom->num_neighbors == 2 && in_5_ring) {
            for (int r = 0; r < mol->num_rings; r++) {
                const ring_t* ring = &mol->rings[r];
                if (ring->size != 5) continue;

                /* Check if this atom is in this ring */
                bool atom_in_ring = false;
                for (int i = 0; i < ring->size; i++) {
                    if (ring->atoms[i] == atom->index) {
                        atom_in_ring = true;
                        break;
                    }
                }
                if (!atom_in_ring) continue;

                /* Check other N's in this ring */
                bool has_3_neighbor_n = false;
                int max_2_neighbor_n = -1;
                for (int i = 0; i < ring->size; i++) {
                    const atom_t* ra = &mol->atoms[ring->atoms[i]];
                    if (ra->element != ELEM_N) continue;
                    if (ra->num_neighbors == 3) {
                        has_3_neighbor_n = true;
                    } else if (ra->num_neighbors == 2) {
                        if (ring->atoms[i] > max_2_neighbor_n) {
                            max_2_neighbor_n = ring->atoms[i];
                        }
                    }
                }

                /* If another N has 3 neighbors, it's the lone pair contributor */
                if (has_3_neighbor_n) {
                    return false;  /* This N participates in double bonds */
                }

                /* Among 2-neighbor N's, the highest-indexed one contributes lone pair */
                if (atom->index == max_2_neighbor_n) {
                    return true;  /* This N contributes lone pair */
                }
            }
        }
    }

    /* Furan-type oxygen: aromatic O with 2 neighbors contributes lone pair */
    if (atom->element == ELEM_O && atom->num_neighbors == 2) {
        return true;
    }

    /* Aromatic sulfur contributes lone pair (thiophene, sulfone, etc.)
     * S in aromatic systems contributes via lone pair regardless of
     * additional substituents like =O in sulfones */
    if (atom->element == ELEM_S) {
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
        /* Also verify no atom has more than one double (from aromatic bonds) */
        for (int a = 0; a < state->mol->num_atoms; a++) {
            int count = 0;
            for (int i = 0; i < state->num_arom_bonds; i++) {
                if (!state->bond_is_double[i]) continue;
                int bi = state->arom_bond_idx[i];
                if (state->mol->bonds[bi].atom1 == a || state->mol->bonds[bi].atom2 == a) {
                    count++;
                }
            }
            if (count > 1) {
                return false;  /* Invalid: atom has multiple aromatic doubles */
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

    /* Check if either atom already has a double bond (from explicit DOUBLE or assigned) */
    bool a1_has_double = state->has_double[a1];
    bool a2_has_double = state->has_double[a2];

    /* Also check bonds we've assigned in this recursion */
    for (int i = 0; i < bond_idx; i++) {
        if (!state->bond_is_double[i]) continue;
        int bi = state->arom_bond_idx[i];
        if (state->mol->bonds[bi].atom1 == a1 || state->mol->bonds[bi].atom2 == a1) a1_has_double = true;
        if (state->mol->bonds[bi].atom1 == a2 || state->mol->bonds[bi].atom2 == a2) a2_has_double = true;
    }

    /* Check if either atom contributes lone pair (cannot accept double bond) */
    bool a1_lone_pair = atom_has_lone_pair_contribution(&state->mol->atoms[a1], state->mol);
    bool a2_lone_pair = atom_has_lone_pair_contribution(&state->mol->atoms[a2], state->mol);

    /* Treat lone pair as having a double (satisfied, can't take another) */
    if (a1_lone_pair) a1_has_double = true;
    if (a2_lone_pair) a2_has_double = true;

    static int debug_depth = 0;
    if (debug_depth < 3) {
        fprintf(stderr, "DEBUG kekule[%d]: bond %d-%d, a1_has_dbl=%d, a2_has_dbl=%d, a1_lp=%d, a2_lp=%d\n",
                bond_idx, a1, a2, a1_has_double, a2_has_double, a1_lone_pair, a2_lone_pair);
    }

    /* If both atoms already have doubles, this must be single */
    if (a1_has_double && a2_has_double) {
        state->bond_is_double[bond_idx] = false;
        return kekule_solve(state, bond_idx + 1);
    }

    /* If one atom already has a double, this must be single
     * (to avoid giving that atom multiple doubles) */
    if (a1_has_double || a2_has_double) {
        state->bond_is_double[bond_idx] = false;
        return kekule_solve(state, bond_idx + 1);
    }

    /* Neither atom has a double yet - try both options */
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

    /* Mark atoms that already have double bonds (from explicit BOND_DOUBLE) */
    state.has_double = calloc(mol->num_atoms, sizeof(bool));
    for (int b = 0; b < mol->num_bonds; b++) {
        if (mol->bonds[b].type == BOND_DOUBLE) {
            int a1 = mol->bonds[b].atom1;
            int a2 = mol->bonds[b].atom2;
            /* Mark BOTH atoms as having their double bond need satisfied */
            state.has_double[a1] = true;
            state.has_double[a2] = true;
        }
    }

    free(in_aromatic);

    /* Initialize bond assignment */
    state.bond_is_double = calloc(state.num_arom_bonds, sizeof(bool));

    /* Solve using backtracking - constraint propagation is built into the solver */
    bool success = kekule_solve(&state, 0);

    fprintf(stderr, "DEBUG kekule: solve returned %s\n", success ? "SUCCESS" : "FAILED");
    if (success) {
        /* Apply the solution */
        for (int i = 0; i < state.num_arom_bonds; i++) {
            int b = state.arom_bond_idx[i];
            mol->bonds[b].type = state.bond_is_double[i] ? BOND_DOUBLE : BOND_SINGLE;
            if (state.bond_is_double[i]) {
                fprintf(stderr, "DEBUG kekule: bond %d (%d-%d) assigned DOUBLE\n",
                        b, mol->bonds[b].atom1, mol->bonds[b].atom2);
            }
        }
    } else {
        /* Failed to find valid Kekule structure - convert all aromatic to single */
        for (int i = 0; i < state.num_arom_bonds; i++) {
            int b = state.arom_bond_idx[i];
            mol->bonds[b].type = BOND_SINGLE;
        }
    }

    free(state.bond_is_double);
    free(state.has_double);
    free(state.arom_bond_idx);
}

/* Generate wedge bonds for tetrahedral chiral centers based on 2D coordinates */
static void assign_stereo_bonds(molecule_t* mol, const mol_coords_t* coords) {
    if (!mol || !coords || !coords->has_2d) return;

    for (int i = 0; i < mol->num_atoms; i++) {
        atom_t* atom = &mol->atoms[i];
        if (atom->chirality == CHIRALITY_NONE) continue;

        /* Need at least 3 neighbors for tetrahedral chirality
         * (4th could be implicit H) */
        if (atom->num_neighbors < 3) continue;

        point2d_t center_pos = coords->coords_2d[i];

        /* Calculate average position of all neighbors (for determining "outward") */
        point2d_t avg_neighbor = {0, 0};
        for (int j = 0; j < atom->num_neighbors; j++) {
            int nb = atom->neighbors[j];
            avg_neighbor.x += coords->coords_2d[nb].x;
            avg_neighbor.y += coords->coords_2d[nb].y;
        }
        avg_neighbor.x /= atom->num_neighbors;
        avg_neighbor.y /= atom->num_neighbors;

        /* Score all bonds and find best two for wedge/dash */
        typedef struct {
            int bond_idx;
            int neighbor;
            double score;
        } bond_score_t;

        bond_score_t candidates[MAX_NEIGHBORS];
        int num_candidates = 0;

        for (int j = 0; j < atom->num_neighbors; j++) {
            int nb = atom->neighbors[j];
            int bond_idx = atom->neighbor_bonds[j];
            if (bond_idx < 0) continue;

            bond_t* bond = &mol->bonds[bond_idx];

            /* Skip if bond already has stereo marking */
            if (bond->stereo_type != BOND_NONE) continue;

            /* Skip double/triple bonds */
            if (bond->type == BOND_DOUBLE || bond->type == BOND_TRIPLE) continue;

            double score = 0.0;

            /* Prefer non-ring bonds (more visible) */
            if (!bond->in_ring) score += 100.0;

            /* Prefer bonds pointing away from molecule center */
            point2d_t nb_pos = coords->coords_2d[nb];
            point2d_t to_nb = {nb_pos.x - center_pos.x, nb_pos.y - center_pos.y};
            point2d_t to_avg = {avg_neighbor.x - center_pos.x, avg_neighbor.y - center_pos.y};

            /* Dot product: positive if pointing same direction as average */
            double dot = to_nb.x * to_avg.x + to_nb.y * to_avg.y;
            score -= dot;  /* Prefer opposite direction (outward) */

            /* Prefer bonds to carbon (methyl groups are common stereo indicators) */
            if (mol->atoms[nb].element == ELEM_C) score += 10.0;

            candidates[num_candidates].bond_idx = bond_idx;
            candidates[num_candidates].neighbor = nb;
            candidates[num_candidates].score = score;
            num_candidates++;
        }

        /* Sort candidates by score (descending) - simple bubble sort for small array */
        for (int a = 0; a < num_candidates - 1; a++) {
            for (int b = a + 1; b < num_candidates; b++) {
                if (candidates[b].score > candidates[a].score) {
                    bond_score_t tmp = candidates[a];
                    candidates[a] = candidates[b];
                    candidates[b] = tmp;
                }
            }
        }

        /* Assign wedge (up) to best bond, dash (down) to second best
         * This shows both "front" and "back" substituents */
        if (num_candidates >= 1) {
            bond_t* bond1 = &mol->bonds[candidates[0].bond_idx];
            /* CW (@) = first substituent goes up (solid wedge)
             * CCW (@@) = first substituent goes down (dashed wedge) */
            bool first_is_up = (atom->chirality == CHIRALITY_CW);
            bond1->stereo_type = first_is_up ? BOND_UP : BOND_DOWN;
            bond1->stereo_atom = candidates[0].neighbor;
        }

        if (num_candidates >= 2) {
            bond_t* bond2 = &mol->bonds[candidates[1].bond_idx];
            /* Second bond gets opposite stereo */
            bool first_is_up = (atom->chirality == CHIRALITY_CW);
            bond2->stereo_type = first_is_up ? BOND_DOWN : BOND_UP;
            bond2->stereo_atom = candidates[1].neighbor;
        }
    }
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

    /* Assign wedge bonds for chiral centers */
    assign_stereo_bonds(work_mol, coords);

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
    char* canonical = smiles_canonicalize(smiles, NULL, error_buf, error_buf_size);
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

    /* Assign wedge bonds for chiral centers */
    assign_stereo_bonds(work_mol, coords);

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
