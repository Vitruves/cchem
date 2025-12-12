/**
 * @file parser.c
 * @brief SMILES parser implementation - builds molecular graph from tokens
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/parser.h"

cchem_status_t parser_init(parser_t* parser) {
    if (!parser) return CCHEM_ERROR_INVALID_INPUT;

    memset(parser, 0, sizeof(parser_t));

    for (int i = 0; i < 100; i++) {
        parser->ring_open[i] = false;
        parser->ring_closures[i].atom_idx = -1;
        parser->ring_closures[i].bond_type = BOND_NONE;
        parser->ring_closures[i].has_bond = false;
    }

    parser->branch_depth = 0;
    parser->current_atom = -1;
    parser->pending_bond = BOND_NONE;
    parser->has_pending_bond = false;
    parser->num_stereo_bonds = 0;
    parser->has_error = false;
    parser->error_position = -1;
    parser->error_msg[0] = '\0';
    parser->mol = NULL;

    return CCHEM_OK;
}

void parser_free(parser_t* parser) {
    /* Parser doesn't own the molecule, just clear state */
    if (parser) {
        memset(parser, 0, sizeof(parser_t));
    }
}

static void parser_set_error(parser_t* parser, const char* msg, int position) {
    parser->has_error = true;
    parser->error_position = position;
    snprintf(parser->error_msg, sizeof(parser->error_msg), "%s", msg);
}

/* Process atom token and add to molecule */
static cchem_status_t parser_process_atom(parser_t* parser, const token_t* token) {
    atom_t atom;
    atom_init(&atom, token->data.atom.element);

    atom.isotope = token->data.atom.isotope;
    atom.charge = token->data.atom.charge;
    atom.h_count = token->data.atom.h_count;
    atom.aromatic = token->data.atom.aromatic;
    atom.chirality = token->data.atom.chirality;
    atom.chirality_class = token->data.atom.chirality_class;
    atom.atom_class = token->data.atom.atom_class;
    atom.input_order = parser->mol->num_atoms;

    int atom_idx = molecule_add_atom_full(parser->mol, &atom);
    if (atom_idx < 0) {
        parser_set_error(parser, "Failed to add atom", token->position);
        return CCHEM_ERROR_MEMORY;
    }

    /* Store neighbor for stereo in original order */
    if (parser->current_atom >= 0) {
        atom_t* prev_atom = molecule_get_atom(parser->mol, parser->current_atom);
        if (prev_atom && prev_atom->num_stereo_neighbors < MAX_NEIGHBORS) {
            prev_atom->stereo_neighbors[prev_atom->num_stereo_neighbors++] = atom_idx;
        }

        atom_t* new_atom = molecule_get_atom(parser->mol, atom_idx);
        if (new_atom) {
            new_atom->stereo_neighbors[new_atom->num_stereo_neighbors++] = parser->current_atom;
        }
    }

    /* Connect to previous atom */
    if (parser->current_atom >= 0) {
        bond_type_t bond_type;

        if (parser->has_pending_bond) {
            bond_type = parser->pending_bond;
            parser->has_pending_bond = false;
        } else {
            /* Default bond type */
            atom_t* prev = molecule_get_atom(parser->mol, parser->current_atom);
            atom_t* curr = molecule_get_atom(parser->mol, atom_idx);

            if (prev && curr && prev->aromatic && curr->aromatic) {
                bond_type = BOND_AROMATIC;
            } else {
                bond_type = BOND_SINGLE;
            }
        }

        int bond_idx = molecule_add_bond(parser->mol, parser->current_atom, atom_idx, bond_type);
        if (bond_idx < 0) {
            parser_set_error(parser, "Failed to add bond", token->position);
            return CCHEM_ERROR_MEMORY;
        }

        /* Track stereo bonds */
        if (bond_is_stereo(bond_type)) {
            if (parser->num_stereo_bonds < MAX_BONDS) {
                parser->stereo_bonds[parser->num_stereo_bonds].atom_idx = atom_idx;
                parser->stereo_bonds[parser->num_stereo_bonds].bond_type = bond_type;
                parser->stereo_bonds[parser->num_stereo_bonds].from_atom = parser->current_atom;
                parser->num_stereo_bonds++;
            }

            bond_t* bond = molecule_get_bond(parser->mol, bond_idx);
            if (bond) {
                bond->stereo_type = bond_type;
                bond->stereo_atom = parser->current_atom;
            }
        }
    }

    parser->current_atom = atom_idx;
    return CCHEM_OK;
}

/* Process ring closure */
static cchem_status_t parser_process_ring_closure(parser_t* parser, const token_t* token) {
    int ring_num = token->data.ring_number;

    if (ring_num < 0 || ring_num >= 100) {
        parser_set_error(parser, "Invalid ring number", token->position);
        return CCHEM_ERROR_RING_CLOSURE;
    }

    if (parser->current_atom < 0) {
        parser_set_error(parser, "Ring closure without atom", token->position);
        return CCHEM_ERROR_RING_CLOSURE;
    }

    /* Track stereo neighbor for ring closure */
    atom_t* curr_atom = molecule_get_atom(parser->mol, parser->current_atom);

    if (parser->ring_open[ring_num]) {
        /* Close the ring */
        int other_atom = parser->ring_closures[ring_num].atom_idx;

        /* Determine bond type */
        bond_type_t bond_type = BOND_SINGLE;

        if (parser->has_pending_bond) {
            bond_type = parser->pending_bond;
            parser->has_pending_bond = false;
        } else if (parser->ring_closures[ring_num].has_bond) {
            bond_type = parser->ring_closures[ring_num].bond_type;
        } else {
            /* Check aromaticity */
            atom_t* other = molecule_get_atom(parser->mol, other_atom);
            if (curr_atom && other && curr_atom->aromatic && other->aromatic) {
                bond_type = BOND_AROMATIC;
            }
        }

        int bond_idx = molecule_add_bond(parser->mol, parser->current_atom, other_atom, bond_type);
        if (bond_idx < 0) {
            parser_set_error(parser, "Failed to add ring bond", token->position);
            return CCHEM_ERROR_MEMORY;
        }

        /* Mark bond as ring bond */
        bond_t* bond = molecule_get_bond(parser->mol, bond_idx);
        if (bond) {
            bond->in_ring = true;
        }

        /* Update ring counts */
        if (curr_atom) curr_atom->ring_count++;
        atom_t* other = molecule_get_atom(parser->mol, other_atom);
        if (other) other->ring_count++;

        /* Handle stereo for ring closure */
        if (parser->ring_closures[ring_num].stereo_type != BOND_NONE) {
            /* Stereo from ring opening */
            if (bond) {
                bond->stereo_type = parser->ring_closures[ring_num].stereo_type;
                bond->stereo_atom = parser->ring_closures[ring_num].stereo_atom;
            }
        } else if (bond_is_stereo(bond_type)) {
            /* Stereo from ring closing (pending bond was stereo) */
            if (bond) {
                bond->stereo_type = bond_type;
                bond->stereo_atom = parser->current_atom;
            }
        }

        /* Add to stereo neighbors */
        if (curr_atom && curr_atom->num_stereo_neighbors < MAX_NEIGHBORS) {
            curr_atom->stereo_neighbors[curr_atom->num_stereo_neighbors++] = other_atom;
        }
        if (other && other->num_stereo_neighbors < MAX_NEIGHBORS) {
            other->stereo_neighbors[other->num_stereo_neighbors++] = parser->current_atom;
        }

        parser->ring_open[ring_num] = false;
        parser->ring_closures[ring_num].atom_idx = -1;
        parser->ring_closures[ring_num].has_bond = false;
        parser->ring_closures[ring_num].stereo_type = BOND_NONE;
    } else {
        /* Open new ring */
        parser->ring_open[ring_num] = true;
        parser->ring_closures[ring_num].atom_idx = parser->current_atom;

        if (parser->has_pending_bond) {
            parser->ring_closures[ring_num].has_bond = true;
            parser->ring_closures[ring_num].bond_type = parser->pending_bond;

            if (bond_is_stereo(parser->pending_bond)) {
                parser->ring_closures[ring_num].stereo_type = parser->pending_bond;
                parser->ring_closures[ring_num].stereo_atom = parser->current_atom;
            }

            parser->has_pending_bond = false;
        } else {
            parser->ring_closures[ring_num].has_bond = false;
            parser->ring_closures[ring_num].stereo_type = BOND_NONE;
        }
    }

    return CCHEM_OK;
}

/* Process bond token */
static cchem_status_t parser_process_bond(parser_t* parser, const token_t* token) {
    bond_type_t bond_type;

    switch (token->type) {
        case TOKEN_BOND_SINGLE:   bond_type = BOND_SINGLE; break;
        case TOKEN_BOND_DOUBLE:   bond_type = BOND_DOUBLE; break;
        case TOKEN_BOND_TRIPLE:   bond_type = BOND_TRIPLE; break;
        case TOKEN_BOND_QUAD:     bond_type = BOND_TRIPLE; break;  /* Treat as triple */
        case TOKEN_BOND_AROMATIC: bond_type = BOND_AROMATIC; break;
        case TOKEN_BOND_UP:       bond_type = BOND_UP; break;
        case TOKEN_BOND_DOWN:     bond_type = BOND_DOWN; break;
        default:
            parser_set_error(parser, "Invalid bond type", token->position);
            return CCHEM_ERROR_PARSE;
    }

    parser->pending_bond = bond_type;
    parser->has_pending_bond = true;

    return CCHEM_OK;
}

/* Process branch open */
static cchem_status_t parser_process_branch_open(parser_t* parser, const token_t* token) {
    if (parser->branch_depth >= MAX_ATOMS) {
        parser_set_error(parser, "Branch depth too deep", token->position);
        return CCHEM_ERROR_PARSE;
    }

    parser->branch_stack[parser->branch_depth++] = parser->current_atom;
    return CCHEM_OK;
}

/* Process branch close */
static cchem_status_t parser_process_branch_close(parser_t* parser, const token_t* token) {
    if (parser->branch_depth <= 0) {
        parser_set_error(parser, "Unmatched branch close", token->position);
        return CCHEM_ERROR_PARSE;
    }

    parser->current_atom = parser->branch_stack[--parser->branch_depth];
    parser->has_pending_bond = false;  /* Clear any pending bond after branch */

    return CCHEM_OK;
}

/* Process dot (disconnection) */
static cchem_status_t parser_process_dot(parser_t* parser, const token_t* token) {
    (void)token;
    parser->current_atom = -1;  /* Next atom starts new fragment */
    parser->has_pending_bond = false;
    return CCHEM_OK;
}

cchem_status_t parser_parse(parser_t* parser, const char* smiles, molecule_t** mol) {
    if (!parser || !smiles || !mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Create new molecule */
    *mol = molecule_create();
    if (!*mol) return CCHEM_ERROR_MEMORY;

    cchem_status_t status = parser_parse_into(parser, smiles, *mol);
    if (status != CCHEM_OK) {
        molecule_free(*mol);
        *mol = NULL;
    }

    return status;
}

cchem_status_t parser_parse_into(parser_t* parser, const char* smiles, molecule_t* mol) {
    if (!parser || !smiles || !mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Reset parser state */
    parser_init(parser);
    parser->mol = mol;

    /* Store original SMILES */
    mol->original_smiles = strdup(smiles);

    /* Initialize lexer */
    cchem_status_t status = lexer_init(&parser->lexer, smiles);
    if (status != CCHEM_OK) {
        parser_set_error(parser, "Failed to initialize lexer", 0);
        return status;
    }

    /* Parse tokens */
    token_t token;
    while (true) {
        status = lexer_next_token(&parser->lexer, &token);
        if (status != CCHEM_OK) {
            parser_set_error(parser, lexer_get_error(&parser->lexer), (int)lexer_get_position(&parser->lexer));
            return status;
        }

        if (token.type == TOKEN_EOF) {
            break;
        }

        switch (token.type) {
            case TOKEN_ATOM_ORGANIC:
            case TOKEN_ATOM_BRACKET:
            case TOKEN_WILDCARD:
                status = parser_process_atom(parser, &token);
                break;

            case TOKEN_RING_CLOSURE:
                status = parser_process_ring_closure(parser, &token);
                break;

            case TOKEN_BOND_SINGLE:
            case TOKEN_BOND_DOUBLE:
            case TOKEN_BOND_TRIPLE:
            case TOKEN_BOND_QUAD:
            case TOKEN_BOND_AROMATIC:
            case TOKEN_BOND_UP:
            case TOKEN_BOND_DOWN:
                status = parser_process_bond(parser, &token);
                break;

            case TOKEN_BRANCH_OPEN:
                status = parser_process_branch_open(parser, &token);
                break;

            case TOKEN_BRANCH_CLOSE:
                status = parser_process_branch_close(parser, &token);
                break;

            case TOKEN_DOT:
                status = parser_process_dot(parser, &token);
                break;

            case TOKEN_ERROR:
                parser_set_error(parser, lexer_get_error(&parser->lexer), token.position);
                return CCHEM_ERROR_PARSE;

            default:
                parser_set_error(parser, "Unexpected token", token.position);
                return CCHEM_ERROR_PARSE;
        }

        if (status != CCHEM_OK) {
            return status;
        }
    }

    /* Check for unclosed rings */
    for (int i = 0; i < 100; i++) {
        if (parser->ring_open[i]) {
            char msg[64];
            snprintf(msg, sizeof(msg), "Unclosed ring %d", i);
            parser_set_error(parser, msg, -1);
            return CCHEM_ERROR_RING_CLOSURE;
        }
    }

    /* Check for unclosed branches */
    if (parser->branch_depth > 0) {
        parser_set_error(parser, "Unclosed branch", -1);
        return CCHEM_ERROR_PARSE;
    }

    /* Calculate implicit hydrogens */
    molecule_calc_implicit_h(mol);

    /* Find fragments */
    molecule_find_fragments(mol);

    return CCHEM_OK;
}

const char* parser_get_error(const parser_t* parser) {
    if (!parser) return "NULL parser";
    return parser->error_msg;
}

int parser_get_error_position(const parser_t* parser) {
    if (!parser) return -1;
    return parser->error_position;
}

/* Convenience function */
molecule_t* smiles_to_molecule(const char* smiles, char* error_buf, size_t error_buf_size) {
    if (!smiles) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "NULL SMILES input");
        }
        return NULL;
    }

    parser_t parser;
    molecule_t* mol = NULL;

    cchem_status_t status = parser_init(&parser);
    if (status != CCHEM_OK) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "Failed to initialize parser");
        }
        return NULL;
    }

    status = parser_parse(&parser, smiles, &mol);
    if (status != CCHEM_OK) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "%s", parser_get_error(&parser));
        }
        return NULL;
    }

    /* Find rings and set ring_count for atoms */
    if (mol && mol->num_atoms > 0) {
        molecule_find_rings(mol);
        /* Perceive aromaticity (converts Kekule to aromatic) */
        molecule_perceive_aromaticity(mol);
    }

    return mol;
}

cchem_status_t smiles_validate(const char* smiles, char* error_buf, size_t error_buf_size) {
    molecule_t* mol = smiles_to_molecule(smiles, error_buf, error_buf_size);
    if (!mol) {
        return CCHEM_ERROR_INVALID_SMILES;
    }

    molecule_free(mol);
    return CCHEM_OK;
}

/* Reuse an existing molecule - avoids malloc/free overhead in batch processing */
cchem_status_t smiles_to_molecule_reuse(molecule_t* mol, const char* smiles,
                                        char* error_buf, size_t error_buf_size) {
    if (!mol || !smiles) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "NULL %s", mol ? "SMILES input" : "molecule");
        }
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Reset molecule state (keeps allocated arrays) */
    molecule_reset(mol);

    parser_t parser;
    cchem_status_t status = parser_init(&parser);
    if (status != CCHEM_OK) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "Failed to initialize parser");
        }
        return status;
    }

    status = parser_parse_into(&parser, smiles, mol);
    if (status != CCHEM_OK) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "%s", parser_get_error(&parser));
        }
        return status;
    }

    /* Find rings and perceive aromaticity */
    if (mol->num_atoms > 0) {
        molecule_find_rings(mol);
        molecule_perceive_aromaticity(mol);
    }

    return CCHEM_OK;
}
