/**
 * @file smiles_writer.c
 * @brief Generate SMILES string from molecular graph
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "cchem/canonicalizer/smiles_writer.h"
#include "cchem/canonicalizer/stereo.h"

/* Forward declarations for ring_info access */
static bool is_ring_opening_target(const smiles_writer_t* writer, int atom_idx, int neighbor);

/*
 * Compute canonical chirality using ring_info for correct neighbor ordering.
 *
 * The SMILES output order for neighbors of a stereocenter is:
 * 1. from_atom (implicit, the atom we came from in DFS)
 * 2. H (if implicit, written in bracket)
 * 3. Ring openings (to unvisited atoms that are ring_info.atom_to where we're atom_from) - sorted by rank
 * 4. Ring closings (to visited atoms that are ring_info.atom_from where we're atom_to) - sorted by rank
 * 5. Branches/chain (to unvisited atoms that are NOT ring opening targets) - sorted by rank
 */
static chirality_t compute_canonical_chirality(const smiles_writer_t* writer, int atom_idx, int from_atom) {
    const molecule_t* mol = writer->mol;
    const atom_t* atom = molecule_get_atom_const(mol, atom_idx);
    if (!atom || atom->chirality == CHIRALITY_NONE) return CHIRALITY_NONE;

    /* Build original neighbor order from parsing */
    int orig_neighbors[4];
    int num_orig = 0;

    /* Determine if original SMILES had a from-atom */
    bool had_from_atom = false;
    if (atom->num_stereo_neighbors > 0) {
        had_from_atom = (from_atom >= 0) || (atom->stereo_neighbors[0] < atom->index);
    }

    if (atom->implicit_h_count > 0 && atom->num_stereo_neighbors > 0) {
        if (had_from_atom) {
            orig_neighbors[num_orig++] = atom->stereo_neighbors[0];
            orig_neighbors[num_orig++] = -1;  /* H at position 1 */
            for (int i = 1; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
                orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
            }
        } else {
            orig_neighbors[num_orig++] = -1;  /* H first */
            for (int i = 0; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
                orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
            }
        }
    } else {
        for (int i = 0; i < atom->num_stereo_neighbors && num_orig < 4; i++) {
            orig_neighbors[num_orig++] = atom->stereo_neighbors[i];
        }
    }

    /* Build canonical neighbor order using ring_info */
    int canon_neighbors[4];
    int num_canon = 0;

    /* 1. From-atom first */
    if (from_atom >= 0) {
        canon_neighbors[num_canon++] = from_atom;
    }

    /* 2. H second (if implicit) */
    if (atom->implicit_h_count > 0) {
        canon_neighbors[num_canon++] = -1;
    }

    /* Collect remaining neighbors into categories */
    int ring_openings[4];   /* Unvisited atoms where we're atom_from in ring_info */
    int num_ring_open = 0;
    int ring_closings[4];   /* Visited atoms */
    int num_ring_close = 0;
    int branches[4];        /* Unvisited atoms that are NOT ring opening targets */
    int num_branches = 0;

    /*
     * IMPORTANT: Ring openings must be collected in ring_info order, NOT atom->neighbors order.
     * This is because write_ring_openings writes them in ring_info order, and we need
     * to match that exact order for stereo calculation.
     */
    for (int r = 0; r < writer->num_ring_info; r++) {
        if (writer->ring_info[r].atom_from == atom_idx) {
            int neighbor = writer->ring_info[r].atom_to;
            if (neighbor != from_atom && !writer->visited[neighbor]) {
                ring_openings[num_ring_open++] = neighbor;
            }
        }
    }

    /* Collect ring closings (visited neighbors) and branches (unvisited non-ring-opening) */
    for (int i = 0; i < atom->num_neighbors; i++) {
        int neighbor = atom->neighbors[i];
        if (neighbor == from_atom) continue;

        if (writer->visited[neighbor]) {
            /* Visited = ring closing */
            ring_closings[num_ring_close++] = neighbor;
        } else {
            /* Unvisited - check if it's a ring opening target (already collected above) */
            if (!is_ring_opening_target(writer, atom_idx, neighbor)) {
                branches[num_branches++] = neighbor;
            }
        }
    }

    /* Sort ring closings and branches by canonical rank.
     * Ring openings are NOT sorted - they stay in ring_info order to match SMILES output. */
    for (int i = 0; i < num_ring_close - 1; i++) {
        for (int j = i + 1; j < num_ring_close; j++) {
            if (mol->atoms[ring_closings[i]].canon_rank > mol->atoms[ring_closings[j]].canon_rank) {
                int tmp = ring_closings[i];
                ring_closings[i] = ring_closings[j];
                ring_closings[j] = tmp;
            }
        }
    }

    for (int i = 0; i < num_branches - 1; i++) {
        for (int j = i + 1; j < num_branches; j++) {
            if (mol->atoms[branches[i]].canon_rank > mol->atoms[branches[j]].canon_rank) {
                int tmp = branches[i];
                branches[i] = branches[j];
                branches[j] = tmp;
            }
        }
    }

    /* 3. Ring openings (to unvisited) come first after H */
    for (int i = 0; i < num_ring_open && num_canon < 4; i++) {
        canon_neighbors[num_canon++] = ring_openings[i];
    }

    /* 4. Ring closings (to visited) come second */
    for (int i = 0; i < num_ring_close && num_canon < 4; i++) {
        canon_neighbors[num_canon++] = ring_closings[i];
    }

    /* 5. Branches/chain (other unvisited) come last */
    for (int i = 0; i < num_branches && num_canon < 4; i++) {
        canon_neighbors[num_canon++] = branches[i];
    }

#if 0  /* Debug output */
    fprintf(stderr, "compute_canonical_chirality: atom %d, from_atom=%d, orig=%d, implicit_h=%d\n",
            atom_idx, from_atom, atom->chirality, atom->implicit_h_count);
    /* Show ring_info order for this atom's ring openings */
    fprintf(stderr, "  ring_info order for openings: ");
    for (int r = 0; r < writer->num_ring_info; r++) {
        if (writer->ring_info[r].atom_from == atom_idx) {
            fprintf(stderr, "%d(ring%d,rank%d) ",
                    writer->ring_info[r].atom_to,
                    writer->ring_info[r].ring_num,
                    mol->atoms[writer->ring_info[r].atom_to].canon_rank);
        }
    }
    fprintf(stderr, "\n  ring_open_targets (ring_info order): ");
    for (int i = 0; i < num_ring_open; i++) fprintf(stderr, "%d(r%d) ", ring_openings[i], mol->atoms[ring_openings[i]].canon_rank);
    fprintf(stderr, "\n  ring_close_targets: ");
    for (int i = 0; i < num_ring_close; i++) fprintf(stderr, "%d(r%d) ", ring_closings[i], mol->atoms[ring_closings[i]].canon_rank);
    fprintf(stderr, "\n  branches: ");
    for (int i = 0; i < num_branches; i++) fprintf(stderr, "%d(r%d) ", branches[i], mol->atoms[branches[i]].canon_rank);
    fprintf(stderr, "\n  orig[%d]: ", num_orig);
    for (int i = 0; i < num_orig; i++) fprintf(stderr, "%d ", orig_neighbors[i]);
    fprintf(stderr, "\n  canon[%d]: ", num_canon);
    for (int i = 0; i < num_canon; i++) fprintf(stderr, "%d ", canon_neighbors[i]);
    fprintf(stderr, "\n");
#endif

    return stereo_permute_chirality(atom->chirality, orig_neighbors, canon_neighbors, num_canon);
}

/* Check if neighbor is a ring opening target from atom_idx */
static bool is_ring_opening_target(const smiles_writer_t* writer, int atom_idx, int neighbor) {
    for (int r = 0; r < writer->num_ring_info; r++) {
        if (writer->ring_info[r].atom_from == atom_idx &&
            writer->ring_info[r].atom_to == neighbor) {
            return true;
        }
    }
    return false;
}

/* Default output options */
const smiles_output_options_t SMILES_OUTPUT_DEFAULT = {
    .canonical = false,
    .kekulize = false,
    .show_explicit_h = false,
    .show_stereo = true,
    .show_isotopes = true,
    .show_charges = true,
    .show_atom_class = false,
    .aromatic_symbols = true,
    .use_brackets = false,
    .max_ring_num = 99
};

const smiles_output_options_t SMILES_OUTPUT_CANONICAL = {
    .canonical = true,
    .kekulize = false,
    .show_explicit_h = false,
    .show_stereo = true,
    .show_isotopes = true,
    .show_charges = true,
    .show_atom_class = false,
    .aromatic_symbols = true,
    .use_brackets = false,
    .max_ring_num = 99
};

#define INITIAL_BUFFER_SIZE 1024

static cchem_status_t buffer_ensure(smiles_writer_t* writer, size_t needed) {
    if (writer->buffer_pos + needed < writer->buffer_size) {
        return CCHEM_OK;
    }

    size_t new_size = writer->buffer_size * 2;
    while (writer->buffer_pos + needed >= new_size) {
        new_size *= 2;
    }

    char* new_buf = (char*)realloc(writer->buffer, new_size);
    if (!new_buf) return CCHEM_ERROR_MEMORY;

    writer->buffer = new_buf;
    writer->buffer_size = new_size;
    return CCHEM_OK;
}

static void buffer_append_char(smiles_writer_t* writer, char c) {
    if (buffer_ensure(writer, 2) == CCHEM_OK) {
        writer->buffer[writer->buffer_pos++] = c;
        writer->buffer[writer->buffer_pos] = '\0';
    }
}

static void buffer_append_str(smiles_writer_t* writer, const char* str) {
    size_t len = strlen(str);
    if (buffer_ensure(writer, len + 1) == CCHEM_OK) {
        memcpy(writer->buffer + writer->buffer_pos, str, len);
        writer->buffer_pos += len;
        writer->buffer[writer->buffer_pos] = '\0';
    }
}

static void buffer_append_int(smiles_writer_t* writer, int val) {
    char buf[16];
    snprintf(buf, sizeof(buf), "%d", val);
    buffer_append_str(writer, buf);
}

smiles_writer_t* smiles_writer_create(const molecule_t* mol,
                                      const smiles_output_options_t* options) {
    if (!mol) return NULL;

    smiles_writer_t* writer = (smiles_writer_t*)calloc(1, sizeof(smiles_writer_t));
    if (!writer) return NULL;

    writer->mol = mol;
    writer->options = options ? *options : SMILES_OUTPUT_DEFAULT;
    writer->atom_order = NULL;

    writer->buffer_size = INITIAL_BUFFER_SIZE;
    writer->buffer = (char*)calloc(writer->buffer_size, sizeof(char));
    writer->buffer_pos = 0;

    writer->visited = (bool*)calloc(mol->num_atoms, sizeof(bool));
    writer->ring_closures = (int*)calloc(mol->num_atoms, sizeof(int));

    if (!writer->buffer || !writer->visited || !writer->ring_closures) {
        smiles_writer_free(writer);
        return NULL;
    }

    writer->buffer[0] = '\0';

    /* Initialize ring number stack */
    writer->next_ring_num = 1;
    writer->ring_num_stack_top = 0;
    for (int i = 1; i <= 99; i++) {
        writer->ring_num_stack[writer->ring_num_stack_top++] = 100 - i;
    }

    writer->num_pending = 0;
    writer->has_error = false;
    writer->error_msg[0] = '\0';

    return writer;
}

void smiles_writer_free(smiles_writer_t* writer) {
    if (!writer) return;

    if (writer->buffer) free(writer->buffer);
    if (writer->visited) free(writer->visited);
    if (writer->ring_closures) free(writer->ring_closures);

    free(writer);
}

void smiles_writer_set_order(smiles_writer_t* writer, const int* order) {
    if (writer) {
        writer->atom_order = order;
    }
}

/* Get lowest available ring number */
static int get_ring_num(smiles_writer_t* writer) {
    if (writer->ring_num_stack_top > 0) {
        return writer->ring_num_stack[--writer->ring_num_stack_top];
    }
    return writer->next_ring_num++;
}

/* Check if atom needs brackets */
bool smiles_atom_needs_brackets(const molecule_t* mol, int atom_idx,
                                const smiles_output_options_t* options) {
    const atom_t* atom = molecule_get_atom_const(mol, atom_idx);
    if (!atom) return true;

    if (options->use_brackets) return true;

    /* Organic subset atoms may not need brackets */
    bool organic_subset = element_is_organic_subset(atom->element);

    /* Needs brackets if not in organic subset */
    if (!organic_subset) return true;

    /* Needs brackets if has isotope */
    if (atom->isotope > 0 && options->show_isotopes) return true;

    /* Needs brackets if has non-default charge */
    if (atom->charge != 0 && options->show_charges) return true;

    /* Needs brackets if has explicit H count specified */
    if (atom->h_count >= 0) return true;

    /* Aromatic heteroatoms with implicit hydrogens need brackets (e.g., [nH] for pyrrole) */
    if (atom->aromatic && atom->implicit_h_count > 0) {
        /* Aromatic N, O, S with hydrogens need brackets */
        if (atom->element == ELEM_N || atom->element == ELEM_O || atom->element == ELEM_S) {
            return true;
        }
    }

    /* Needs brackets if has chirality */
    if (atom->chirality != CHIRALITY_NONE && options->show_stereo) return true;

    /* Needs brackets if has atom class */
    if (atom->atom_class > 0 && options->show_atom_class) return true;

    return false;
}

/* Get element symbol (lowercase for aromatic) */
const char* smiles_get_element_symbol(element_t elem, bool aromatic) {
    /* Use thread-local storage to avoid race conditions in parallel processing */
    static _Thread_local char aromatic_symbol[4];

    const char* symbol = element_to_symbol(elem);

    if (aromatic && strlen(symbol) > 0) {
        aromatic_symbol[0] = tolower(symbol[0]);
        aromatic_symbol[1] = symbol[1];
        aromatic_symbol[2] = '\0';
        return aromatic_symbol;
    }

    return symbol;
}

void smiles_format_charge(int charge, char* buf, size_t buf_size) {
    if (charge == 0) {
        buf[0] = '\0';
        return;
    }

    if (charge == 1) {
        snprintf(buf, buf_size, "+");
    } else if (charge == -1) {
        snprintf(buf, buf_size, "-");
    } else if (charge > 0) {
        snprintf(buf, buf_size, "+%d", charge);
    } else {
        snprintf(buf, buf_size, "%d", charge);
    }
}

void smiles_format_h_count(int h_count, char* buf, size_t buf_size) {
    if (h_count <= 0) {
        buf[0] = '\0';
    } else if (h_count == 1) {
        snprintf(buf, buf_size, "H");
    } else {
        snprintf(buf, buf_size, "H%d", h_count);
    }
}

/* Write a single atom */
cchem_status_t smiles_write_atom(smiles_writer_t* writer, int atom_idx, int from_atom) {
    const atom_t* atom = molecule_get_atom_const(writer->mol, atom_idx);
    if (!atom) return CCHEM_ERROR_INVALID_INPUT;

    bool need_brackets = smiles_atom_needs_brackets(writer->mol, atom_idx, &writer->options);

    if (need_brackets) {
        buffer_append_char(writer, '[');

        /* Isotope */
        if (atom->isotope > 0 && writer->options.show_isotopes) {
            buffer_append_int(writer, atom->isotope);
        }

        /* Element symbol */
        const char* symbol = smiles_get_element_symbol(atom->element,
                              atom->aromatic && writer->options.aromatic_symbols);
        buffer_append_str(writer, symbol);

        /* Chirality */
        if (atom->chirality != CHIRALITY_NONE && writer->options.show_stereo) {
            chirality_t chiral = atom->chirality;
            if (writer->mol->is_canonical) {
                /* Use our function that correctly handles ring opening order */
                chiral = compute_canonical_chirality(writer, atom_idx, from_atom);
            }

            switch (chiral) {
                case CHIRALITY_CW:  buffer_append_char(writer, '@'); break;
                case CHIRALITY_CCW: buffer_append_str(writer, "@@"); break;
                case CHIRALITY_TH1: buffer_append_str(writer, "@TH1"); break;
                case CHIRALITY_TH2: buffer_append_str(writer, "@TH2"); break;
                default: break;
            }
        }

        /* Hydrogen count */
        if (atom->h_count >= 0 || atom->implicit_h_count > 0) {
            int h = (atom->h_count >= 0) ? atom->h_count : atom->implicit_h_count;
            if (h > 0) {
                char h_buf[8];
                smiles_format_h_count(h, h_buf, sizeof(h_buf));
                buffer_append_str(writer, h_buf);
            }
        }

        /* Charge */
        if (atom->charge != 0 && writer->options.show_charges) {
            char charge_buf[8];
            smiles_format_charge(atom->charge, charge_buf, sizeof(charge_buf));
            buffer_append_str(writer, charge_buf);
        }

        /* Atom class */
        if (atom->atom_class > 0 && writer->options.show_atom_class) {
            buffer_append_char(writer, ':');
            buffer_append_int(writer, atom->atom_class);
        }

        buffer_append_char(writer, ']');
    } else {
        /* Simple organic subset atom */
        const char* symbol = smiles_get_element_symbol(atom->element,
                              atom->aromatic && writer->options.aromatic_symbols);
        buffer_append_str(writer, symbol);
    }

    return CCHEM_OK;
}

/* Write bond symbol */
cchem_status_t smiles_write_bond(smiles_writer_t* writer, bond_type_t type,
                                 bool aromatic_context) {
    /* Don't write single bonds in normal context */
    if (type == BOND_SINGLE && !aromatic_context) return CCHEM_OK;

    /* Don't write aromatic bonds between aromatic atoms */
    if (type == BOND_AROMATIC && aromatic_context) return CCHEM_OK;

    char c = bond_to_smiles_char(type);
    if (c != '\0') {
        buffer_append_char(writer, c);
    }

    return CCHEM_OK;
}

/* Write ring closure */
cchem_status_t smiles_write_ring_closure(smiles_writer_t* writer, int ring_num,
                                         bond_type_t bond_type) {
    /* Write bond type if not single/aromatic */
    if (bond_type == BOND_DOUBLE) {
        buffer_append_char(writer, '=');
    } else if (bond_type == BOND_TRIPLE) {
        buffer_append_char(writer, '#');
    }

    if (ring_num < 10) {
        buffer_append_char(writer, '0' + ring_num);
    } else {
        buffer_append_char(writer, '%');
        buffer_append_char(writer, '0' + (ring_num / 10));
        buffer_append_char(writer, '0' + (ring_num % 10));
    }

    return CCHEM_OK;
}

/* Check if we need to write ring opening at this atom */
static void write_ring_openings(smiles_writer_t* writer, int atom_idx) {
    for (int i = 0; i < writer->num_ring_info; i++) {
        if (writer->ring_info[i].atom_from == atom_idx && !writer->ring_info[i].written_from) {
            smiles_write_ring_closure(writer, writer->ring_info[i].ring_num, BOND_SINGLE);
            writer->ring_info[i].written_from = true;
        }
    }
}

/* Check if we need to write ring closing at this atom */
static void write_ring_closings(smiles_writer_t* writer, int atom_idx) {
    for (int i = 0; i < writer->num_ring_info; i++) {
        if (writer->ring_info[i].atom_to == atom_idx && !writer->ring_info[i].written_to) {
            /* Write bond type only if not single/aromatic */
            bond_type_t bt = writer->ring_info[i].bond_type;
            bond_type_t st = writer->ring_info[i].stereo_type;

            if (bt == BOND_DOUBLE) {
                buffer_append_char(writer, '=');
            } else if (bt == BOND_TRIPLE) {
                buffer_append_char(writer, '#');
            }

            /* Write stereo marker for E/Z double bonds */
            if (st != BOND_NONE && writer->options.show_stereo) {
                /* Determine correct direction based on which atom we're at */
                if (writer->ring_info[i].stereo_atom == atom_idx) {
                    buffer_append_char(writer, bond_to_smiles_char(st));
                } else {
                    /* Flip direction if we're at the other atom */
                    if (st == BOND_UP) {
                        buffer_append_char(writer, bond_to_smiles_char(BOND_DOWN));
                    } else if (st == BOND_DOWN) {
                        buffer_append_char(writer, bond_to_smiles_char(BOND_UP));
                    }
                }
            }

            smiles_write_ring_closure(writer, writer->ring_info[i].ring_num, BOND_SINGLE);
            writer->ring_info[i].written_to = true;
        }
    }
}

/* First pass: identify ring closures */
static void identify_rings_dfs(smiles_writer_t* writer, int atom_idx, int from_atom, int* dfs_order, int* order_idx) {
    const atom_t* atom = molecule_get_atom_const(writer->mol, atom_idx);
    if (!atom) return;

    writer->visited[atom_idx] = true;
    dfs_order[atom_idx] = (*order_idx)++;

    /* Collect and sort neighbors by canonical rank (same as write_dfs) */
    typedef struct {
        int neighbor_idx;
        int bond_idx;
        int canon_rank;
    } nb_info_t;

    nb_info_t neighbors[16];
    int num_neighbors = 0;

    for (int i = 0; i < atom->num_neighbors && num_neighbors < 16; i++) {
        int neighbor = atom->neighbors[i];
        if (neighbor == from_atom) continue;

        neighbors[num_neighbors].neighbor_idx = neighbor;
        neighbors[num_neighbors].bond_idx = atom->neighbor_bonds[i];
        neighbors[num_neighbors].canon_rank = writer->mol->atoms[neighbor].canon_rank;
        num_neighbors++;
    }

    /* Sort by canonical rank */
    for (int i = 0; i < num_neighbors - 1; i++) {
        for (int j = i + 1; j < num_neighbors; j++) {
            if (neighbors[i].canon_rank > neighbors[j].canon_rank) {
                nb_info_t tmp = neighbors[i];
                neighbors[i] = neighbors[j];
                neighbors[j] = tmp;
            }
        }
    }

    /* Check for ring closures */
    for (int i = 0; i < num_neighbors; i++) {
        int neighbor = neighbors[i].neighbor_idx;

        if (writer->visited[neighbor]) {
            /* This is a back edge - ring closure */
            int bond_idx = neighbors[i].bond_idx;
            const bond_t* bond = molecule_get_bond_const(writer->mol, bond_idx);

            if (writer->num_ring_info < 100) {
                int ring_num = get_ring_num(writer);
                writer->ring_info[writer->num_ring_info].atom_from = neighbor;  /* Earlier atom */
                writer->ring_info[writer->num_ring_info].atom_to = atom_idx;    /* Current atom */
                writer->ring_info[writer->num_ring_info].ring_num = ring_num;
                writer->ring_info[writer->num_ring_info].bond_type = bond ? bond->type : BOND_SINGLE;
                writer->ring_info[writer->num_ring_info].stereo_type = bond ? bond->stereo_type : BOND_NONE;
                writer->ring_info[writer->num_ring_info].stereo_atom = bond ? bond->stereo_atom : -1;
                writer->ring_info[writer->num_ring_info].written_from = false;
                writer->ring_info[writer->num_ring_info].written_to = false;
                writer->num_ring_info++;
            }
        }
    }

    /* Recurse to unvisited neighbors in canonical order */
    for (int i = 0; i < num_neighbors; i++) {
        int neighbor = neighbors[i].neighbor_idx;
        if (writer->visited[neighbor]) continue;

        identify_rings_dfs(writer, neighbor, atom_idx, dfs_order, order_idx);
    }
}

/* Structure for sorting neighbors by canonical rank */
typedef struct {
    int neighbor_idx;
    int bond_idx;
    int canon_rank;
} neighbor_info_t;

static int neighbor_compare(const void* a, const void* b) {
    const neighbor_info_t* na = (const neighbor_info_t*)a;
    const neighbor_info_t* nb = (const neighbor_info_t*)b;
    return na->canon_rank - nb->canon_rank;
}

/* DFS traversal to generate SMILES */
static void write_dfs(smiles_writer_t* writer, int atom_idx, int from_atom, bond_type_t from_bond) {
    const atom_t* atom = molecule_get_atom_const(writer->mol, atom_idx);
    if (!atom) return;

    writer->visited[atom_idx] = true;

    /* Write bond from parent (if not first atom) */
    if (from_atom >= 0) {
        const atom_t* from = molecule_get_atom_const(writer->mol, from_atom);
        bool aromatic_context = (atom->aromatic && from->aromatic);

        /* Handle stereo bonds */
        if (bond_is_stereo(from_bond) && writer->options.show_stereo) {
            buffer_append_char(writer, bond_to_smiles_char(from_bond));
        } else {
            smiles_write_bond(writer, from_bond, aromatic_context);
        }
    }

    /* Write atom (pass from_atom for chirality calculation) */
    smiles_write_atom(writer, atom_idx, from_atom);

    /* Write any ring openings at this atom */
    write_ring_openings(writer, atom_idx);

    /* Write any ring closings at this atom */
    write_ring_closings(writer, atom_idx);

    /* Collect unvisited neighbors and sort by canonical rank */
    neighbor_info_t neighbors[16];  /* Max 16 neighbors */
    int num_unvisited = 0;

    for (int i = 0; i < atom->num_neighbors && num_unvisited < 16; i++) {
        int neighbor = atom->neighbors[i];

        if (neighbor == from_atom) continue;
        if (writer->visited[neighbor]) continue;

        /* Skip neighbors that are ring closure targets from this atom */
        bool is_ring_target = false;
        for (int r = 0; r < writer->num_ring_info; r++) {
            if (writer->ring_info[r].atom_from == atom_idx && writer->ring_info[r].atom_to == neighbor) {
                is_ring_target = true;
                break;
            }
        }
        if (is_ring_target) continue;

        neighbors[num_unvisited].neighbor_idx = neighbor;
        neighbors[num_unvisited].bond_idx = atom->neighbor_bonds[i];
        neighbors[num_unvisited].canon_rank = writer->mol->atoms[neighbor].canon_rank;
        num_unvisited++;
    }

    /* Sort neighbors by canonical rank */
    if (num_unvisited > 1) {
        qsort(neighbors, num_unvisited, sizeof(neighbor_info_t), neighbor_compare);
    }

    /*
     * Visit neighbors: write branches first (in parentheses), then main chain last.
     * In SMILES, branches must appear BEFORE the main chain continuation because
     * the parser attaches branches to the current atom on the stack.
     *
     * If we have neighbors sorted by rank [a, b, c] (lowest rank first),
     * the last one (c) becomes the main chain continuation, and [a, b] are branches.
     * We write: (a)(b)c so that branches attach to the current atom before we move to c.
     *
     * Exception: if there's only one neighbor, it's the main chain (no parentheses).
     */
    if (num_unvisited == 1) {
        /* Single neighbor - just continue main chain */
        int neighbor = neighbors[0].neighbor_idx;
        int bond_idx = neighbors[0].bond_idx;

        const bond_t* bond = molecule_get_bond_const(writer->mol, bond_idx);
        bond_type_t bt = bond ? bond->type : BOND_SINGLE;

        if (bond && bond->stereo_type != BOND_NONE && writer->options.show_stereo) {
            if (bond->stereo_atom == atom_idx) {
                bt = bond->stereo_type;
            } else {
                if (bond->stereo_type == BOND_UP) bt = BOND_DOWN;
                else if (bond->stereo_type == BOND_DOWN) bt = BOND_UP;
            }
        }

        write_dfs(writer, neighbor, atom_idx, bt);
    } else if (num_unvisited > 1) {
        /* Multiple neighbors: write branches first, main chain last */
        /* Main chain is the LAST neighbor (highest rank among unvisited) */
        int main_chain_idx = num_unvisited - 1;

        /* Write all branches first (all except the last one) */
        for (int i = 0; i < main_chain_idx; i++) {
            int neighbor = neighbors[i].neighbor_idx;
            int bond_idx = neighbors[i].bond_idx;

            const bond_t* bond = molecule_get_bond_const(writer->mol, bond_idx);
            bond_type_t bt = bond ? bond->type : BOND_SINGLE;

            if (bond && bond->stereo_type != BOND_NONE && writer->options.show_stereo) {
                if (bond->stereo_atom == atom_idx) {
                    bt = bond->stereo_type;
                } else {
                    if (bond->stereo_type == BOND_UP) bt = BOND_DOWN;
                    else if (bond->stereo_type == BOND_DOWN) bt = BOND_UP;
                }
            }

            buffer_append_char(writer, '(');
            write_dfs(writer, neighbor, atom_idx, bt);
            buffer_append_char(writer, ')');
        }

        /* Write main chain last (no parentheses) */
        {
            int neighbor = neighbors[main_chain_idx].neighbor_idx;
            int bond_idx = neighbors[main_chain_idx].bond_idx;

            const bond_t* bond = molecule_get_bond_const(writer->mol, bond_idx);
            bond_type_t bt = bond ? bond->type : BOND_SINGLE;

            if (bond && bond->stereo_type != BOND_NONE && writer->options.show_stereo) {
                if (bond->stereo_atom == atom_idx) {
                    bt = bond->stereo_type;
                } else {
                    if (bond->stereo_type == BOND_UP) bt = BOND_DOWN;
                    else if (bond->stereo_type == BOND_DOWN) bt = BOND_UP;
                }
            }

            write_dfs(writer, neighbor, atom_idx, bt);
        }
    }
}

cchem_status_t smiles_writer_write(smiles_writer_t* writer) {
    if (!writer || !writer->mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Reset state */
    writer->buffer_pos = 0;
    writer->buffer[0] = '\0';
    writer->num_pending = 0;

    /* Reset ring closure tracking */
    writer->num_ring_info = 0;

    /* Reset ring number stack */
    writer->next_ring_num = 1;
    writer->ring_num_stack_top = 0;
    for (int i = 1; i <= 99; i++) {
        writer->ring_num_stack[writer->ring_num_stack_top++] = 100 - i;
    }

    int num_atoms = writer->mol->num_atoms;
    if (num_atoms <= 0) {
        return CCHEM_OK;
    }

    memset(writer->visited, 0, (size_t)num_atoms * sizeof(bool));

    /* Determine traversal order */
    const int* order = writer->atom_order;
    int* default_order = NULL;

    if (!order) {
        default_order = (int*)calloc(num_atoms, sizeof(int));
        if (!default_order) return CCHEM_ERROR_MEMORY;
        for (int i = 0; i < num_atoms; i++) {
            default_order[i] = i;
        }
        order = default_order;
    }

    /* Allocate DFS order tracking array */
    int* dfs_order = (int*)calloc(num_atoms, sizeof(int));
    if (!dfs_order) {
        if (default_order) free(default_order);
        return CCHEM_ERROR_MEMORY;
    }
    memset(dfs_order, -1, (size_t)num_atoms * sizeof(int));

    /* FIRST PASS: Identify all ring closures */
    int order_idx = 0;
    for (int i = 0; i < num_atoms; i++) {
        int start_atom = order[i];
        if (writer->visited[start_atom]) continue;
        identify_rings_dfs(writer, start_atom, -1, dfs_order, &order_idx);
    }

    free(dfs_order);

    /* Reset visited for second pass */
    memset(writer->visited, 0, (size_t)num_atoms * sizeof(bool));

    /* SECOND PASS: Write SMILES with ring closures */
    bool first_fragment = true;

    for (int i = 0; i < num_atoms; i++) {
        int start_atom = order[i];

        if (writer->visited[start_atom]) continue;

        if (!first_fragment) {
            buffer_append_char(writer, '.');
        }
        first_fragment = false;

        write_dfs(writer, start_atom, -1, BOND_NONE);
    }

    if (default_order) free(default_order);

    return CCHEM_OK;
}

const char* smiles_writer_get_smiles(const smiles_writer_t* writer) {
    if (!writer) return NULL;
    return writer->buffer;
}

char* smiles_writer_get_smiles_copy(const smiles_writer_t* writer) {
    if (!writer || !writer->buffer) return NULL;
    return strdup(writer->buffer);
}

const char* smiles_writer_get_error(const smiles_writer_t* writer) {
    if (!writer) return "NULL writer";
    return writer->error_msg;
}

/* High-level function */
char* molecule_to_smiles(const molecule_t* mol, const smiles_output_options_t* options) {
    if (!mol) return NULL;

    smiles_writer_t* writer = smiles_writer_create(mol, options);
    if (!writer) return NULL;

    /* Use canonical order if available */
    if (mol->canon_order && options && options->canonical) {
        smiles_writer_set_order(writer, mol->canon_order);
    }

    cchem_status_t status = smiles_writer_write(writer);
    if (status != CCHEM_OK) {
        smiles_writer_free(writer);
        return NULL;
    }

    char* result = smiles_writer_get_smiles_copy(writer);
    smiles_writer_free(writer);

    return result;
}

char* molecule_to_canonical_smiles_str(const molecule_t* mol) {
    return molecule_to_smiles(mol, &SMILES_OUTPUT_CANONICAL);
}
