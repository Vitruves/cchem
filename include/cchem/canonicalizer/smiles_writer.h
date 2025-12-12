/**
 * @file smiles_writer.h
 * @brief Generate SMILES string from molecular graph
 *
 * Generates canonical SMILES using DFS traversal with:
 * - Optimal ring closure numbering
 * - Minimal bracket atom usage
 * - Proper stereochemistry encoding
 */

#ifndef CCHEM_CANONICALIZER_SMILES_WRITER_H
#define CCHEM_CANONICALIZER_SMILES_WRITER_H

#include "types.h"
#include "molecule.h"
#include "canon.h"

/* SMILES output options */
typedef struct {
    bool canonical;              /* Use canonical ordering */
    bool kekulize;               /* Output Kekule form */
    bool show_explicit_h;        /* Always show explicit H */
    bool show_stereo;            /* Include stereochemistry */
    bool show_isotopes;          /* Include isotopes */
    bool show_charges;           /* Include charges */
    bool show_atom_class;        /* Include atom class [C:1] */
    bool aromatic_symbols;       /* Use lowercase for aromatic atoms */
    bool use_brackets;           /* Force brackets on all atoms */
    int max_ring_num;            /* Maximum ring closure number to use */
} smiles_output_options_t;

/* Default SMILES output options */
extern const smiles_output_options_t SMILES_OUTPUT_DEFAULT;

/* Canonical output options */
extern const smiles_output_options_t SMILES_OUTPUT_CANONICAL;

/* Ring closure info for DFS tracking */
typedef struct {
    int atom_from;     /* Atom where ring opens */
    int atom_to;       /* Atom where ring closes */
    int ring_num;      /* Assigned ring number */
    bond_type_t bond_type;
    bond_type_t stereo_type;  /* Stereo direction (UP/DOWN) for E/Z stereo */
    int stereo_atom;          /* Atom that defines the stereo direction */
    bool written_from; /* Written at opening atom */
    bool written_to;   /* Written at closing atom */
} ring_closure_info_t;

/* SMILES writer context */
typedef struct {
    const molecule_t* mol;       /* Source molecule */
    const int* atom_order;       /* Atom traversal order */
    smiles_output_options_t options;

    /* Output buffer */
    char* buffer;
    size_t buffer_size;
    size_t buffer_pos;

    /* Traversal state */
    bool* visited;               /* Visited atoms */
    int* ring_closures;          /* Ring closure assignments */
    int next_ring_num;           /* Next available ring number */
    int ring_num_stack[100];     /* Available ring numbers */
    int ring_num_stack_top;

    /* Ring closure tracking - moved from static to instance for thread safety */
    ring_closure_info_t ring_info[100];
    int num_ring_info;

    /* Legacy pending closures (unused, kept for compatibility) */
    struct {
        int atom_idx;
        int ring_num;
        bond_type_t bond_type;
    } pending_closures[100];
    int num_pending;

    /* Error handling */
    char error_msg[256];
    bool has_error;
} smiles_writer_t;

/* Create SMILES writer */
smiles_writer_t* smiles_writer_create(const molecule_t* mol,
                                      const smiles_output_options_t* options);

/* Free SMILES writer */
void smiles_writer_free(smiles_writer_t* writer);

/* Set atom order for canonical output */
void smiles_writer_set_order(smiles_writer_t* writer, const int* order);

/* Generate SMILES string */
cchem_status_t smiles_writer_write(smiles_writer_t* writer);

/* Get generated SMILES (owned by writer) */
const char* smiles_writer_get_smiles(const smiles_writer_t* writer);

/* Get copy of generated SMILES (caller must free) */
char* smiles_writer_get_smiles_copy(const smiles_writer_t* writer);

/* Get error message */
const char* smiles_writer_get_error(const smiles_writer_t* writer);

/* High-level: molecule to SMILES string */
char* molecule_to_smiles(const molecule_t* mol,
                         const smiles_output_options_t* options);

/* High-level: molecule to canonical SMILES */
char* molecule_to_canonical_smiles_str(const molecule_t* mol);

/* Internal: write single atom (from_atom is the atom we came from in DFS, -1 if starting) */
cchem_status_t smiles_write_atom(smiles_writer_t* writer, int atom_idx, int from_atom);

/* Internal: write bond symbol */
cchem_status_t smiles_write_bond(smiles_writer_t* writer, bond_type_t type,
                                 bool aromatic_context);

/* Internal: write ring closure */
cchem_status_t smiles_write_ring_closure(smiles_writer_t* writer, int ring_num,
                                         bond_type_t bond_type);

/* Check if atom needs brackets */
bool smiles_atom_needs_brackets(const molecule_t* mol, int atom_idx,
                                const smiles_output_options_t* options);

/* Get canonical element symbol (lowercase if aromatic) */
const char* smiles_get_element_symbol(element_t elem, bool aromatic);

/* Format charge string for bracket atom */
void smiles_format_charge(int charge, char* buf, size_t buf_size);

/* Format hydrogen count for bracket atom */
void smiles_format_h_count(int h_count, char* buf, size_t buf_size);

#endif /* CCHEM_CANONICALIZER_SMILES_WRITER_H */
