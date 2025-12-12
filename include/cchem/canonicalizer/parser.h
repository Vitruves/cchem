/**
 * @file parser.h
 * @brief SMILES parser - builds molecular graph from tokens
 */

#ifndef CCHEM_CANONICALIZER_PARSER_H
#define CCHEM_CANONICALIZER_PARSER_H

#include "types.h"
#include "lexer.h"
#include "molecule.h"

/* Parser state */
typedef struct {
    lexer_t lexer;              /* Lexer instance */
    molecule_t* mol;            /* Molecule being built */

    /* Ring closure tracking */
    struct {
        int atom_idx;           /* Atom where ring opened */
        bond_type_t bond_type;  /* Bond type if specified */
        bool has_bond;          /* Was bond type specified */
        int stereo_atom;        /* For stereo bonds */
        bond_type_t stereo_type;
    } ring_closures[100];       /* Ring closure stack (0-99) */
    bool ring_open[100];        /* Which ring numbers are open */

    /* Branch stack */
    int branch_stack[MAX_ATOMS]; /* Stack of atom indices */
    int branch_depth;            /* Current branch depth */

    /* Current state */
    int current_atom;           /* Last atom added */
    bond_type_t pending_bond;   /* Bond type for next connection */
    bool has_pending_bond;      /* Is there a pending bond */

    /* Stereochemistry tracking */
    struct {
        int atom_idx;
        bond_type_t bond_type;
        int from_atom;
    } stereo_bonds[MAX_BONDS];
    int num_stereo_bonds;

    /* Error handling */
    char error_msg[512];
    bool has_error;
    int error_position;
} parser_t;

/* Initialize parser */
cchem_status_t parser_init(parser_t* parser);

/* Free parser resources */
void parser_free(parser_t* parser);

/* Parse SMILES string into molecule */
cchem_status_t parser_parse(parser_t* parser, const char* smiles, molecule_t** mol);

/* Parse SMILES string into existing molecule */
cchem_status_t parser_parse_into(parser_t* parser, const char* smiles, molecule_t* mol);

/* Get error message */
const char* parser_get_error(const parser_t* parser);

/* Get error position */
int parser_get_error_position(const parser_t* parser);

/* Convenience function: parse SMILES directly */
molecule_t* smiles_to_molecule(const char* smiles, char* error_buf, size_t error_buf_size);

/* Parse SMILES into existing molecule (reuse allocated memory) */
cchem_status_t smiles_to_molecule_reuse(molecule_t* mol, const char* smiles,
                                        char* error_buf, size_t error_buf_size);

/* Validate SMILES syntax without building molecule */
cchem_status_t smiles_validate(const char* smiles, char* error_buf, size_t error_buf_size);

#endif /* CCHEM_CANONICALIZER_PARSER_H */
