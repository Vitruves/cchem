/**
 * @file lexer.h
 * @brief SMILES lexer - tokenizes SMILES strings
 *
 * Full SMILES grammar support including:
 * - Organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I)
 * - Bracket atoms with isotope, element, chirality, hcount, charge, class
 * - Bond symbols (-, =, #, $, :, /, \)
 * - Ring closures (digits 0-9, %nn)
 * - Branches ( )
 * - Dot disconnection .
 * - Aromatic atoms (b, c, n, o, p, s)
 */

#ifndef CCHEM_CANONICALIZER_LEXER_H
#define CCHEM_CANONICALIZER_LEXER_H

#include "types.h"
#include "element.h"
#include "bond.h"

/* Token types */
typedef enum {
    TOKEN_EOF = 0,
    TOKEN_ERROR,

    /* Atoms */
    TOKEN_ATOM_ORGANIC,     /* B, C, N, O, P, S, F, Cl, Br, I (and aromatic) */
    TOKEN_ATOM_BRACKET,     /* [xxx] complete bracket atom */

    /* Bonds */
    TOKEN_BOND_SINGLE,      /* - */
    TOKEN_BOND_DOUBLE,      /* = */
    TOKEN_BOND_TRIPLE,      /* # */
    TOKEN_BOND_QUAD,        /* $ */
    TOKEN_BOND_AROMATIC,    /* : */
    TOKEN_BOND_UP,          /* / */
    TOKEN_BOND_DOWN,        /* \ */

    /* Structure */
    TOKEN_BRANCH_OPEN,      /* ( */
    TOKEN_BRANCH_CLOSE,     /* ) */
    TOKEN_DOT,              /* . */
    TOKEN_RING_CLOSURE,     /* digit or %nn */

    /* Wildcards */
    TOKEN_WILDCARD          /* * */
} token_type_t;

/* Token data for atoms */
typedef struct {
    element_t element;
    int isotope;
    int charge;
    int h_count;           /* -1 means not specified */
    chirality_t chirality;
    int chirality_class;   /* For OH, TB extended chirality */
    bool aromatic;
    int atom_class;        /* [C:1] -> class = 1 */
    bool in_bracket;       /* Was in brackets */
    char symbol[4];        /* Original symbol */
} token_atom_data_t;

/* Token structure */
typedef struct {
    token_type_t type;
    int position;           /* Position in input string */
    int length;             /* Length of token in input */

    union {
        token_atom_data_t atom;  /* Atom data */
        int ring_number;         /* Ring closure number */
        bond_type_t bond_type;   /* Bond type */
    } data;
} token_t;

/* Lexer state */
typedef struct {
    const char* input;      /* Input SMILES string */
    size_t input_len;       /* Length of input */
    size_t pos;             /* Current position */
    size_t line;            /* Line number (always 1 for SMILES) */
    size_t column;          /* Column number */

    /* Current token */
    token_t current;

    /* Error handling */
    char error_msg[256];
    bool has_error;
} lexer_t;

/* Initialize lexer with input string */
cchem_status_t lexer_init(lexer_t* lexer, const char* input);

/* Reset lexer to beginning */
void lexer_reset(lexer_t* lexer);

/* Get next token */
cchem_status_t lexer_next_token(lexer_t* lexer, token_t* token);

/* Check if at end of input */
bool lexer_at_end(const lexer_t* lexer);

/* Get current position */
size_t lexer_get_position(const lexer_t* lexer);

/* Get error message */
const char* lexer_get_error(const lexer_t* lexer);

#endif /* CCHEM_CANONICALIZER_LEXER_H */
