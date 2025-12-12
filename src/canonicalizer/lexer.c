/**
 * @file lexer.c
 * @brief SMILES lexer implementation
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include "cchem/canonicalizer/lexer.h"

/* Lookup tables for fast character classification (avoid branches) */
static const uint8_t organic_subset_table[128] = {
    ['B']=1,['C']=1,['N']=1,['O']=1,['P']=1,['S']=1,['F']=1,['I']=1,
    ['b']=1,['c']=1,['n']=1,['o']=1,['p']=1,['s']=1
};
static const uint8_t bond_char_table[128] = {
    ['-']=1,['=']=1,['#']=1,['$']=1,[':']=1,['/']=1,['\\']=1
};

/* Helper: check if character is start of organic subset atom */
static inline bool is_organic_subset_start(char c) {
    return (unsigned char)c < 128 && organic_subset_table[(unsigned char)c];
}

/* Helper: check if character is bond symbol */
static inline bool is_bond_char(char c) {
    return (unsigned char)c < 128 && bond_char_table[(unsigned char)c];
}

/* Helper: peek character at offset */
static char lexer_peek_at(const lexer_t* lexer, size_t offset) {
    if (lexer->pos + offset >= lexer->input_len) return '\0';
    return lexer->input[lexer->pos + offset];
}

/* Helper: current character */
static char lexer_current(const lexer_t* lexer) {
    return lexer_peek_at(lexer, 0);
}

/* Helper: advance position */
static void lexer_advance(lexer_t* lexer) {
    if (lexer->pos < lexer->input_len) {
        lexer->pos++;
        lexer->column++;
    }
}

/* Helper: advance by n characters */
static inline void lexer_advance_by(lexer_t* lexer, size_t n) {
    size_t remaining = lexer->input_len - lexer->pos;
    size_t advance = (n < remaining) ? n : remaining;
    lexer->pos += advance;
    lexer->column += advance;
}

/* Parse element from symbol */
static element_t parse_element(const char* symbol, bool* aromatic) {
    *aromatic = false;

    /* Check for aromatic lowercase */
    if (islower(symbol[0])) {
        *aromatic = true;
        char upper[4];
        upper[0] = toupper(symbol[0]);
        upper[1] = symbol[1];
        upper[2] = symbol[2];
        upper[3] = '\0';

        /* Special case: 'c' not 'C' */
        if (strcmp(upper, "C") == 0) return ELEM_C;
        if (strcmp(upper, "N") == 0) return ELEM_N;
        if (strcmp(upper, "O") == 0) return ELEM_O;
        if (strcmp(upper, "S") == 0) return ELEM_S;
        if (strcmp(upper, "P") == 0) return ELEM_P;
        if (strcmp(upper, "B") == 0) return ELEM_B;
        if (strcmp(upper, "Se") == 0) return ELEM_Se;
        if (strcmp(upper, "As") == 0) return ELEM_As;
    }

    return element_from_symbol(symbol);
}

/* Parse organic subset atom */
static cchem_status_t lexer_parse_organic_atom(lexer_t* lexer, token_t* token) {
    char c = lexer_current(lexer);
    char next = lexer_peek_at(lexer, 1);
    bool aromatic = islower(c);

    token->type = TOKEN_ATOM_ORGANIC;
    token->position = (int)lexer->pos;
    token->data.atom.in_bracket = false;
    token->data.atom.isotope = 0;
    token->data.atom.charge = 0;
    token->data.atom.h_count = -1;
    token->data.atom.chirality = CHIRALITY_NONE;
    token->data.atom.chirality_class = 0;
    token->data.atom.atom_class = 0;
    token->data.atom.aromatic = aromatic;

    /* Two-letter elements: Br, Cl */
    if ((c == 'B' && next == 'r') || (c == 'C' && next == 'l')) {
        token->data.atom.symbol[0] = c;
        token->data.atom.symbol[1] = next;
        token->data.atom.symbol[2] = '\0';
        token->length = 2;
        lexer_advance_by(lexer, 2);

        if (c == 'B') {
            token->data.atom.element = ELEM_Br;
        } else {
            token->data.atom.element = ELEM_Cl;
        }
        return CCHEM_OK;
    }

    /* Single letter elements */
    token->data.atom.symbol[0] = c;
    token->data.atom.symbol[1] = '\0';
    token->length = 1;
    lexer_advance(lexer);

    char upper = toupper(c);
    switch (upper) {
        case 'B': token->data.atom.element = ELEM_B; break;
        case 'C': token->data.atom.element = ELEM_C; break;
        case 'N': token->data.atom.element = ELEM_N; break;
        case 'O': token->data.atom.element = ELEM_O; break;
        case 'P': token->data.atom.element = ELEM_P; break;
        case 'S': token->data.atom.element = ELEM_S; break;
        case 'F': token->data.atom.element = ELEM_F; break;
        case 'I': token->data.atom.element = ELEM_I; break;
        default:
            snprintf(lexer->error_msg, sizeof(lexer->error_msg),
                    "Invalid organic atom '%c' at position %zu", c, lexer->pos);
            lexer->has_error = true;
            return CCHEM_ERROR_INVALID_SMILES;
    }

    return CCHEM_OK;
}

/* Parse bracket atom [xxx] */
static cchem_status_t lexer_parse_bracket_atom(lexer_t* lexer, token_t* token) {
    token->type = TOKEN_ATOM_BRACKET;
    token->position = (int)lexer->pos;
    token->data.atom.in_bracket = true;
    token->data.atom.isotope = 0;
    token->data.atom.charge = 0;
    token->data.atom.h_count = -1;
    token->data.atom.chirality = CHIRALITY_NONE;
    token->data.atom.chirality_class = 0;
    token->data.atom.atom_class = 0;
    token->data.atom.aromatic = false;
    token->data.atom.element = ELEM_UNKNOWN;
    token->data.atom.symbol[0] = '\0';

    size_t start_pos = lexer->pos;
    lexer_advance(lexer);  /* Skip '[' */

    /* Parse isotope number (optional) */
    if (isdigit(lexer_current(lexer))) {
        int isotope = 0;
        while (isdigit(lexer_current(lexer))) {
            isotope = isotope * 10 + (lexer_current(lexer) - '0');
            lexer_advance(lexer);
        }
        token->data.atom.isotope = isotope;
    }

    /* Parse element symbol */
    char symbol[4] = {0};
    int sym_len = 0;
    char c = lexer_current(lexer);

    if (c == '*') {
        /* Wildcard atom */
        symbol[0] = '*';
        sym_len = 1;
        token->data.atom.element = ELEM_UNKNOWN;
        lexer_advance(lexer);
    } else if (isupper(c)) {
        /* Normal element - starts with uppercase */
        symbol[sym_len++] = c;
        lexer_advance(lexer);

        /* Second letter if lowercase */
        if (islower(lexer_current(lexer))) {
            symbol[sym_len++] = lexer_current(lexer);
            lexer_advance(lexer);
        }
        symbol[sym_len] = '\0';

        token->data.atom.element = element_from_symbol(symbol);
        if (token->data.atom.element == ELEM_UNKNOWN) {
            snprintf(lexer->error_msg, sizeof(lexer->error_msg),
                    "Unknown element '%s' at position %zu", symbol, start_pos);
            lexer->has_error = true;
            return CCHEM_ERROR_INVALID_SMILES;
        }
    } else if (islower(c)) {
        /* Aromatic element - starts with lowercase */
        token->data.atom.aromatic = true;
        symbol[sym_len++] = c;
        lexer_advance(lexer);

        /* Check for 'se' or 'as' */
        if ((c == 's' && lexer_current(lexer) == 'e') ||
            (c == 'a' && lexer_current(lexer) == 's')) {
            symbol[sym_len++] = lexer_current(lexer);
            lexer_advance(lexer);
        }
        symbol[sym_len] = '\0';

        /* Convert to element */
        bool aromatic;
        token->data.atom.element = parse_element(symbol, &aromatic);
        if (token->data.atom.element == ELEM_UNKNOWN) {
            snprintf(lexer->error_msg, sizeof(lexer->error_msg),
                    "Unknown aromatic element '%s' at position %zu", symbol, start_pos);
            lexer->has_error = true;
            return CCHEM_ERROR_INVALID_SMILES;
        }
    } else {
        snprintf(lexer->error_msg, sizeof(lexer->error_msg),
                "Expected element symbol at position %zu", lexer->pos);
        lexer->has_error = true;
        return CCHEM_ERROR_INVALID_SMILES;
    }

    strncpy(token->data.atom.symbol, symbol, sizeof(token->data.atom.symbol) - 1);

    /* Parse chirality (optional) */
    if (lexer_current(lexer) == '@') {
        lexer_advance(lexer);
        if (lexer_current(lexer) == '@') {
            token->data.atom.chirality = CHIRALITY_CCW;
            lexer_advance(lexer);
        } else if (lexer_current(lexer) == 'T') {
            /* @TH1, @TH2, @TB1-20 */
            lexer_advance(lexer);
            if (lexer_current(lexer) == 'H') {
                lexer_advance(lexer);
                int cls = 0;
                while (isdigit(lexer_current(lexer))) {
                    cls = cls * 10 + (lexer_current(lexer) - '0');
                    lexer_advance(lexer);
                }
                token->data.atom.chirality = (cls == 1) ? CHIRALITY_TH1 : CHIRALITY_TH2;
                token->data.atom.chirality_class = cls;
            } else if (lexer_current(lexer) == 'B') {
                lexer_advance(lexer);
                int cls = 0;
                while (isdigit(lexer_current(lexer))) {
                    cls = cls * 10 + (lexer_current(lexer) - '0');
                    lexer_advance(lexer);
                }
                token->data.atom.chirality = CHIRALITY_TB1;
                token->data.atom.chirality_class = cls;
            }
        } else if (lexer_current(lexer) == 'A' && lexer_peek_at(lexer, 1) == 'L') {
            /* @AL1, @AL2 */
            lexer_advance_by(lexer, 2);
            int cls = 0;
            while (isdigit(lexer_current(lexer))) {
                cls = cls * 10 + (lexer_current(lexer) - '0');
                lexer_advance(lexer);
            }
            token->data.atom.chirality = (cls == 1) ? CHIRALITY_AL1 : CHIRALITY_AL2;
            token->data.atom.chirality_class = cls;
        } else if (lexer_current(lexer) == 'S' && lexer_peek_at(lexer, 1) == 'P') {
            /* @SP1, @SP2, @SP3 */
            lexer_advance_by(lexer, 2);
            int cls = 0;
            while (isdigit(lexer_current(lexer))) {
                cls = cls * 10 + (lexer_current(lexer) - '0');
                lexer_advance(lexer);
            }
            token->data.atom.chirality = CHIRALITY_SP1 + cls - 1;
            token->data.atom.chirality_class = cls;
        } else if (lexer_current(lexer) == 'O' && lexer_peek_at(lexer, 1) == 'H') {
            /* @OH1-30 */
            lexer_advance_by(lexer, 2);
            int cls = 0;
            while (isdigit(lexer_current(lexer))) {
                cls = cls * 10 + (lexer_current(lexer) - '0');
                lexer_advance(lexer);
            }
            token->data.atom.chirality = CHIRALITY_OH1;
            token->data.atom.chirality_class = cls;
        } else {
            token->data.atom.chirality = CHIRALITY_CW;
        }
    }

    /* Parse hydrogen count (optional) */
    if (lexer_current(lexer) == 'H') {
        lexer_advance(lexer);
        if (isdigit(lexer_current(lexer))) {
            int h = 0;
            while (isdigit(lexer_current(lexer))) {
                h = h * 10 + (lexer_current(lexer) - '0');
                lexer_advance(lexer);
            }
            token->data.atom.h_count = h;
        } else {
            token->data.atom.h_count = 1;
        }
    }

    /* Parse charge (optional) */
    if (lexer_current(lexer) == '+' || lexer_current(lexer) == '-') {
        int sign = (lexer_current(lexer) == '+') ? 1 : -1;
        lexer_advance(lexer);

        int charge = 1;
        if (isdigit(lexer_current(lexer))) {
            charge = 0;
            while (isdigit(lexer_current(lexer))) {
                charge = charge * 10 + (lexer_current(lexer) - '0');
                lexer_advance(lexer);
            }
        } else {
            /* Count consecutive +/- signs */
            while (lexer_current(lexer) == (sign > 0 ? '+' : '-')) {
                charge++;
                lexer_advance(lexer);
            }
        }
        token->data.atom.charge = sign * charge;
    }

    /* Parse atom class (optional) :n */
    if (lexer_current(lexer) == ':') {
        lexer_advance(lexer);
        int cls = 0;
        while (isdigit(lexer_current(lexer))) {
            cls = cls * 10 + (lexer_current(lexer) - '0');
            lexer_advance(lexer);
        }
        token->data.atom.atom_class = cls;
    }

    /* Expect closing bracket */
    if (lexer_current(lexer) != ']') {
        snprintf(lexer->error_msg, sizeof(lexer->error_msg),
                "Expected ']' at position %zu, got '%c'", lexer->pos, lexer_current(lexer));
        lexer->has_error = true;
        return CCHEM_ERROR_INVALID_SMILES;
    }
    lexer_advance(lexer);

    token->length = (int)(lexer->pos - start_pos);
    return CCHEM_OK;
}

/* Parse ring closure digit or %nn */
static cchem_status_t lexer_parse_ring_closure(lexer_t* lexer, token_t* token) {
    token->type = TOKEN_RING_CLOSURE;
    token->position = (int)lexer->pos;

    char c = lexer_current(lexer);
    if (c == '%') {
        /* Two-digit ring number %nn */
        lexer_advance(lexer);
        if (!isdigit(lexer_current(lexer)) || !isdigit(lexer_peek_at(lexer, 1))) {
            snprintf(lexer->error_msg, sizeof(lexer->error_msg),
                    "Expected two digits after '%%' at position %zu", lexer->pos);
            lexer->has_error = true;
            return CCHEM_ERROR_INVALID_SMILES;
        }
        int num = (lexer_current(lexer) - '0') * 10;
        lexer_advance(lexer);
        num += (lexer_current(lexer) - '0');
        lexer_advance(lexer);
        token->data.ring_number = num;
        token->length = 3;
    } else {
        /* Single digit */
        token->data.ring_number = c - '0';
        token->length = 1;
        lexer_advance(lexer);
    }

    return CCHEM_OK;
}

cchem_status_t lexer_init(lexer_t* lexer, const char* input) {
    if (!lexer || !input) return CCHEM_ERROR_INVALID_INPUT;

    lexer->input = input;
    lexer->input_len = strlen(input);
    lexer->pos = 0;
    lexer->line = 1;
    lexer->column = 1;
    lexer->has_error = false;
    lexer->error_msg[0] = '\0';

    memset(&lexer->current, 0, sizeof(token_t));

    return CCHEM_OK;
}

void lexer_reset(lexer_t* lexer) {
    if (!lexer) return;

    lexer->pos = 0;
    lexer->column = 1;
    lexer->has_error = false;
    lexer->error_msg[0] = '\0';
}

cchem_status_t lexer_next_token(lexer_t* lexer, token_t* token) {
    if (!lexer || !token) return CCHEM_ERROR_INVALID_INPUT;

    memset(token, 0, sizeof(token_t));

    /* Skip whitespace */
    while (lexer->pos < lexer->input_len && isspace(lexer_current(lexer))) {
        lexer_advance(lexer);
    }

    /* Check for end of input */
    if (lexer->pos >= lexer->input_len) {
        token->type = TOKEN_EOF;
        token->position = (int)lexer->pos;
        token->length = 0;
        return CCHEM_OK;
    }

    char c = lexer_current(lexer);

    /* Branch open/close */
    if (c == '(') {
        token->type = TOKEN_BRANCH_OPEN;
        token->position = (int)lexer->pos;
        token->length = 1;
        lexer_advance(lexer);
        return CCHEM_OK;
    }

    if (c == ')') {
        token->type = TOKEN_BRANCH_CLOSE;
        token->position = (int)lexer->pos;
        token->length = 1;
        lexer_advance(lexer);
        return CCHEM_OK;
    }

    /* Dot (disconnection) */
    if (c == '.') {
        token->type = TOKEN_DOT;
        token->position = (int)lexer->pos;
        token->length = 1;
        lexer_advance(lexer);
        return CCHEM_OK;
    }

    /* Bond symbols */
    if (is_bond_char(c)) {
        token->position = (int)lexer->pos;
        token->length = 1;

        switch (c) {
            case '-': token->type = TOKEN_BOND_SINGLE; break;
            case '=': token->type = TOKEN_BOND_DOUBLE; break;
            case '#': token->type = TOKEN_BOND_TRIPLE; break;
            case '$': token->type = TOKEN_BOND_QUAD; break;
            case ':': token->type = TOKEN_BOND_AROMATIC; break;
            case '/': token->type = TOKEN_BOND_UP; break;
            case '\\': token->type = TOKEN_BOND_DOWN; break;
        }
        token->data.bond_type = bond_from_smiles_char(c);
        lexer_advance(lexer);
        return CCHEM_OK;
    }

    /* Ring closures */
    if (isdigit(c) || c == '%') {
        return lexer_parse_ring_closure(lexer, token);
    }

    /* Bracket atom */
    if (c == '[') {
        return lexer_parse_bracket_atom(lexer, token);
    }

    /* Wildcard */
    if (c == '*') {
        token->type = TOKEN_WILDCARD;
        token->position = (int)lexer->pos;
        token->length = 1;
        token->data.atom.element = ELEM_UNKNOWN;
        token->data.atom.symbol[0] = '*';
        token->data.atom.symbol[1] = '\0';
        lexer_advance(lexer);
        return CCHEM_OK;
    }

    /* Organic subset atom */
    if (is_organic_subset_start(c)) {
        return lexer_parse_organic_atom(lexer, token);
    }

    /* Unknown character */
    snprintf(lexer->error_msg, sizeof(lexer->error_msg),
            "Unexpected character '%c' at position %zu", c, lexer->pos);
    lexer->has_error = true;
    token->type = TOKEN_ERROR;
    token->position = (int)lexer->pos;
    return CCHEM_ERROR_INVALID_SMILES;
}

bool lexer_at_end(const lexer_t* lexer) {
    if (!lexer) return true;
    return lexer->pos >= lexer->input_len;
}

size_t lexer_get_position(const lexer_t* lexer) {
    if (!lexer) return 0;
    return lexer->pos;
}

const char* lexer_get_error(const lexer_t* lexer) {
    if (!lexer) return "NULL lexer";
    return lexer->error_msg;
}
