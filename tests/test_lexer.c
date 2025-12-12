/**
 * @file test_lexer.c
 * @brief Tests for SMILES lexer
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "cchem/cchem.h"

#define TEST(name) static void test_##name(void)
#define RUN_TEST(name) do { \
    printf("  Testing %s... ", #name); \
    test_##name(); \
    printf("OK\n"); \
} while(0)

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_EQ(a, b) do { \
    if ((a) != (b)) { \
        printf("FAILED: %s != %s (line %d)\n", #a, #b, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_STR_EQ(a, b) do { \
    if (strcmp((a), (b)) != 0) { \
        printf("FAILED: %s != %s (line %d)\n", a, b, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(x) ASSERT_EQ(!!(x), 1)
#define ASSERT_FALSE(x) ASSERT_EQ(!!(x), 0)

TEST(empty_input) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, ""), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_EOF);
    tests_passed++;
}

TEST(single_carbon) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "C"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_ORGANIC);
    ASSERT_EQ(token.data.atom.element, ELEM_C);
    ASSERT_FALSE(token.data.atom.aromatic);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_EOF);
    tests_passed++;
}

TEST(aromatic_carbon) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "c"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_ORGANIC);
    ASSERT_EQ(token.data.atom.element, ELEM_C);
    ASSERT_TRUE(token.data.atom.aromatic);
    tests_passed++;
}

TEST(two_letter_elements) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "ClBr"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_ORGANIC);
    ASSERT_EQ(token.data.atom.element, ELEM_Cl);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_ORGANIC);
    ASSERT_EQ(token.data.atom.element, ELEM_Br);
    tests_passed++;
}

TEST(bond_tokens) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "-=#:"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BOND_SINGLE);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BOND_DOUBLE);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BOND_TRIPLE);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BOND_AROMATIC);
    tests_passed++;
}

TEST(stereo_bonds) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "/\\"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BOND_UP);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BOND_DOWN);
    tests_passed++;
}

TEST(branches) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "()"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BRANCH_OPEN);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_BRANCH_CLOSE);
    tests_passed++;
}

TEST(ring_closures) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "1234567890%10%99"), CCHEM_OK);

    token_t token;
    for (int i = 1; i <= 9; i++) {
        ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
        ASSERT_EQ(token.type, TOKEN_RING_CLOSURE);
        ASSERT_EQ(token.data.ring_number, i);
    }

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_RING_CLOSURE);
    ASSERT_EQ(token.data.ring_number, 0);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_RING_CLOSURE);
    ASSERT_EQ(token.data.ring_number, 10);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_RING_CLOSURE);
    ASSERT_EQ(token.data.ring_number, 99);
    tests_passed++;
}

TEST(bracket_atom_simple) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "[Na]"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_BRACKET);
    ASSERT_EQ(token.data.atom.element, ELEM_Na);
    ASSERT_TRUE(token.data.atom.in_bracket);
    tests_passed++;
}

TEST(bracket_atom_isotope) {
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "[13C]"), CCHEM_OK);

    token_t token;
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_BRACKET);
    ASSERT_EQ(token.data.atom.element, ELEM_C);
    ASSERT_EQ(token.data.atom.isotope, 13);
    tests_passed++;
}

TEST(bracket_atom_charge) {
    lexer_t lexer;
    token_t token;

    /* Positive charge */
    ASSERT_EQ(lexer_init(&lexer, "[Na+]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.charge, 1);

    /* Double positive */
    ASSERT_EQ(lexer_init(&lexer, "[Ca++]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.charge, 2);

    /* Numeric charge */
    ASSERT_EQ(lexer_init(&lexer, "[Fe+3]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.charge, 3);

    /* Negative charge */
    ASSERT_EQ(lexer_init(&lexer, "[O-]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.charge, -1);

    /* Double negative */
    ASSERT_EQ(lexer_init(&lexer, "[O--]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.charge, -2);
    tests_passed++;
}

TEST(bracket_atom_hydrogen) {
    lexer_t lexer;
    token_t token;

    /* No H */
    ASSERT_EQ(lexer_init(&lexer, "[CH4]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.h_count, 4);

    /* H with no number = 1 */
    ASSERT_EQ(lexer_init(&lexer, "[NH]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.h_count, 1);
    tests_passed++;
}

TEST(bracket_atom_chirality) {
    lexer_t lexer;
    token_t token;

    /* @ (clockwise) */
    ASSERT_EQ(lexer_init(&lexer, "[C@H]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.chirality, CHIRALITY_CW);

    /* @@ (counter-clockwise) */
    ASSERT_EQ(lexer_init(&lexer, "[C@@H]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.chirality, CHIRALITY_CCW);
    tests_passed++;
}

TEST(bracket_atom_class) {
    lexer_t lexer;
    token_t token;

    ASSERT_EQ(lexer_init(&lexer, "[C:1]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.atom_class, 1);

    ASSERT_EQ(lexer_init(&lexer, "[N:99]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.data.atom.atom_class, 99);
    tests_passed++;
}

TEST(bracket_atom_full) {
    lexer_t lexer;
    token_t token;

    /* Full bracket atom: isotope, element, chirality, H, charge, class */
    ASSERT_EQ(lexer_init(&lexer, "[13C@@H2+:5]"), CCHEM_OK);
    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_BRACKET);
    ASSERT_EQ(token.data.atom.element, ELEM_C);
    ASSERT_EQ(token.data.atom.isotope, 13);
    ASSERT_EQ(token.data.atom.chirality, CHIRALITY_CCW);
    ASSERT_EQ(token.data.atom.h_count, 2);
    ASSERT_EQ(token.data.atom.charge, 1);
    ASSERT_EQ(token.data.atom.atom_class, 5);
    tests_passed++;
}

TEST(dot_disconnection) {
    lexer_t lexer;
    token_t token;

    ASSERT_EQ(lexer_init(&lexer, "C.C"), CCHEM_OK);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_ORGANIC);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_DOT);

    ASSERT_EQ(lexer_next_token(&lexer, &token), CCHEM_OK);
    ASSERT_EQ(token.type, TOKEN_ATOM_ORGANIC);
    tests_passed++;
}

TEST(complex_smiles) {
    /* Benzene */
    lexer_t lexer;
    ASSERT_EQ(lexer_init(&lexer, "c1ccccc1"), CCHEM_OK);

    token_t token;
    int count = 0;
    while (lexer_next_token(&lexer, &token) == CCHEM_OK && token.type != TOKEN_EOF) {
        count++;
    }
    ASSERT_EQ(count, 8); /* 6 carbons + 2 ring closures */
    tests_passed++;
}

int main(void) {
    printf("Running lexer tests...\n\n");

    RUN_TEST(empty_input);
    RUN_TEST(single_carbon);
    RUN_TEST(aromatic_carbon);
    RUN_TEST(two_letter_elements);
    RUN_TEST(bond_tokens);
    RUN_TEST(stereo_bonds);
    RUN_TEST(branches);
    RUN_TEST(ring_closures);
    RUN_TEST(bracket_atom_simple);
    RUN_TEST(bracket_atom_isotope);
    RUN_TEST(bracket_atom_charge);
    RUN_TEST(bracket_atom_hydrogen);
    RUN_TEST(bracket_atom_chirality);
    RUN_TEST(bracket_atom_class);
    RUN_TEST(bracket_atom_full);
    RUN_TEST(dot_disconnection);
    RUN_TEST(complex_smiles);

    printf("\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
