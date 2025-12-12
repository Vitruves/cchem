/**
 * @file test_parser.c
 * @brief Tests for SMILES parser
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "cchem/cchem.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) static void test_##name(void)
#define RUN_TEST(name) do { \
    printf("  Testing %s... ", #name); \
    test_##name(); \
    printf("OK\n"); \
} while(0)

#define ASSERT_EQ(a, b) do { \
    if ((a) != (b)) { \
        printf("FAILED: %s != %s (line %d)\n", #a, #b, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(x) ASSERT_EQ(!!(x), 1)
#define ASSERT_FALSE(x) ASSERT_EQ(!!(x), 0)
#define ASSERT_NOT_NULL(x) ASSERT_TRUE((x) != NULL)
#define ASSERT_NULL(x) ASSERT_TRUE((x) == NULL)

TEST(single_atom) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("C", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 1);
    ASSERT_EQ(mol->num_bonds, 0);
    ASSERT_EQ(mol->atoms[0].element, ELEM_C);

    molecule_free(mol);
    tests_passed++;
}

TEST(simple_chain) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("CCCC", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 4);
    ASSERT_EQ(mol->num_bonds, 3);

    /* Check all atoms are carbon */
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(mol->atoms[i].element, ELEM_C);
    }

    /* Check connectivity */
    ASSERT_EQ(mol->atoms[0].num_neighbors, 1);
    ASSERT_EQ(mol->atoms[1].num_neighbors, 2);
    ASSERT_EQ(mol->atoms[2].num_neighbors, 2);
    ASSERT_EQ(mol->atoms[3].num_neighbors, 1);

    molecule_free(mol);
    tests_passed++;
}

TEST(branched_chain) {
    char error_buf[256];
    /* Isobutane: CC(C)C */
    molecule_t* mol = smiles_to_molecule("CC(C)C", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 4);
    ASSERT_EQ(mol->num_bonds, 3);

    /* Central carbon should have 3 neighbors */
    ASSERT_EQ(mol->atoms[1].num_neighbors, 3);

    molecule_free(mol);
    tests_passed++;
}

TEST(double_bond) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("C=C", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 2);
    ASSERT_EQ(mol->num_bonds, 1);
    ASSERT_EQ(mol->bonds[0].type, BOND_DOUBLE);

    molecule_free(mol);
    tests_passed++;
}

TEST(triple_bond) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("C#C", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 2);
    ASSERT_EQ(mol->num_bonds, 1);
    ASSERT_EQ(mol->bonds[0].type, BOND_TRIPLE);

    molecule_free(mol);
    tests_passed++;
}

TEST(benzene_kekulized) {
    char error_buf[256];
    /* Kekulized benzene */
    molecule_t* mol = smiles_to_molecule("C1=CC=CC=C1", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 6);
    ASSERT_EQ(mol->num_bonds, 6);

    /* All atoms should be in ring */
    for (int i = 0; i < 6; i++) {
        ASSERT_TRUE(mol->atoms[i].ring_count > 0);
    }

    molecule_free(mol);
    tests_passed++;
}

TEST(benzene_aromatic) {
    char error_buf[256];
    /* Aromatic benzene */
    molecule_t* mol = smiles_to_molecule("c1ccccc1", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 6);
    ASSERT_EQ(mol->num_bonds, 6);

    /* All atoms should be aromatic */
    for (int i = 0; i < 6; i++) {
        ASSERT_TRUE(mol->atoms[i].aromatic);
    }

    molecule_free(mol);
    tests_passed++;
}

TEST(charged_atoms) {
    char error_buf[256];
    /* Sodium cation */
    molecule_t* mol = smiles_to_molecule("[Na+]", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 1);
    ASSERT_EQ(mol->atoms[0].element, ELEM_Na);
    ASSERT_EQ(mol->atoms[0].charge, 1);

    molecule_free(mol);

    /* Chloride anion */
    mol = smiles_to_molecule("[Cl-]", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].element, ELEM_Cl);
    ASSERT_EQ(mol->atoms[0].charge, -1);

    molecule_free(mol);
    tests_passed++;
}

TEST(disconnected_fragments) {
    char error_buf[256];
    /* NaCl */
    molecule_t* mol = smiles_to_molecule("[Na+].[Cl-]", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 2);
    ASSERT_EQ(mol->num_bonds, 0);
    ASSERT_EQ(mol->num_fragments, 2);

    molecule_free(mol);
    tests_passed++;
}

TEST(cyclopropane) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("C1CC1", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 3);
    ASSERT_EQ(mol->num_bonds, 3);

    /* All atoms in ring */
    for (int i = 0; i < 3; i++) {
        ASSERT_TRUE(mol->atoms[i].ring_count > 0);
    }

    molecule_free(mol);
    tests_passed++;
}

TEST(naphthalene) {
    char error_buf[256];
    /* Naphthalene - fused rings */
    molecule_t* mol = smiles_to_molecule("c1ccc2ccccc2c1", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 10);
    ASSERT_EQ(mol->num_bonds, 11);

    molecule_free(mol);
    tests_passed++;
}

TEST(isotope) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[13C]", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].isotope, 13);

    molecule_free(mol);
    tests_passed++;
}

TEST(explicit_hydrogen) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[CH4]", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].h_count, 4);

    molecule_free(mol);
    tests_passed++;
}

TEST(chirality) {
    char error_buf[256];
    /* L-alanine */
    molecule_t* mol = smiles_to_molecule("C[C@H](N)C(=O)O", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    /* Check that chirality is preserved */
    ASSERT_EQ(mol->atoms[1].chirality, CHIRALITY_CW);

    molecule_free(mol);
    tests_passed++;
}

TEST(invalid_smiles) {
    char error_buf[256];

    /* Unclosed ring */
    molecule_t* mol = smiles_to_molecule("C1CC", error_buf, sizeof(error_buf));
    ASSERT_NULL(mol);

    /* Unclosed branch */
    mol = smiles_to_molecule("C(C", error_buf, sizeof(error_buf));
    ASSERT_NULL(mol);

    /* Invalid character */
    mol = smiles_to_molecule("CXC", error_buf, sizeof(error_buf));
    ASSERT_NULL(mol);

    tests_passed++;
}

TEST(ethanol) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("CCO", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 3);
    ASSERT_EQ(mol->num_bonds, 2);

    ASSERT_EQ(mol->atoms[0].element, ELEM_C);
    ASSERT_EQ(mol->atoms[1].element, ELEM_C);
    ASSERT_EQ(mol->atoms[2].element, ELEM_O);

    molecule_free(mol);
    tests_passed++;
}

TEST(acetic_acid) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("CC(=O)O", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 4);
    ASSERT_EQ(mol->num_bonds, 3);

    /* Check for double bond */
    bool has_double = false;
    for (int i = 0; i < mol->num_bonds; i++) {
        if (mol->bonds[i].type == BOND_DOUBLE) {
            has_double = true;
            break;
        }
    }
    ASSERT_TRUE(has_double);

    molecule_free(mol);
    tests_passed++;
}

TEST(multiple_rings) {
    char error_buf[256];
    /* Bicyclo compound */
    molecule_t* mol = smiles_to_molecule("C1CC2CCC1C2", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 7);

    molecule_free(mol);
    tests_passed++;
}

TEST(stereo_double_bond) {
    char error_buf[256];
    /* (E)-2-butene */
    molecule_t* mol = smiles_to_molecule("C/C=C/C", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 4);
    ASSERT_EQ(mol->num_bonds, 3);

    molecule_free(mol);
    tests_passed++;
}

TEST(halogens) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("CBrClFI", error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 5);
    ASSERT_EQ(mol->atoms[1].element, ELEM_Br);
    ASSERT_EQ(mol->atoms[2].element, ELEM_Cl);
    ASSERT_EQ(mol->atoms[3].element, ELEM_F);
    ASSERT_EQ(mol->atoms[4].element, ELEM_I);

    molecule_free(mol);
    tests_passed++;
}

int main(void) {
    printf("Running parser tests...\n\n");

    RUN_TEST(single_atom);
    RUN_TEST(simple_chain);
    RUN_TEST(branched_chain);
    RUN_TEST(double_bond);
    RUN_TEST(triple_bond);
    RUN_TEST(benzene_kekulized);
    RUN_TEST(benzene_aromatic);
    RUN_TEST(charged_atoms);
    RUN_TEST(disconnected_fragments);
    RUN_TEST(cyclopropane);
    RUN_TEST(naphthalene);
    RUN_TEST(isotope);
    RUN_TEST(explicit_hydrogen);
    RUN_TEST(chirality);
    RUN_TEST(invalid_smiles);
    RUN_TEST(ethanol);
    RUN_TEST(acetic_acid);
    RUN_TEST(multiple_rings);
    RUN_TEST(stereo_double_bond);
    RUN_TEST(halogens);

    printf("\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
