/**
 * @file test_molecule.c
 * @brief Tests for molecule data structure
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
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

#define ASSERT_NEAR(a, b, eps) do { \
    if (fabs((a) - (b)) > (eps)) { \
        printf("FAILED: %s != %s (line %d)\n", #a, #b, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(x) ASSERT_EQ(!!(x), 1)
#define ASSERT_FALSE(x) ASSERT_EQ(!!(x), 0)
#define ASSERT_NOT_NULL(x) ASSERT_TRUE((x) != NULL)

TEST(create_molecule) {
    molecule_t* mol = molecule_create();
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_atoms, 0);
    ASSERT_EQ(mol->num_bonds, 0);

    molecule_free(mol);
    tests_passed++;
}

TEST(add_atoms) {
    molecule_t* mol = molecule_create();

    int idx1 = molecule_add_atom(mol, ELEM_C);
    int idx2 = molecule_add_atom(mol, ELEM_O);
    int idx3 = molecule_add_atom(mol, ELEM_N);

    ASSERT_EQ(idx1, 0);
    ASSERT_EQ(idx2, 1);
    ASSERT_EQ(idx3, 2);
    ASSERT_EQ(mol->num_atoms, 3);

    ASSERT_EQ(mol->atoms[0].element, ELEM_C);
    ASSERT_EQ(mol->atoms[1].element, ELEM_O);
    ASSERT_EQ(mol->atoms[2].element, ELEM_N);

    molecule_free(mol);
    tests_passed++;
}

TEST(add_bonds) {
    molecule_t* mol = molecule_create();

    molecule_add_atom(mol, ELEM_C);
    molecule_add_atom(mol, ELEM_C);
    molecule_add_atom(mol, ELEM_C);

    int bond1 = molecule_add_bond(mol, 0, 1, BOND_SINGLE);
    int bond2 = molecule_add_bond(mol, 1, 2, BOND_DOUBLE);

    ASSERT_EQ(bond1, 0);
    ASSERT_EQ(bond2, 1);
    ASSERT_EQ(mol->num_bonds, 2);

    ASSERT_EQ(mol->bonds[0].type, BOND_SINGLE);
    ASSERT_EQ(mol->bonds[1].type, BOND_DOUBLE);

    /* Check connectivity */
    ASSERT_EQ(mol->atoms[0].num_neighbors, 1);
    ASSERT_EQ(mol->atoms[1].num_neighbors, 2);
    ASSERT_EQ(mol->atoms[2].num_neighbors, 1);

    molecule_free(mol);
    tests_passed++;
}

TEST(get_bond_between) {
    molecule_t* mol = molecule_create();

    molecule_add_atom(mol, ELEM_C);
    molecule_add_atom(mol, ELEM_O);
    molecule_add_bond(mol, 0, 1, BOND_DOUBLE);

    bond_t* bond = molecule_get_bond_between(mol, 0, 1);
    ASSERT_NOT_NULL(bond);
    ASSERT_EQ(bond->type, BOND_DOUBLE);

    /* Both directions should work */
    bond = molecule_get_bond_between(mol, 1, 0);
    ASSERT_NOT_NULL(bond);

    /* Non-existent bond */
    molecule_add_atom(mol, ELEM_C);
    bond = molecule_get_bond_between(mol, 0, 2);
    ASSERT_TRUE(bond == NULL);

    molecule_free(mol);
    tests_passed++;
}

TEST(clone_molecule) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("c1ccccc1", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    molecule_t* clone = molecule_clone(mol);
    ASSERT_NOT_NULL(clone);

    ASSERT_EQ(clone->num_atoms, mol->num_atoms);
    ASSERT_EQ(clone->num_bonds, mol->num_bonds);

    /* Atoms should match */
    for (int i = 0; i < mol->num_atoms; i++) {
        ASSERT_EQ(clone->atoms[i].element, mol->atoms[i].element);
        ASSERT_EQ(clone->atoms[i].aromatic, mol->atoms[i].aromatic);
    }

    molecule_free(mol);
    molecule_free(clone);
    tests_passed++;
}

TEST(find_fragments) {
    char error_buf[256];

    /* Single fragment */
    molecule_t* mol = smiles_to_molecule("CCCC", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_fragments, 1);
    molecule_free(mol);

    /* Two fragments */
    mol = smiles_to_molecule("CC.CC", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_fragments, 2);
    molecule_free(mol);

    /* Three fragments */
    mol = smiles_to_molecule("C.C.C", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->num_fragments, 3);
    molecule_free(mol);

    tests_passed++;
}

TEST(find_rings) {
    char error_buf[256];

    /* No rings */
    molecule_t* mol = smiles_to_molecule("CCCC", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    molecule_find_rings(mol);
    ASSERT_EQ(mol->num_rings, 0);
    molecule_free(mol);

    /* Single ring (benzene) */
    mol = smiles_to_molecule("c1ccccc1", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    molecule_find_rings(mol);
    ASSERT_EQ(mol->num_rings, 1);
    ASSERT_EQ(mol->rings[0].size, 6);
    molecule_free(mol);

    /* Three-membered ring */
    mol = smiles_to_molecule("C1CC1", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    molecule_find_rings(mol);
    ASSERT_EQ(mol->num_rings, 1);
    ASSERT_EQ(mol->rings[0].size, 3);
    molecule_free(mol);

    tests_passed++;
}

TEST(is_connected) {
    char error_buf[256];

    /* Connected molecule */
    molecule_t* mol = smiles_to_molecule("CCCC", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_TRUE(molecule_is_connected(mol));
    molecule_free(mol);

    /* Disconnected molecule */
    mol = smiles_to_molecule("C.C", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_FALSE(molecule_is_connected(mol));
    molecule_free(mol);

    tests_passed++;
}

TEST(heavy_atom_count) {
    char error_buf[256];

    molecule_t* mol = smiles_to_molecule("CCO", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(molecule_num_heavy_atoms(mol), 3);
    molecule_free(mol);

    /* With explicit hydrogen */
    mol = smiles_to_molecule("[H]C([H])([H])O[H]", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    /* Should have 2 heavy atoms (C and O), rest are H */
    ASSERT_EQ(molecule_num_heavy_atoms(mol), 2);
    molecule_free(mol);

    tests_passed++;
}

TEST(molecular_weight) {
    char error_buf[256];

    /* Methane CH4: 12.01 + 4*1.008 = 16.042 */
    molecule_t* mol = smiles_to_molecule("C", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    molecule_calc_implicit_h(mol);
    double mw = molecule_calc_weight(mol);
    ASSERT_NEAR(mw, 16.04, 0.1);
    molecule_free(mol);

    /* Water H2O: 2*1.008 + 16.00 = 18.016 */
    mol = smiles_to_molecule("O", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    molecule_calc_implicit_h(mol);
    mw = molecule_calc_weight(mol);
    ASSERT_NEAR(mw, 18.0, 0.1);
    molecule_free(mol);

    tests_passed++;
}

TEST(validate_molecule) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("c1ccccc1", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    ASSERT_EQ(molecule_validate(mol), CCHEM_OK);

    molecule_free(mol);
    tests_passed++;
}

TEST(implicit_hydrogen) {
    char error_buf[256];

    /* Methane - carbon should have 4 implicit H */
    molecule_t* mol = smiles_to_molecule("C", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].implicit_h_count, 4);
    molecule_free(mol);

    /* Methanol - OH has 1 implicit H */
    mol = smiles_to_molecule("CO", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[1].implicit_h_count, 1); /* O atom */
    molecule_free(mol);

    /* Ammonia - N has 3 implicit H */
    mol = smiles_to_molecule("N", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].implicit_h_count, 3);
    molecule_free(mol);

    tests_passed++;
}

TEST(atom_neighbors) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("CC(C)C", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    /* Central carbon (index 1) should have 3 neighbors */
    ASSERT_EQ(mol->atoms[1].num_neighbors, 3);
    ASSERT_TRUE(atom_is_bonded_to(&mol->atoms[1], 0));
    ASSERT_TRUE(atom_is_bonded_to(&mol->atoms[1], 2));
    ASSERT_TRUE(atom_is_bonded_to(&mol->atoms[1], 3));

    /* End carbons should have 1 neighbor */
    ASSERT_EQ(mol->atoms[0].num_neighbors, 1);
    ASSERT_EQ(mol->atoms[2].num_neighbors, 1);
    ASSERT_EQ(mol->atoms[3].num_neighbors, 1);

    molecule_free(mol);
    tests_passed++;
}

TEST(bond_order) {
    molecule_t* mol = molecule_create();

    molecule_add_atom(mol, ELEM_C);
    molecule_add_atom(mol, ELEM_C);

    molecule_add_bond(mol, 0, 1, BOND_SINGLE);
    ASSERT_EQ(bond_get_int_order(&mol->bonds[0]), 1);

    mol->bonds[0].type = BOND_DOUBLE;
    ASSERT_EQ(bond_get_int_order(&mol->bonds[0]), 2);

    mol->bonds[0].type = BOND_TRIPLE;
    ASSERT_EQ(bond_get_int_order(&mol->bonds[0]), 3);

    mol->bonds[0].type = BOND_AROMATIC;
    ASSERT_NEAR(bond_get_order(&mol->bonds[0]), 1.5, 0.01);

    molecule_free(mol);
    tests_passed++;
}

int main(void) {
    printf("Running molecule tests...\n\n");

    RUN_TEST(create_molecule);
    RUN_TEST(add_atoms);
    RUN_TEST(add_bonds);
    RUN_TEST(get_bond_between);
    RUN_TEST(clone_molecule);
    RUN_TEST(find_fragments);
    RUN_TEST(find_rings);
    RUN_TEST(is_connected);
    RUN_TEST(heavy_atom_count);
    RUN_TEST(molecular_weight);
    RUN_TEST(validate_molecule);
    RUN_TEST(implicit_hydrogen);
    RUN_TEST(atom_neighbors);
    RUN_TEST(bond_order);

    printf("\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
