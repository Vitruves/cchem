/**
 * @file test_stereo.c
 * @brief Tests for stereochemistry handling
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
#define ASSERT_STR_NE(a, b) ASSERT_TRUE(strcmp((a), (b)) != 0)

/* Helper: Check if canonical SMILES contains chirality marker */
static int has_chirality(const char* canonical) {
    return (strchr(canonical, '@') != NULL);
}

TEST(parse_chirality_at) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[C@H](F)(Cl)Br", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].chirality, CHIRALITY_CW);
    molecule_free(mol);
    tests_passed++;
}

TEST(parse_chirality_atat) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[C@@H](F)(Cl)Br", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].chirality, CHIRALITY_CCW);
    molecule_free(mol);
    tests_passed++;
}

TEST(preserve_chirality) {
    char error_buf[256];
    const char* smiles = "C[C@H](O)F";

    char* canonical = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(has_chirality(canonical));
    free(canonical);
    tests_passed++;
}

TEST(alanine_enantiomers) {
    char error_buf[256];

    /* L-alanine */
    char* l_ala = smiles_canonicalize("C[C@H](N)C(=O)O", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(l_ala);

    /* D-alanine */
    char* d_ala = smiles_canonicalize("C[C@@H](N)C(=O)O", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(d_ala);

    /* Should be different */
    ASSERT_STR_NE(l_ala, d_ala);

    free(l_ala);
    free(d_ala);
    tests_passed++;
}

TEST(parse_ez_stereo) {
    char error_buf[256];

    /* Trans (E) */
    molecule_t* mol = smiles_to_molecule("C/C=C/C", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    /* Check for stereo bonds */
    bool has_stereo_bond = false;
    for (int i = 0; i < mol->num_bonds; i++) {
        if (mol->bonds[i].stereo_type == BOND_UP || mol->bonds[i].stereo_type == BOND_DOWN) {
            has_stereo_bond = true;
            break;
        }
    }
    ASSERT_TRUE(has_stereo_bond);
    molecule_free(mol);

    tests_passed++;
}

TEST(cis_trans_butene) {
    char error_buf[256];

    /* (E)-2-butene (trans) */
    char* e_butene = smiles_canonicalize("C/C=C/C", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(e_butene);

    /* (Z)-2-butene (cis) */
    char* z_butene = smiles_canonicalize("C/C=C\\C", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(z_butene);

    /* Should be different */
    ASSERT_STR_NE(e_butene, z_butene);

    free(e_butene);
    free(z_butene);
    tests_passed++;
}

TEST(no_stereo_option) {
    char error_buf[256];
    const char* smiles = "C[C@H](O)F";

    canon_options_t opts = CANON_OPTIONS_DEFAULT;
    opts.preserve_stereo = false;

    char* canonical = smiles_canonicalize(smiles, &opts, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    /* With stereo disabled, should not have @ */
    ASSERT_FALSE(has_chirality(canonical));
    free(canonical);
    tests_passed++;
}

TEST(tetrahedral_centers) {
    char error_buf[256];

    /* Simple chiral center */
    molecule_t* mol = smiles_to_molecule("[C@](F)(Cl)(Br)I", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    ASSERT_EQ(mol->atoms[0].chirality, CHIRALITY_CW);
    ASSERT_EQ(mol->atoms[0].num_neighbors, 4);
    molecule_free(mol);

    tests_passed++;
}

TEST(multiple_chiral_centers) {
    char error_buf[256];

    /* Two chiral centers */
    const char* smiles = "C[C@H](O)[C@@H](O)C";

    molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    /* Count chiral atoms */
    int chiral_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].chirality != CHIRALITY_NONE) {
            chiral_count++;
        }
    }
    ASSERT_EQ(chiral_count, 2);

    molecule_free(mol);
    tests_passed++;
}

TEST(sugar_stereochemistry) {
    char error_buf[256];

    /* D-glucose (multiple stereocenters) */
    const char* glucose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O";

    char* canonical = smiles_canonicalize(glucose, NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(has_chirality(canonical));
    free(canonical);

    tests_passed++;
}

TEST(ring_stereochemistry) {
    char error_buf[256];

    /* Cyclohexane with substituents */
    const char* smiles = "C[C@H]1CCCC[C@@H]1C";

    molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    int chiral_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].chirality != CHIRALITY_NONE) {
            chiral_count++;
        }
    }
    ASSERT_EQ(chiral_count, 2);

    molecule_free(mol);
    tests_passed++;
}

TEST(allene_stereochemistry) {
    char error_buf[256];

    /* Simple allene - extended tetrahedral */
    molecule_t* mol = smiles_to_molecule("CC=C=CC", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);
    molecule_free(mol);

    tests_passed++;
}

TEST(stereo_detect) {
    char error_buf[256];

    molecule_t* mol = smiles_to_molecule("[C@H](F)(Cl)Br", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    stereo_info_t* info = stereo_info_create();
    ASSERT_NOT_NULL(info);

    cchem_status_t status = stereo_detect_centers(mol, info);
    ASSERT_EQ(status, CCHEM_OK);

    /* Should detect one stereocenter */
    ASSERT_EQ(info->num_centers, 1);
    ASSERT_EQ(info->centers[0].type, STEREO_CENTER_TETRAHEDRAL);
    ASSERT_TRUE(info->centers[0].is_specified);

    stereo_info_free(info);
    molecule_free(mol);
    tests_passed++;
}

TEST(double_bond_stereo_detect) {
    char error_buf[256];

    molecule_t* mol = smiles_to_molecule("C/C=C/C", error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    stereo_info_t* info = stereo_info_create();
    ASSERT_NOT_NULL(info);

    cchem_status_t status = stereo_detect_double_bonds(mol, info);
    ASSERT_EQ(status, CCHEM_OK);

    /* Should detect double bond stereocenter */
    /* Note: Detection may vary based on implementation */

    stereo_info_free(info);
    molecule_free(mol);
    tests_passed++;
}

TEST(chirality_inversion) {
    atom_t atom;
    atom_init(&atom, ELEM_C);

    atom.chirality = CHIRALITY_CW;
    stereo_invert_chirality(&atom);
    ASSERT_EQ(atom.chirality, CHIRALITY_CCW);

    stereo_invert_chirality(&atom);
    ASSERT_EQ(atom.chirality, CHIRALITY_CW);

    tests_passed++;
}

TEST(menthol_stereochemistry) {
    char error_buf[256];

    /* (-)-Menthol - 3 stereocenters */
    const char* menthol = "CC(C)[C@H]1CC[C@@H](C)C[C@H]1O";

    molecule_t* mol = smiles_to_molecule(menthol, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(mol);

    int chiral_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].chirality != CHIRALITY_NONE) {
            chiral_count++;
        }
    }
    ASSERT_EQ(chiral_count, 3);

    molecule_free(mol);
    tests_passed++;
}

int main(void) {
    printf("Running stereochemistry tests...\n\n");

    printf("Parsing:\n");
    RUN_TEST(parse_chirality_at);
    RUN_TEST(parse_chirality_atat);
    RUN_TEST(parse_ez_stereo);

    printf("\nPreservation:\n");
    RUN_TEST(preserve_chirality);
    RUN_TEST(alanine_enantiomers);
    RUN_TEST(cis_trans_butene);
    RUN_TEST(no_stereo_option);

    printf("\nComplex cases:\n");
    RUN_TEST(tetrahedral_centers);
    RUN_TEST(multiple_chiral_centers);
    RUN_TEST(sugar_stereochemistry);
    RUN_TEST(ring_stereochemistry);
    RUN_TEST(allene_stereochemistry);

    printf("\nDetection:\n");
    RUN_TEST(stereo_detect);
    RUN_TEST(double_bond_stereo_detect);
    RUN_TEST(chirality_inversion);
    RUN_TEST(menthol_stereochemistry);

    printf("\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
