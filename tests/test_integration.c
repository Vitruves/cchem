/**
 * @file test_integration.c
 * @brief Integration tests with real molecules (common drugs and compounds)
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

#define ASSERT_TRUE(x) do { \
    if (!(x)) { \
        printf("FAILED: %s (line %d)\n", #x, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_NOT_NULL(x) ASSERT_TRUE((x) != NULL)

/* Helper: Test that SMILES parses and canonicalizes successfully */
static int test_smiles(const char* name, const char* smiles) {
    char error_buf[256];

    molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
    if (!mol) {
        printf("\n    FAILED to parse %s: %s\n", name, error_buf);
        return 0;
    }

    char* canonical = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));
    if (!canonical) {
        printf("\n    FAILED to canonicalize %s: %s\n", name, error_buf);
        molecule_free(mol);
        return 0;
    }

    /* Parse canonical form to verify it's valid */
    molecule_t* mol2 = smiles_to_molecule(canonical, error_buf, sizeof(error_buf));
    if (!mol2) {
        printf("\n    FAILED to re-parse canonical %s: %s\n", name, error_buf);
        molecule_free(mol);
        free(canonical);
        return 0;
    }

    /* Check atom count matches */
    if (mol->num_atoms != mol2->num_atoms) {
        printf("\n    FAILED: atom count mismatch for %s (%d vs %d)\n",
               name, mol->num_atoms, mol2->num_atoms);
        molecule_free(mol);
        molecule_free(mol2);
        free(canonical);
        return 0;
    }

    molecule_free(mol);
    molecule_free(mol2);
    free(canonical);

    return 1;
}

/* Simple molecules */
TEST(water) {
    ASSERT_TRUE(test_smiles("Water", "O"));
    tests_passed++;
}

TEST(methane) {
    ASSERT_TRUE(test_smiles("Methane", "C"));
    tests_passed++;
}

TEST(ethanol) {
    ASSERT_TRUE(test_smiles("Ethanol", "CCO"));
    tests_passed++;
}

TEST(acetic_acid) {
    ASSERT_TRUE(test_smiles("Acetic acid", "CC(=O)O"));
    tests_passed++;
}

TEST(benzene) {
    ASSERT_TRUE(test_smiles("Benzene", "c1ccccc1"));
    tests_passed++;
}

TEST(cyclohexane) {
    ASSERT_TRUE(test_smiles("Cyclohexane", "C1CCCCC1"));
    tests_passed++;
}

/* Common drugs */
TEST(aspirin) {
    ASSERT_TRUE(test_smiles("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"));
    tests_passed++;
}

TEST(caffeine) {
    ASSERT_TRUE(test_smiles("Caffeine", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"));
    tests_passed++;
}

TEST(acetaminophen) {
    ASSERT_TRUE(test_smiles("Acetaminophen/Paracetamol", "CC(=O)Nc1ccc(O)cc1"));
    tests_passed++;
}

TEST(ibuprofen) {
    ASSERT_TRUE(test_smiles("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"));
    tests_passed++;
}

TEST(nicotine) {
    ASSERT_TRUE(test_smiles("Nicotine", "CN1CCCC1c2cccnc2"));
    tests_passed++;
}

TEST(morphine) {
    ASSERT_TRUE(test_smiles("Morphine", "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"));
    tests_passed++;
}

TEST(penicillin_g) {
    ASSERT_TRUE(test_smiles("Penicillin G",
        "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"));
    tests_passed++;
}

/* Amino acids */
TEST(glycine) {
    ASSERT_TRUE(test_smiles("Glycine", "NCC(=O)O"));
    tests_passed++;
}

TEST(alanine) {
    ASSERT_TRUE(test_smiles("Alanine", "CC(N)C(=O)O"));
    tests_passed++;
}

TEST(phenylalanine) {
    ASSERT_TRUE(test_smiles("Phenylalanine", "NC(Cc1ccccc1)C(=O)O"));
    tests_passed++;
}

TEST(tryptophan) {
    ASSERT_TRUE(test_smiles("Tryptophan", "NC(Cc1c[nH]c2ccccc12)C(=O)O"));
    tests_passed++;
}

/* Nucleotides */
TEST(adenine) {
    ASSERT_TRUE(test_smiles("Adenine", "Nc1ncnc2[nH]cnc12"));
    tests_passed++;
}

TEST(guanine) {
    ASSERT_TRUE(test_smiles("Guanine", "Nc1nc2[nH]cnc2c(=O)[nH]1"));
    tests_passed++;
}

TEST(cytosine) {
    ASSERT_TRUE(test_smiles("Cytosine", "Nc1cc[nH]c(=O)n1"));
    tests_passed++;
}

TEST(thymine) {
    ASSERT_TRUE(test_smiles("Thymine", "Cc1c[nH]c(=O)[nH]c1=O"));
    tests_passed++;
}

/* Sugars */
TEST(glucose) {
    ASSERT_TRUE(test_smiles("Glucose", "OCC1OC(O)C(O)C(O)C1O"));
    tests_passed++;
}

TEST(fructose) {
    ASSERT_TRUE(test_smiles("Fructose", "OCC1(O)OCC(O)C(O)C1O"));
    tests_passed++;
}

/* Steroids */
TEST(cholesterol) {
    ASSERT_TRUE(test_smiles("Cholesterol",
        "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"));
    tests_passed++;
}

TEST(testosterone) {
    ASSERT_TRUE(test_smiles("Testosterone",
        "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"));
    tests_passed++;
}

/* Complex molecules */
TEST(taxol) {
    ASSERT_TRUE(test_smiles("Taxol (simplified)",
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"));
    tests_passed++;
}

TEST(chlorophyll_core) {
    /* Porphyrin ring system */
    ASSERT_TRUE(test_smiles("Porphyrin",
        "c1cc2nc1cc3[nH]c(cc4nc(cc5[nH]c(cc2)c5)cc4)cc3"));
    tests_passed++;
}

/* Salts and charged species */
TEST(sodium_chloride) {
    ASSERT_TRUE(test_smiles("Sodium chloride", "[Na+].[Cl-]"));
    tests_passed++;
}

TEST(ammonium) {
    ASSERT_TRUE(test_smiles("Ammonium", "[NH4+]"));
    tests_passed++;
}

TEST(acetate) {
    ASSERT_TRUE(test_smiles("Acetate", "CC(=O)[O-]"));
    tests_passed++;
}

/* Isotopes */
TEST(heavy_water) {
    ASSERT_TRUE(test_smiles("Heavy water (D2O)", "[2H]O[2H]"));
    tests_passed++;
}

TEST(carbon_13_methane) {
    ASSERT_TRUE(test_smiles("13C-Methane", "[13CH4]"));
    tests_passed++;
}

/* Round-trip test */
TEST(round_trip) {
    const char* molecules[] = {
        "C", "CC", "CCC", "CCCC",
        "c1ccccc1", "c1ccc2ccccc2c1",
        "CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O",
        "CN1CCCC1c2cccnc2",
        "[Na+].[Cl-]",
        NULL
    };

    char error_buf[256];

    for (int i = 0; molecules[i] != NULL; i++) {
        char* canon1 = smiles_canonicalize(molecules[i], NULL, error_buf, sizeof(error_buf));
        ASSERT_NOT_NULL(canon1);

        /* Canonicalize the canonical form */
        char* canon2 = smiles_canonicalize(canon1, NULL, error_buf, sizeof(error_buf));
        ASSERT_NOT_NULL(canon2);

        /* Should be identical */
        if (strcmp(canon1, canon2) != 0) {
            printf("\n    Round-trip failed for %s: %s -> %s\n",
                   molecules[i], canon1, canon2);
            free(canon1);
            free(canon2);
            tests_failed++;
            return;
        }

        free(canon1);
        free(canon2);
    }

    tests_passed++;
}

/* Performance test (optional) */
TEST(performance) {
    char error_buf[256];
    const char* smiles = "c1ccc2c(c1)ccc3c2ccc4c3cccc4";  /* Pyrene */

    /* Run many canonicalizations */
    for (int i = 0; i < 100; i++) {
        char* canonical = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));
        ASSERT_NOT_NULL(canonical);
        free(canonical);
    }

    tests_passed++;
}

int main(void) {
    printf("Running integration tests with real molecules...\n\n");

    /* Simple molecules */
    printf("Simple molecules:\n");
    RUN_TEST(water);
    RUN_TEST(methane);
    RUN_TEST(ethanol);
    RUN_TEST(acetic_acid);
    RUN_TEST(benzene);
    RUN_TEST(cyclohexane);

    /* Drugs */
    printf("\nCommon drugs:\n");
    RUN_TEST(aspirin);
    RUN_TEST(caffeine);
    RUN_TEST(acetaminophen);
    RUN_TEST(ibuprofen);
    RUN_TEST(nicotine);
    RUN_TEST(morphine);
    RUN_TEST(penicillin_g);

    /* Amino acids */
    printf("\nAmino acids:\n");
    RUN_TEST(glycine);
    RUN_TEST(alanine);
    RUN_TEST(phenylalanine);
    RUN_TEST(tryptophan);

    /* Nucleotides */
    printf("\nNucleobases:\n");
    RUN_TEST(adenine);
    RUN_TEST(guanine);
    RUN_TEST(cytosine);
    RUN_TEST(thymine);

    /* Sugars */
    printf("\nSugars:\n");
    RUN_TEST(glucose);
    RUN_TEST(fructose);

    /* Steroids */
    printf("\nSteroids:\n");
    RUN_TEST(cholesterol);
    RUN_TEST(testosterone);

    /* Complex molecules */
    printf("\nComplex molecules:\n");
    RUN_TEST(taxol);
    RUN_TEST(chlorophyll_core);

    /* Salts */
    printf("\nSalts and ions:\n");
    RUN_TEST(sodium_chloride);
    RUN_TEST(ammonium);
    RUN_TEST(acetate);

    /* Isotopes */
    printf("\nIsotopes:\n");
    RUN_TEST(heavy_water);
    RUN_TEST(carbon_13_methane);

    /* Round-trip */
    printf("\nRound-trip tests:\n");
    RUN_TEST(round_trip);

    /* Performance */
    printf("\nPerformance:\n");
    RUN_TEST(performance);

    printf("\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
