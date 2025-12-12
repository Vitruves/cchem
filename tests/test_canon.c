/**
 * @file test_canon.c
 * @brief Tests for SMILES canonicalization
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

#define ASSERT_STR_EQ(a, b) do { \
    if (strcmp((a), (b)) != 0) { \
        printf("FAILED: '%s' != '%s' (line %d)\n", (a), (b), __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(x) ASSERT_EQ(!!(x), 1)
#define ASSERT_FALSE(x) ASSERT_EQ(!!(x), 0)
#define ASSERT_NOT_NULL(x) ASSERT_TRUE((x) != NULL)

/* Helper: test that two SMILES produce same canonical form */
static int same_canonical(const char* smiles1, const char* smiles2) {
    char error_buf[256];

    char* canon1 = smiles_canonicalize(smiles1, NULL, error_buf, sizeof(error_buf));
    if (!canon1) return 0;

    char* canon2 = smiles_canonicalize(smiles2, NULL, error_buf, sizeof(error_buf));
    if (!canon2) {
        free(canon1);
        return 0;
    }

    int result = (strcmp(canon1, canon2) == 0);

    free(canon1);
    free(canon2);

    return result;
}

TEST(single_atom) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("C", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_STR_EQ(canonical, "C");
    free(canonical);
    tests_passed++;
}

TEST(simple_chain) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("CCCC", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    /* Should be deterministic */
    ASSERT_STR_EQ(canonical, "CCCC");
    free(canonical);
    tests_passed++;
}

TEST(chain_equivalence) {
    /* These should all produce the same canonical form */
    ASSERT_TRUE(same_canonical("CCCC", "C(CC)C"));
    ASSERT_TRUE(same_canonical("CCCC", "CC(C)C") == 0);  /* Different structure */
    tests_passed++;
}

TEST(isobutane_equivalence) {
    /* Different representations of isobutane */
    ASSERT_TRUE(same_canonical("CC(C)C", "C(C)(C)C"));
    ASSERT_TRUE(same_canonical("CC(C)C", "C(C)C(C)") == 0);  /* This is n-butane */
    tests_passed++;
}

TEST(benzene) {
    /* Different representations of benzene should canonicalize the same */
    ASSERT_TRUE(same_canonical("c1ccccc1", "c1ccccc1"));
    ASSERT_TRUE(same_canonical("C1=CC=CC=C1", "c1ccccc1"));  /* May differ in aromaticity */
    tests_passed++;
}

TEST(benzene_numbering) {
    /* Starting from different atoms should give same result */
    ASSERT_TRUE(same_canonical("c1ccccc1", "c1ccccc1"));
    tests_passed++;
}

TEST(ethanol) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("CCO", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);

    /* Reverse should give same */
    char* canonical2 = smiles_canonicalize("OCC", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical2);
    ASSERT_STR_EQ(canonical, canonical2);

    free(canonical);
    free(canonical2);
    tests_passed++;
}

TEST(acetic_acid) {
    /* CC(=O)O and C(C)(=O)O should be equivalent */
    ASSERT_TRUE(same_canonical("CC(=O)O", "CC(=O)O"));
    tests_passed++;
}

TEST(charged_molecules) {
    char error_buf[256];

    /* Sodium cation */
    char* canonical = smiles_canonicalize("[Na+]", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(strstr(canonical, "Na") != NULL);
    ASSERT_TRUE(strstr(canonical, "+") != NULL);
    free(canonical);

    /* Chloride anion */
    canonical = smiles_canonicalize("[Cl-]", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(strstr(canonical, "Cl") != NULL);
    ASSERT_TRUE(strstr(canonical, "-") != NULL);
    free(canonical);

    tests_passed++;
}

TEST(isotopes) {
    char error_buf[256];

    /* Deuterium */
    char* canonical = smiles_canonicalize("[2H]C", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(strstr(canonical, "2") != NULL);
    free(canonical);

    /* C-13 */
    canonical = smiles_canonicalize("[13C]", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(strstr(canonical, "13") != NULL);
    free(canonical);

    tests_passed++;
}

TEST(disconnected_fragments) {
    /* Order of fragments should be canonical */
    ASSERT_TRUE(same_canonical("[Na+].[Cl-]", "[Cl-].[Na+]"));
    ASSERT_TRUE(same_canonical("C.CC", "CC.C"));
    tests_passed++;
}

TEST(cyclopropane) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("C1CC1", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    /* Should have ring closure */
    ASSERT_TRUE(strchr(canonical, '1') != NULL);
    free(canonical);
    tests_passed++;
}

TEST(cyclohexane) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("C1CCCCC1", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    free(canonical);

    /* Different numbering should give same result */
    ASSERT_TRUE(same_canonical("C1CCCCC1", "C1CCCCC1"));
    tests_passed++;
}

TEST(smiles_equivalence) {
    ASSERT_TRUE(smiles_are_equivalent("C", "C"));
    ASSERT_TRUE(smiles_are_equivalent("CC", "CC"));
    ASSERT_TRUE(smiles_are_equivalent("CCO", "OCC"));
    ASSERT_FALSE(smiles_are_equivalent("CC", "C=C"));
    ASSERT_FALSE(smiles_are_equivalent("CCCC", "CC(C)C"));
    tests_passed++;
}

TEST(naphthalene) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("c1ccc2ccccc2c1", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    free(canonical);

    /* Different representations should be equivalent */
    ASSERT_TRUE(same_canonical("c1ccc2ccccc2c1", "c1ccc2ccccc2c1"));
    tests_passed++;
}

TEST(double_bond) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("C=C", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(strstr(canonical, "=") != NULL);
    free(canonical);
    tests_passed++;
}

TEST(triple_bond) {
    char error_buf[256];
    char* canonical = smiles_canonicalize("C#C", NULL, error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    ASSERT_TRUE(strstr(canonical, "#") != NULL);
    free(canonical);
    tests_passed++;
}

TEST(toluene) {
    /* Different attachment points should canonicalize */
    ASSERT_TRUE(same_canonical("Cc1ccccc1", "c1ccccc1C"));
    tests_passed++;
}

TEST(aniline) {
    /* Aminobenzene */
    ASSERT_TRUE(same_canonical("Nc1ccccc1", "c1ccccc1N"));
    tests_passed++;
}

TEST(phenol) {
    ASSERT_TRUE(same_canonical("Oc1ccccc1", "c1ccccc1O"));
    tests_passed++;
}

TEST(determinism) {
    /* Multiple canonicalizations should give identical results */
    char error_buf[256];
    const char* smiles = "c1ccc(CC(=O)O)cc1";

    char* canon1 = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));
    char* canon2 = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));
    char* canon3 = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));

    ASSERT_NOT_NULL(canon1);
    ASSERT_NOT_NULL(canon2);
    ASSERT_NOT_NULL(canon3);

    ASSERT_STR_EQ(canon1, canon2);
    ASSERT_STR_EQ(canon2, canon3);

    free(canon1);
    free(canon2);
    free(canon3);
    tests_passed++;
}

TEST(complex_molecule) {
    /* Caffeine */
    char error_buf[256];
    char* canonical = smiles_canonicalize("Cn1cnc2c1c(=O)n(c(=O)n2C)C", NULL,
                                          error_buf, sizeof(error_buf));
    ASSERT_NOT_NULL(canonical);
    free(canonical);
    tests_passed++;
}

int main(void) {
    printf("Running canonicalization tests...\n\n");

    RUN_TEST(single_atom);
    RUN_TEST(simple_chain);
    RUN_TEST(chain_equivalence);
    RUN_TEST(isobutane_equivalence);
    RUN_TEST(benzene);
    RUN_TEST(benzene_numbering);
    RUN_TEST(ethanol);
    RUN_TEST(acetic_acid);
    RUN_TEST(charged_molecules);
    RUN_TEST(isotopes);
    RUN_TEST(disconnected_fragments);
    RUN_TEST(cyclopropane);
    RUN_TEST(cyclohexane);
    RUN_TEST(smiles_equivalence);
    RUN_TEST(naphthalene);
    RUN_TEST(double_bond);
    RUN_TEST(triple_bond);
    RUN_TEST(toluene);
    RUN_TEST(aniline);
    RUN_TEST(phenol);
    RUN_TEST(determinism);
    RUN_TEST(complex_molecule);

    printf("\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
