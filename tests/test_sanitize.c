/**
 * @file test_sanitize.c
 * @brief Comprehensive tests for molecular sanitization functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cchem/cchem.h"

/* Test framework macros */
static int test_count = 0;
static int pass_count = 0;
static int fail_count = 0;

#define TEST(name) \
    static void test_##name(void); \
    static void run_test_##name(void) { \
        test_count++; \
        printf("  [%d] %s... ", test_count, #name); \
        fflush(stdout); \
        test_##name(); \
    } \
    static void test_##name(void)

#define ASSERT(condition) do { \
    if (!(condition)) { \
        printf("FAIL\n    Assertion failed: %s\n    at %s:%d\n", \
               #condition, __FILE__, __LINE__); \
        fail_count++; \
        return; \
    } \
} while(0)

#define ASSERT_STR_EQ(expected, actual) do { \
    if (strcmp(expected, actual) != 0) { \
        printf("FAIL\n    Expected: %s\n    Actual:   %s\n    at %s:%d\n", \
               expected, actual, __FILE__, __LINE__); \
        fail_count++; \
        return; \
    } \
} while(0)

#define PASS() do { \
    printf("PASS\n"); \
    pass_count++; \
} while(0)

#define RUN_TEST(name) run_test_##name()

/* =========================================================================
 * Unsalt Tests
 * ========================================================================= */

TEST(unsalt_remove_sodium_chloride) {
    char error_buf[256];
    char* result = smiles_sanitize("[Na+].[Cl-].c1ccccc1", SANITIZE_UNSALT,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    ASSERT_STR_EQ("c1ccccc1", result);
    free(result);
    PASS();
}

TEST(unsalt_remove_sodium_from_acetate) {
    char error_buf[256];
    /* Sodium acetate -> acetic acid (after unsalt only, no neutralize) */
    molecule_t* mol = smiles_to_molecule("[Na+].CC(=O)[O-]", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    sanitize_options_t opts = SANITIZE_OPTIONS_DEFAULT;
    opts.flags = SANITIZE_UNSALT;

    cchem_status_t status = molecule_sanitize(mol, &opts, error_buf, sizeof(error_buf));
    ASSERT(status == CCHEM_OK);

    /* After unsalt, we should have acetate anion (still charged) */
    ASSERT(mol->num_atoms > 0);
    molecule_free(mol);
    PASS();
}

TEST(unsalt_keep_largest_organic_fragment) {
    char error_buf[256];
    /* Two organic fragments - keep the larger one */
    char* result = smiles_sanitize("[Na+].C.CCCCC", SANITIZE_UNSALT,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Result should be pentane (5 carbons) - canonical form may vary */
    ASSERT(strlen(result) >= 5);  /* At least 5 characters for CCCCC */
    ASSERT(strchr(result, 'N') == NULL);  /* No sodium */
    /* Check it's pentane by verifying it has 5 carbons */
    int c_count = 0;
    for (int i = 0; result[i]; i++) {
        if (result[i] == 'C') c_count++;
    }
    ASSERT(c_count == 5);
    free(result);
    PASS();
}

TEST(unsalt_single_fragment_unchanged) {
    char error_buf[256];
    char* result = smiles_sanitize("c1ccccc1", SANITIZE_UNSALT,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    ASSERT_STR_EQ("c1ccccc1", result);
    free(result);
    PASS();
}

/* =========================================================================
 * Aromatize Tests
 * ========================================================================= */

TEST(aromatize_benzene_kekule_to_aromatic) {
    char error_buf[256];
    char* result = smiles_sanitize("C1=CC=CC=C1", SANITIZE_AROMATIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    ASSERT_STR_EQ("c1ccccc1", result);
    free(result);
    PASS();
}

TEST(aromatize_pyridine) {
    char error_buf[256];
    char* result = smiles_sanitize("C1=CC=NC=C1", SANITIZE_AROMATIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should contain lowercase aromatic atoms */
    ASSERT(strchr(result, 'c') != NULL);
    free(result);
    PASS();
}

TEST(aromatize_naphthalene) {
    char error_buf[256];
    char* result = smiles_sanitize("C1=CC2=CC=CC=C2C=C1", SANITIZE_AROMATIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* All carbons should be aromatic */
    for (int i = 0; result[i]; i++) {
        if (result[i] == 'C') {
            /* Uppercase C shouldn't appear in aromatic naphthalene */
            ASSERT(0);
        }
    }
    free(result);
    PASS();
}

TEST(aromatize_already_aromatic_unchanged) {
    char error_buf[256];
    char* result = smiles_sanitize("c1ccccc1", SANITIZE_AROMATIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    ASSERT_STR_EQ("c1ccccc1", result);
    free(result);
    PASS();
}

/* =========================================================================
 * Kekulize Tests
 * ========================================================================= */

TEST(kekulize_benzene_aromatic_to_kekule) {
    char error_buf[256];
    char* result = smiles_sanitize("c1ccccc1", SANITIZE_KEKULIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should have alternating double/single bonds (uppercase C and = signs) */
    ASSERT(strchr(result, 'C') != NULL);
    ASSERT(strchr(result, '=') != NULL);
    free(result);
    PASS();
}

TEST(kekulize_pyridine) {
    char error_buf[256];
    char* result = smiles_sanitize("c1ccncc1", SANITIZE_KEKULIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should have uppercase C and N with = signs */
    ASSERT(strchr(result, 'C') != NULL || strchr(result, 'N') != NULL);
    free(result);
    PASS();
}

/* =========================================================================
 * Neutralize Tests
 * ========================================================================= */

TEST(neutralize_phenolate_to_phenol) {
    char error_buf[256];
    char* result = smiles_sanitize("c1ccc([O-])cc1", SANITIZE_NEUTRALIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    ASSERT_STR_EQ("c1ccccc1O", result);
    free(result);
    PASS();
}

TEST(neutralize_ammonium_to_amine) {
    char error_buf[256];
    char* result = smiles_sanitize("CC[NH3+]", SANITIZE_NEUTRALIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have + charge */
    ASSERT(strchr(result, '+') == NULL);
    free(result);
    PASS();
}

TEST(neutralize_carboxylate) {
    char error_buf[256];
    char* result = smiles_sanitize("CC(=O)[O-]", SANITIZE_NEUTRALIZE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have - charge */
    ASSERT(strchr(result, '-') == NULL);
    free(result);
    PASS();
}

TEST(neutralize_preserve_quaternary_nitrogen) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("C[N+](C)(C)C", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    sanitize_options_t opts = SANITIZE_OPTIONS_DEFAULT;
    opts.flags = SANITIZE_NEUTRALIZE;
    opts.preserve_quaternary_n = true;

    cchem_status_t status = molecule_sanitize(mol, &opts, error_buf, sizeof(error_buf));
    ASSERT(status == CCHEM_OK);

    /* Quaternary N should still have positive charge */
    bool found_positive_n = false;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_N && mol->atoms[i].charge > 0) {
            found_positive_n = true;
            break;
        }
    }
    ASSERT(found_positive_n);
    molecule_free(mol);
    PASS();
}

/* =========================================================================
 * Remove Stereo Tests
 * ========================================================================= */

TEST(remove_stereo_ez_double_bond) {
    char error_buf[256];
    char* result = smiles_sanitize("C/C=C/C", SANITIZE_REMOVE_STEREO,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have / or \ */
    ASSERT(strchr(result, '/') == NULL);
    ASSERT(strchr(result, '\\') == NULL);
    free(result);
    PASS();
}

TEST(remove_stereo_chiral_center) {
    char error_buf[256];
    char* result = smiles_sanitize("[C@H](F)(Cl)Br", SANITIZE_REMOVE_STEREO,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have @ */
    ASSERT(strchr(result, '@') == NULL);
    free(result);
    PASS();
}

TEST(remove_stereo_complex_molecule) {
    char error_buf[256];
    /* Molecule with both types of stereochemistry */
    char* result = smiles_sanitize("[C@@H](O)(F)C/C=C/C",
                                    SANITIZE_REMOVE_STEREO,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    ASSERT(strchr(result, '@') == NULL);
    ASSERT(strchr(result, '/') == NULL);
    ASSERT(strchr(result, '\\') == NULL);
    free(result);
    PASS();
}

/* =========================================================================
 * Remove Isotopes Tests
 * ========================================================================= */

TEST(remove_isotopes_carbon_13) {
    char error_buf[256];
    char* result = smiles_sanitize("[13C]CO", SANITIZE_REMOVE_ISOTOPES,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have 13 isotope label */
    ASSERT(strstr(result, "13") == NULL);
    free(result);
    PASS();
}

TEST(remove_isotopes_deuterium) {
    char error_buf[256];
    char* result = smiles_sanitize("[2H]C([2H])([2H])C", SANITIZE_REMOVE_ISOTOPES,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have 2H */
    ASSERT(strstr(result, "2H") == NULL);
    free(result);
    PASS();
}

TEST(remove_isotopes_multiple) {
    char error_buf[256];
    char* result = smiles_sanitize("[13C][14C][15N]", SANITIZE_REMOVE_ISOTOPES,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have any isotope labels */
    ASSERT(strstr(result, "13") == NULL);
    ASSERT(strstr(result, "14") == NULL);
    ASSERT(strstr(result, "15") == NULL);
    free(result);
    PASS();
}

/* =========================================================================
 * Remove H Tests
 * ========================================================================= */

TEST(remove_h_explicit_hydrogens) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[H]C([H])([H])O[H]", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    int initial_h_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) initial_h_count++;
    }
    ASSERT(initial_h_count > 0);

    sanitize_options_t opts = SANITIZE_OPTIONS_DEFAULT;
    opts.flags = SANITIZE_REMOVE_H;
    opts.remove_all_h = true;

    cchem_status_t status = molecule_sanitize(mol, &opts, error_buf, sizeof(error_buf));
    ASSERT(status == CCHEM_OK);

    int final_h_count = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_H) final_h_count++;
    }
    ASSERT(final_h_count < initial_h_count);

    molecule_free(mol);
    PASS();
}

/* =========================================================================
 * Complete Sanitization Tests
 * ========================================================================= */

TEST(complete_sanitize_drug_salt) {
    char error_buf[256];
    /* Aspirin sodium salt */
    char* result = smiles_sanitize("[Na+].CC(=O)Oc1ccccc1C(=O)[O-]",
                                    SANITIZE_COMPLETE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have Na+ or charges */
    ASSERT(strstr(result, "Na") == NULL);
    ASSERT(strchr(result, '+') == NULL);
    ASSERT(strchr(result, '-') == NULL);
    /* Should be aromatic */
    ASSERT(strchr(result, 'c') != NULL);
    free(result);
    PASS();
}

TEST(complete_sanitize_potassium_phenolate) {
    char error_buf[256];
    char* result = smiles_sanitize("[K+].c1ccc([O-])cc1",
                                    SANITIZE_COMPLETE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    ASSERT_STR_EQ("c1ccccc1O", result);
    free(result);
    PASS();
}

TEST(complete_sanitize_multiple_salts) {
    char error_buf[256];
    char* result = smiles_sanitize("[Na+].[Na+].[O-]C(=O)CCC(=O)[O-]",
                                    SANITIZE_COMPLETE,
                                    error_buf, sizeof(error_buf));
    ASSERT(result != NULL);
    /* Should not have Na+ */
    ASSERT(strstr(result, "Na") == NULL);
    /* Should not have charges */
    ASSERT(strchr(result, '+') == NULL);
    ASSERT(strchr(result, '-') == NULL);
    free(result);
    PASS();
}

/* =========================================================================
 * Fragment Functions Tests
 * ========================================================================= */

TEST(fragment_info_single_fragment) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("c1ccccc1", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    int num_frags, largest_idx, largest_atoms;
    cchem_status_t status = molecule_get_fragment_info(mol, &num_frags,
                                                        &largest_idx, &largest_atoms);
    ASSERT(status == CCHEM_OK);
    ASSERT(num_frags == 1);
    ASSERT(largest_idx == 0);
    ASSERT(largest_atoms == 6);

    molecule_free(mol);
    PASS();
}

TEST(fragment_info_multiple_fragments) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[Na+].[Cl-].c1ccccc1", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    int num_frags, largest_idx, largest_atoms;
    cchem_status_t status = molecule_get_fragment_info(mol, &num_frags,
                                                        &largest_idx, &largest_atoms);
    ASSERT(status == CCHEM_OK);
    ASSERT(num_frags == 3);
    ASSERT(largest_atoms == 6);  /* Benzene has 6 heavy atoms */

    molecule_free(mol);
    PASS();
}

TEST(fragment_is_salt_sodium) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[Na+].c1ccccc1", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    molecule_find_fragments(mol);

    /* Find the sodium fragment */
    bool found_salt = false;
    for (int f = 0; f < mol->num_fragments; f++) {
        if (molecule_fragment_is_salt(mol, f)) {
            found_salt = true;
            break;
        }
    }
    ASSERT(found_salt);

    molecule_free(mol);
    PASS();
}

TEST(fragment_is_organic) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("[Na+].c1ccccc1", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    molecule_find_fragments(mol);

    /* Find the organic fragment */
    bool found_organic = false;
    for (int f = 0; f < mol->num_fragments; f++) {
        if (molecule_fragment_is_organic(mol, f)) {
            found_organic = true;
            break;
        }
    }
    ASSERT(found_organic);

    molecule_free(mol);
    PASS();
}

/* =========================================================================
 * Tautomer Tests
 * ========================================================================= */

TEST(tautomer_enumerate_basic) {
    char error_buf[256];
    molecule_t* mol = smiles_to_molecule("CC(=O)N", error_buf, sizeof(error_buf));
    ASSERT(mol != NULL);

    tautomer_result_t result;
    tautomer_options_t opts = TAUTOMER_OPTIONS_DEFAULT;

    cchem_status_t status = tautomer_enumerate(mol, &opts, &result);
    ASSERT(status == CCHEM_OK);
    ASSERT(result.num_tautomers >= 1);
    ASSERT(result.smiles != NULL);
    ASSERT(result.smiles[0] != NULL);

    tautomer_result_free(&result);
    molecule_free(mol);
    PASS();
}

TEST(tautomer_canonical) {
    char error_buf[256];
    char* canonical = tautomer_canonical("CC(=O)N", NULL,
                                          error_buf, sizeof(error_buf));
    ASSERT(canonical != NULL);
    /* Should return a valid SMILES */
    ASSERT(strlen(canonical) > 0);
    free(canonical);
    PASS();
}

/* =========================================================================
 * Flag Parsing Tests
 * ========================================================================= */

TEST(parse_flags_single) {
    sanitize_flags_t flags;
    cchem_status_t status = sanitize_parse_flags("unsalt", &flags);
    ASSERT(status == CCHEM_OK);
    ASSERT(flags == SANITIZE_UNSALT);
    PASS();
}

TEST(parse_flags_multiple) {
    sanitize_flags_t flags;
    cchem_status_t status = sanitize_parse_flags("unsalt,aromatize,neutralize", &flags);
    ASSERT(status == CCHEM_OK);
    ASSERT(flags & SANITIZE_UNSALT);
    ASSERT(flags & SANITIZE_AROMATIZE);
    ASSERT(flags & SANITIZE_NEUTRALIZE);
    PASS();
}

TEST(parse_flags_preset_complete) {
    sanitize_flags_t flags;
    cchem_status_t status = sanitize_parse_flags("complete", &flags);
    ASSERT(status == CCHEM_OK);
    ASSERT(flags == SANITIZE_COMPLETE);
    PASS();
}

TEST(parse_flags_preset_all) {
    sanitize_flags_t flags;
    cchem_status_t status = sanitize_parse_flags("all", &flags);
    ASSERT(status == CCHEM_OK);
    ASSERT(flags == SANITIZE_ALL);
    PASS();
}

TEST(parse_flags_invalid) {
    sanitize_flags_t flags;
    cchem_status_t status = sanitize_parse_flags("invalid_flag", &flags);
    ASSERT(status != CCHEM_OK);
    PASS();
}

TEST(flags_to_string) {
    char buf[256];
    sanitize_flags_to_string(SANITIZE_UNSALT | SANITIZE_AROMATIZE, buf, sizeof(buf));
    ASSERT(strstr(buf, "unsalt") != NULL);
    ASSERT(strstr(buf, "aromatize") != NULL);
    PASS();
}

/* =========================================================================
 * Main
 * ========================================================================= */

int main(void) {
    printf("\n=== Sanitization Tests ===\n\n");

    printf("Unsalt Tests:\n");
    RUN_TEST(unsalt_remove_sodium_chloride);
    RUN_TEST(unsalt_remove_sodium_from_acetate);
    RUN_TEST(unsalt_keep_largest_organic_fragment);
    RUN_TEST(unsalt_single_fragment_unchanged);

    printf("\nAromatize Tests:\n");
    RUN_TEST(aromatize_benzene_kekule_to_aromatic);
    RUN_TEST(aromatize_pyridine);
    RUN_TEST(aromatize_naphthalene);
    RUN_TEST(aromatize_already_aromatic_unchanged);

    printf("\nKekulize Tests:\n");
    RUN_TEST(kekulize_benzene_aromatic_to_kekule);
    RUN_TEST(kekulize_pyridine);

    printf("\nNeutralize Tests:\n");
    RUN_TEST(neutralize_phenolate_to_phenol);
    RUN_TEST(neutralize_ammonium_to_amine);
    RUN_TEST(neutralize_carboxylate);
    RUN_TEST(neutralize_preserve_quaternary_nitrogen);

    printf("\nRemove Stereo Tests:\n");
    RUN_TEST(remove_stereo_ez_double_bond);
    RUN_TEST(remove_stereo_chiral_center);
    RUN_TEST(remove_stereo_complex_molecule);

    printf("\nRemove Isotopes Tests:\n");
    RUN_TEST(remove_isotopes_carbon_13);
    RUN_TEST(remove_isotopes_deuterium);
    RUN_TEST(remove_isotopes_multiple);

    printf("\nRemove H Tests:\n");
    RUN_TEST(remove_h_explicit_hydrogens);

    printf("\nComplete Sanitization Tests:\n");
    RUN_TEST(complete_sanitize_drug_salt);
    RUN_TEST(complete_sanitize_potassium_phenolate);
    RUN_TEST(complete_sanitize_multiple_salts);

    printf("\nFragment Function Tests:\n");
    RUN_TEST(fragment_info_single_fragment);
    RUN_TEST(fragment_info_multiple_fragments);
    RUN_TEST(fragment_is_salt_sodium);
    RUN_TEST(fragment_is_organic);

    printf("\nTautomer Tests:\n");
    RUN_TEST(tautomer_enumerate_basic);
    RUN_TEST(tautomer_canonical);

    printf("\nFlag Parsing Tests:\n");
    RUN_TEST(parse_flags_single);
    RUN_TEST(parse_flags_multiple);
    RUN_TEST(parse_flags_preset_complete);
    RUN_TEST(parse_flags_preset_all);
    RUN_TEST(parse_flags_invalid);
    RUN_TEST(flags_to_string);

    printf("\n=== Results ===\n");
    printf("Total: %d, Passed: %d, Failed: %d\n\n", test_count, pass_count, fail_count);

    return fail_count > 0 ? 1 : 0;
}
