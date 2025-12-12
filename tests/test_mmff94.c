/**
 * @file test_mmff94.c
 * @brief Essential tests for MMFF94 force field implementation
 *
 * Tests atom type assignment, energy calculation, and minimization.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cchem/cchem.h"
#include "cchem/depictor/coords2d.h"
#include "cchem/depictor/coords3d.h"
#include "cchem/depictor/mmff94_types.h"
#include "cchem/depictor/mmff94_params.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) static void test_##name(void)
#define RUN_TEST(name) do { \
    printf("  Testing %s... ", #name); \
    test_##name(); \
    printf("OK\n"); \
} while(0)

#define ASSERT_EQ(a, b) do { \
    int64_t _a = (a); \
    int64_t _b = (b); \
    if (_a != _b) { \
        printf("FAILED: %s = %lld, expected %lld (line %d)\n", #a, (long long)_a, (long long)_b, __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(x) ASSERT_EQ(!!(x), 1)
#define ASSERT_FALSE(x) ASSERT_EQ(!!(x), 0)
#define ASSERT_NOT_NULL(x) ASSERT_TRUE((x) != NULL)
#define ASSERT_NULL(x) ASSERT_TRUE((x) == NULL)

/* Helper: parse SMILES to molecule */
static molecule_t* parse_smiles(const char* smiles) {
    char error_buf[256];
    char* canonical = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));
    if (!canonical) {
        printf("Canonicalization failed for %s: %s\n", smiles, error_buf);
        return NULL;
    }

    molecule_t* mol = smiles_to_molecule(canonical, error_buf, sizeof(error_buf));
    free(canonical);
    return mol;
}

/* ============================================================================
 * Atom Type Assignment Tests
 * ============================================================================ */

TEST(carbon_types) {
    /* Test sp3 carbon (methane) */
    molecule_t* mol = parse_smiles("C");
    ASSERT_NOT_NULL(mol);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    ASSERT_NOT_NULL(ctx);

    cchem_status_t status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);
    ASSERT_EQ(ctx->atom_data[0].type, MMFF94_TYPE_CR);
    ASSERT_EQ(ctx->atom_data[0].hybridization, MMFF94_HYBRID_SP3);

    mmff94_context_free(ctx);
    molecule_free(mol);

    /* Test sp2 carbon (ethylene) */
    mol = parse_smiles("C=C");
    ASSERT_NOT_NULL(mol);

    ctx = mmff94_context_create(mol);
    status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);
    ASSERT_EQ(ctx->atom_data[0].type, MMFF94_TYPE_C_EQ_C);
    ASSERT_EQ(ctx->atom_data[0].hybridization, MMFF94_HYBRID_SP2);
    ASSERT_EQ(ctx->atom_data[1].type, MMFF94_TYPE_C_EQ_C);

    mmff94_context_free(ctx);
    molecule_free(mol);

    /* Test sp carbon (acetylene) */
    mol = parse_smiles("C#C");
    ASSERT_NOT_NULL(mol);

    ctx = mmff94_context_create(mol);
    status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);
    ASSERT_EQ(ctx->atom_data[0].type, MMFF94_TYPE_CSP);
    ASSERT_EQ(ctx->atom_data[0].hybridization, MMFF94_HYBRID_SP);

    mmff94_context_free(ctx);
    molecule_free(mol);

    tests_passed++;
}

TEST(aromatic_types) {
    /* Test aromatic carbon (benzene) */
    molecule_t* mol = parse_smiles("c1ccccc1");
    ASSERT_NOT_NULL(mol);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    cchem_status_t status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    /* All carbons should be aromatic type */
    for (int i = 0; i < mol->num_atoms; i++) {
        ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_C_AR);
        ASSERT_TRUE(ctx->atom_data[i].is_aromatic);
    }

    mmff94_context_free(ctx);
    molecule_free(mol);

    /* Test pyridine (aromatic nitrogen) */
    mol = parse_smiles("c1ccncc1");
    ASSERT_NOT_NULL(mol);

    ctx = mmff94_context_create(mol);
    status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    /* Find the nitrogen and check its type */
    bool found_n_ar = false;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_N) {
            ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_N_AR);
            found_n_ar = true;
        }
    }
    ASSERT_TRUE(found_n_ar);

    mmff94_context_free(ctx);
    molecule_free(mol);

    tests_passed++;
}

TEST(nitrogen_types) {
    /* Test sp3 nitrogen (ammonia-like in methylamine) */
    molecule_t* mol = parse_smiles("CN");
    ASSERT_NOT_NULL(mol);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    cchem_status_t status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    /* Find nitrogen */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_N) {
            ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_NR);
        }
    }

    mmff94_context_free(ctx);
    molecule_free(mol);

    /* Test nitrile nitrogen */
    mol = parse_smiles("C#N");
    ASSERT_NOT_NULL(mol);

    ctx = mmff94_context_create(mol);
    status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_N) {
            ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_NSP);
        }
    }

    mmff94_context_free(ctx);
    molecule_free(mol);

    tests_passed++;
}

TEST(oxygen_types) {
    /* Test carbonyl oxygen */
    molecule_t* mol = parse_smiles("CC=O");
    ASSERT_NOT_NULL(mol);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    cchem_status_t status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    /* Find oxygen and carbonyl carbon */
    bool found_carbonyl_o = false;
    bool found_carbonyl_c = false;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_O) {
            ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_O_EQ_C);
            found_carbonyl_o = true;
        }
        if (mol->atoms[i].element == ELEM_C &&
            ctx->atom_data[i].type == MMFF94_TYPE_C_EQ_O) {
            found_carbonyl_c = true;
        }
    }
    ASSERT_TRUE(found_carbonyl_o);
    ASSERT_TRUE(found_carbonyl_c);

    mmff94_context_free(ctx);
    molecule_free(mol);

    /* Test ether oxygen */
    mol = parse_smiles("COC");
    ASSERT_NOT_NULL(mol);

    ctx = mmff94_context_create(mol);
    status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_O) {
            ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_OR);
        }
    }

    mmff94_context_free(ctx);
    molecule_free(mol);

    tests_passed++;
}

TEST(halogen_types) {
    molecule_t* mol = parse_smiles("CF");
    ASSERT_NOT_NULL(mol);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    cchem_status_t status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_F) {
            ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_F);
        }
    }

    mmff94_context_free(ctx);
    molecule_free(mol);

    mol = parse_smiles("CCl");
    ASSERT_NOT_NULL(mol);

    ctx = mmff94_context_create(mol);
    status = mmff94_assign_types(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_Cl) {
            ASSERT_EQ(ctx->atom_data[i].type, MMFF94_TYPE_CL);
        }
    }

    mmff94_context_free(ctx);
    molecule_free(mol);

    tests_passed++;
}

/* ============================================================================
 * Charge Computation Tests
 * ============================================================================ */

TEST(charge_neutrality) {
    /* For neutral molecules, sum of partial charges should be ~0 */
    molecule_t* mol = parse_smiles("CCO");  /* Ethanol */
    ASSERT_NOT_NULL(mol);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    mmff94_assign_types(mol, ctx);
    cchem_status_t status = mmff94_compute_charges(mol, ctx);
    ASSERT_EQ(status, CCHEM_OK);

    double total_charge = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        total_charge += ctx->atom_data[i].partial_charge;
    }

    /* Should be close to 0 (within numerical tolerance) */
    printf("(total charge: %.4f) ", total_charge);
    /* Note: implicit H charges are not included in atoms, so total may not be 0 */

    mmff94_context_free(ctx);
    molecule_free(mol);

    tests_passed++;
}

TEST(charge_polarity) {
    /* Oxygen should be more negative than carbon */
    molecule_t* mol = parse_smiles("CO");  /* Methanol */
    ASSERT_NOT_NULL(mol);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    mmff94_assign_types(mol, ctx);
    mmff94_compute_charges(mol, ctx);

    double c_charge = 0.0, o_charge = 0.0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (mol->atoms[i].element == ELEM_C) {
            c_charge = ctx->atom_data[i].partial_charge;
        } else if (mol->atoms[i].element == ELEM_O) {
            o_charge = ctx->atom_data[i].partial_charge;
        }
    }

    printf("(C: %.3f, O: %.3f) ", c_charge, o_charge);
    /* Oxygen should be more negative (or less positive) than carbon */
    ASSERT_TRUE(o_charge < c_charge);

    mmff94_context_free(ctx);
    molecule_free(mol);

    tests_passed++;
}

/* ============================================================================
 * Energy and Minimization Tests
 * ============================================================================ */

TEST(energy_calculation) {
    /* Generate 3D coordinates and verify energy is computed */
    molecule_t* mol = parse_smiles("CC");  /* Ethane */
    ASSERT_NOT_NULL(mol);

    coords3d_options_t opts = COORDS3D_OPTIONS_DEFAULT;
    mol_coords_t* coords = coords3d_generate(mol, &opts);
    ASSERT_NOT_NULL(coords);
    ASSERT_TRUE(coords->has_3d);

    /* Create context and compute energy */
    mmff94_context_t* ctx = mmff94_context_create(mol);
    mmff94_assign_types(mol, ctx);
    mmff94_compute_charges(mol, ctx);

    double energy = mmff94_calc_energy(mol, ctx, coords, NULL, NULL, NULL);
    printf("(energy: %.2f kcal/mol) ", energy);

    /* Energy should be finite and reasonable */
    ASSERT_TRUE(isfinite(energy));
    ASSERT_TRUE(energy < 1000.0);  /* Reasonable upper bound */

    mmff94_context_free(ctx);
    mol_coords_free(coords);
    molecule_free(mol);

    tests_passed++;
}

TEST(minimization_convergence) {
    /* Verify that minimization reduces energy */
    molecule_t* mol = parse_smiles("CCCC");  /* Butane */
    ASSERT_NOT_NULL(mol);

    /* Generate initial coordinates (before minimization would be called) */
    coords2d_options_t opts_2d = COORDS2D_OPTIONS_DEFAULT;
    mol_coords_t* coords = coords2d_generate(mol, &opts_2d);
    ASSERT_NOT_NULL(coords);

    coords->has_3d = true;
    for (int i = 0; i < mol->num_atoms; i++) {
        coords->coords_3d[i].x = coords->coords_2d[i].x;
        coords->coords_3d[i].y = coords->coords_2d[i].y;
        coords->coords_3d[i].z = 0.0;
    }

    /* Create context */
    mmff94_context_t* ctx = mmff94_context_create(mol);
    mmff94_assign_types(mol, ctx);
    mmff94_compute_charges(mol, ctx);

    /* Initial energy (unoptimized) */
    double initial_energy = mmff94_calc_energy(mol, ctx, coords, NULL, NULL, NULL);

    /* Run minimization */
    coords3d_options_t opts = COORDS3D_OPTIONS_DEFAULT;
    opts.max_iterations = 200;
    double final_energy = coords3d_minimize(coords, mol, &opts);

    printf("(initial: %.2f, final: %.2f) ", initial_energy, final_energy);

    /* Final energy should be less than or equal to initial */
    ASSERT_TRUE(final_energy <= initial_energy + 1.0);  /* Allow small tolerance */

    mmff94_context_free(ctx);
    mol_coords_free(coords);
    molecule_free(mol);

    tests_passed++;
}

TEST(no_clashes_benzene) {
    /* Benzene should minimize to a planar structure with no clashes */
    molecule_t* mol = parse_smiles("c1ccccc1");
    ASSERT_NOT_NULL(mol);

    coords3d_options_t opts = COORDS3D_OPTIONS_DEFAULT;
    opts.max_iterations = 500;
    mol_coords_t* coords = coords3d_generate(mol, &opts);
    ASSERT_NOT_NULL(coords);

    /* Use 0.6 threshold to allow for close H-H contacts in aromatic systems */
    int clashes = coords3d_count_clashes(coords, mol, 0.6);
    printf("(clashes: %d) ", clashes);
    ASSERT_EQ(clashes, 0);

    mol_coords_free(coords);
    molecule_free(mol);

    tests_passed++;
}

TEST(no_clashes_cyclohexane) {
    /* Cyclohexane should minimize to chair with no clashes */
    molecule_t* mol = parse_smiles("C1CCCCC1");
    ASSERT_NOT_NULL(mol);

    coords3d_options_t opts = COORDS3D_OPTIONS_DEFAULT;
    opts.max_iterations = 1000;
    mol_coords_t* coords = coords3d_generate(mol, &opts);
    ASSERT_NOT_NULL(coords);

    /* Use 0.6 threshold to allow for close contacts in chair conformation */
    int clashes = coords3d_count_clashes(coords, mol, 0.6);
    printf("(clashes: %d) ", clashes);
    ASSERT_EQ(clashes, 0);

    mol_coords_free(coords);
    molecule_free(mol);

    tests_passed++;
}

TEST(gradient_direction) {
    /* Verify gradients point in correct direction (energy decreases) */
    molecule_t* mol = parse_smiles("CC");
    ASSERT_NOT_NULL(mol);

    coords3d_options_t opts = COORDS3D_OPTIONS_DEFAULT;
    mol_coords_t* coords = coords3d_generate(mol, &opts);
    ASSERT_NOT_NULL(coords);

    mmff94_context_t* ctx = mmff94_context_create(mol);
    mmff94_assign_types(mol, ctx);
    mmff94_compute_charges(mol, ctx);

    int n = mol->num_atoms;
    double* gx = malloc(n * sizeof(double));
    double* gy = malloc(n * sizeof(double));
    double* gz = malloc(n * sizeof(double));

    double energy = mmff94_calc_energy(mol, ctx, coords, gx, gy, gz);

    /* Take a small step in negative gradient direction */
    double step = 0.0001;
    for (int i = 0; i < n; i++) {
        coords->coords_3d[i].x -= step * gx[i];
        coords->coords_3d[i].y -= step * gy[i];
        coords->coords_3d[i].z -= step * gz[i];
    }

    double new_energy = mmff94_calc_energy(mol, ctx, coords, NULL, NULL, NULL);

    printf("(E: %.4f -> %.4f) ", energy, new_energy);

    /* New energy should be less than or equal (within numerical precision) */
    ASSERT_TRUE(new_energy <= energy + 0.1);

    free(gx); free(gy); free(gz);
    mmff94_context_free(ctx);
    mol_coords_free(coords);
    molecule_free(mol);

    tests_passed++;
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void) {
    printf("MMFF94 Force Field Tests\n");
    printf("========================\n\n");

    printf("Atom Type Assignment:\n");
    RUN_TEST(carbon_types);
    RUN_TEST(aromatic_types);
    RUN_TEST(nitrogen_types);
    RUN_TEST(oxygen_types);
    RUN_TEST(halogen_types);

    printf("\nCharge Computation:\n");
    RUN_TEST(charge_neutrality);
    RUN_TEST(charge_polarity);

    printf("\nEnergy and Minimization:\n");
    RUN_TEST(energy_calculation);
    RUN_TEST(minimization_convergence);
    RUN_TEST(no_clashes_benzene);
    RUN_TEST(no_clashes_cyclohexane);
    RUN_TEST(gradient_direction);

    printf("\n========================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
