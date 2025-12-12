/**
 * @file test_descriptors.c
 * @brief Comprehensive tests for molecular descriptors
 *
 * Tests all 52 count-based descriptors with known molecules
 * to verify correct computation from parsed SMILES.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cchem/cchem.h"
#include "cchem/descriptors.h"

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
#define ASSERT_NOT_NULL(x) ASSERT_TRUE((x) != NULL)

/* Helper: compute descriptor from SMILES */
static int64_t get_desc_value(const char* smiles, const char* desc_name) {
    char error_buf[256];
    descriptor_value_t value = {0};

    cchem_status_t status = descriptor_compute_from_smiles(smiles, desc_name, &value,
                                                           error_buf, sizeof(error_buf));
    if (status != CCHEM_OK) {
        printf("Error computing %s for %s: %s\n", desc_name, smiles, error_buf);
        return -999;
    }
    return value.i;
}

/* Helper: parse molecule for direct testing */
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

/* ===========================================================================
 * Descriptor Registration Tests
 * =========================================================================== */

TEST(descriptor_init) {
    descriptors_init();
    ASSERT_TRUE(descriptor_count() > 0);
    printf("(%d descriptors registered) ", descriptor_count());
    tests_passed++;
}

TEST(descriptor_lookup) {
    const descriptor_def_t* def = descriptor_get("CarbonCount");
    ASSERT_NOT_NULL(def);
    ASSERT_TRUE(strcmp(def->name, "CarbonCount") == 0);

    /* Case insensitive */
    def = descriptor_get("carboncount");
    ASSERT_NOT_NULL(def);

    def = descriptor_get("CARBONCOUNT");
    ASSERT_NOT_NULL(def);

    /* Non-existent */
    def = descriptor_get("NonExistent");
    ASSERT_TRUE(def == NULL);

    tests_passed++;
}

TEST(descriptor_list) {
    const descriptor_def_t* descs[MAX_DESCRIPTORS];
    int count = descriptor_get_all(descs, MAX_DESCRIPTORS);
    ASSERT_TRUE(count >= 50);  /* Should have at least 50 descriptors */

    /* Verify all have valid compute functions */
    for (int i = 0; i < count; i++) {
        ASSERT_NOT_NULL(descs[i]->compute);
        ASSERT_TRUE(strlen(descs[i]->name) > 0);
    }

    tests_passed++;
}

/* ===========================================================================
 * Element Count Tests
 * =========================================================================== */

TEST(carbon_count) {
    /* Methane: 1 C */
    ASSERT_EQ(get_desc_value("C", "CarbonCount"), 1);
    /* Ethane: 2 C */
    ASSERT_EQ(get_desc_value("CC", "CarbonCount"), 2);
    /* Benzene: 6 C */
    ASSERT_EQ(get_desc_value("c1ccccc1", "CarbonCount"), 6);
    /* No carbons */
    ASSERT_EQ(get_desc_value("[Na+]", "CarbonCount"), 0);

    tests_passed++;
}

TEST(hydrogen_count) {
    /* Methane CH4: 4 H */
    ASSERT_EQ(get_desc_value("C", "HydrogenCount"), 4);
    /* Ethane C2H6: 6 H */
    ASSERT_EQ(get_desc_value("CC", "HydrogenCount"), 6);
    /* Benzene C6H6: 6 H */
    ASSERT_EQ(get_desc_value("c1ccccc1", "HydrogenCount"), 6);
    /* Ethene C2H4: 4 H */
    ASSERT_EQ(get_desc_value("C=C", "HydrogenCount"), 4);
    /* Ethyne C2H2: 2 H */
    ASSERT_EQ(get_desc_value("C#C", "HydrogenCount"), 2);

    tests_passed++;
}

TEST(oxygen_count) {
    ASSERT_EQ(get_desc_value("C", "OxygenCount"), 0);
    ASSERT_EQ(get_desc_value("CO", "OxygenCount"), 1);
    ASSERT_EQ(get_desc_value("CC(=O)O", "OxygenCount"), 2);
    ASSERT_EQ(get_desc_value("O", "OxygenCount"), 1);

    tests_passed++;
}

TEST(nitrogen_count) {
    ASSERT_EQ(get_desc_value("C", "NitrogenCount"), 0);
    ASSERT_EQ(get_desc_value("CN", "NitrogenCount"), 1);
    ASSERT_EQ(get_desc_value("c1ccncc1", "NitrogenCount"), 1);  /* Pyridine */
    ASSERT_EQ(get_desc_value("c1cncnc1", "NitrogenCount"), 2);  /* Pyrimidine */

    tests_passed++;
}

TEST(halogen_counts) {
    ASSERT_EQ(get_desc_value("CCl", "ChlorineCount"), 1);
    ASSERT_EQ(get_desc_value("ClC(Cl)(Cl)Cl", "ChlorineCount"), 4);

    ASSERT_EQ(get_desc_value("CF", "FluorineCount"), 1);
    ASSERT_EQ(get_desc_value("CBr", "BromineCount"), 1);
    ASSERT_EQ(get_desc_value("CI", "IodineCount"), 1);

    /* Multiple halogens */
    ASSERT_EQ(get_desc_value("ClCBr", "ChlorineCount"), 1);
    ASSERT_EQ(get_desc_value("ClCBr", "BromineCount"), 1);

    tests_passed++;
}

TEST(sulfur_phosphorus_count) {
    ASSERT_EQ(get_desc_value("CSC", "SulfurCount"), 1);  /* DMS */
    ASSERT_EQ(get_desc_value("c1ccsc1", "SulfurCount"), 1);  /* Thiophene */
    ASSERT_EQ(get_desc_value("OP(=O)(O)O", "PhosphorusCount"), 1);  /* Phosphoric acid */

    tests_passed++;
}

/* ===========================================================================
 * Heteroatom Tests
 * =========================================================================== */

TEST(heteroatom_count) {
    /* Methane: no heteroatoms */
    ASSERT_EQ(get_desc_value("C", "HeteroatomCount"), 0);
    /* Methanol: 1 O */
    ASSERT_EQ(get_desc_value("CO", "HeteroatomCount"), 1);
    /* Pyridine: 1 N */
    ASSERT_EQ(get_desc_value("c1ccncc1", "HeteroatomCount"), 1);
    /* Acetic acid: 2 O */
    ASSERT_EQ(get_desc_value("CC(=O)O", "HeteroatomCount"), 2);
    /* Chlorobenzene: 1 Cl */
    ASSERT_EQ(get_desc_value("Clc1ccccc1", "HeteroatomCount"), 1);

    tests_passed++;
}

TEST(halogen_total_count) {
    ASSERT_EQ(get_desc_value("C", "HalogenCount"), 0);
    ASSERT_EQ(get_desc_value("CCl", "HalogenCount"), 1);
    ASSERT_EQ(get_desc_value("ClC(Cl)Cl", "HalogenCount"), 3);  /* Chloroform */
    ASSERT_EQ(get_desc_value("ClCBr", "HalogenCount"), 2);
    ASSERT_EQ(get_desc_value("FC(F)(F)Cl", "HalogenCount"), 4);

    tests_passed++;
}

/* ===========================================================================
 * Atom Count Tests
 * =========================================================================== */

TEST(atom_count) {
    /* AtomCount = heavy atoms only */
    ASSERT_EQ(get_desc_value("C", "AtomCount"), 1);
    ASSERT_EQ(get_desc_value("CC", "AtomCount"), 2);
    ASSERT_EQ(get_desc_value("c1ccccc1", "AtomCount"), 6);
    ASSERT_EQ(get_desc_value("CCO", "AtomCount"), 3);

    /* HeavyAtomCount should equal AtomCount */
    ASSERT_EQ(get_desc_value("CCO", "HeavyAtomCount"), 3);

    tests_passed++;
}

TEST(aromatic_atom_count) {
    ASSERT_EQ(get_desc_value("C", "AromaticAtomCount"), 0);
    ASSERT_EQ(get_desc_value("c1ccccc1", "AromaticAtomCount"), 6);  /* Benzene */
    ASSERT_EQ(get_desc_value("c1ccncc1", "AromaticAtomCount"), 6);  /* Pyridine */
    ASSERT_EQ(get_desc_value("c1ccc2ccccc2c1", "AromaticAtomCount"), 10);  /* Naphthalene */
    ASSERT_EQ(get_desc_value("Cc1ccccc1", "AromaticAtomCount"), 6);  /* Toluene (methyl not aromatic) */

    tests_passed++;
}

/* ===========================================================================
 * Carbon Hybridization Tests
 * =========================================================================== */

TEST(sp3_carbon_count) {
    /* Methane: sp3 */
    ASSERT_EQ(get_desc_value("C", "SP3CarbonCount"), 1);
    /* Ethane: 2 sp3 */
    ASSERT_EQ(get_desc_value("CC", "SP3CarbonCount"), 2);
    /* Ethene: 0 sp3 */
    ASSERT_EQ(get_desc_value("C=C", "SP3CarbonCount"), 0);
    /* Ethanol: 2 sp3 */
    ASSERT_EQ(get_desc_value("CCO", "SP3CarbonCount"), 2);
    /* Toluene: 1 sp3 (methyl) */
    ASSERT_EQ(get_desc_value("Cc1ccccc1", "SP3CarbonCount"), 1);

    tests_passed++;
}

TEST(sp2_carbon_count) {
    /* Methane: 0 sp2 */
    ASSERT_EQ(get_desc_value("C", "SP2CarbonCount"), 0);
    /* Ethene: 2 sp2 */
    ASSERT_EQ(get_desc_value("C=C", "SP2CarbonCount"), 2);
    /* Benzene: 6 sp2 */
    ASSERT_EQ(get_desc_value("c1ccccc1", "SP2CarbonCount"), 6);
    /* Acetone: 1 sp2 (carbonyl C) */
    ASSERT_EQ(get_desc_value("CC(=O)C", "SP2CarbonCount"), 1);

    tests_passed++;
}

TEST(sp_carbon_count) {
    /* Ethyne: 2 sp */
    ASSERT_EQ(get_desc_value("C#C", "SPCarbonCount"), 2);
    /* Propyne: 2 sp */
    ASSERT_EQ(get_desc_value("CC#C", "SPCarbonCount"), 2);
    /* Acetonitrile: 1 sp */
    ASSERT_EQ(get_desc_value("CC#N", "SPCarbonCount"), 1);
    /* No sp carbons */
    ASSERT_EQ(get_desc_value("CC", "SPCarbonCount"), 0);

    tests_passed++;
}

/* ===========================================================================
 * Bond Count Tests
 * =========================================================================== */

TEST(total_bond_count) {
    ASSERT_EQ(get_desc_value("C", "TotalBondCount"), 0);  /* No explicit bonds */
    ASSERT_EQ(get_desc_value("CC", "TotalBondCount"), 1);
    ASSERT_EQ(get_desc_value("CCC", "TotalBondCount"), 2);
    ASSERT_EQ(get_desc_value("c1ccccc1", "TotalBondCount"), 6);

    tests_passed++;
}

TEST(single_bond_count) {
    ASSERT_EQ(get_desc_value("CC", "SingleBondCount"), 1);
    ASSERT_EQ(get_desc_value("CCC", "SingleBondCount"), 2);
    ASSERT_EQ(get_desc_value("C=C", "SingleBondCount"), 0);

    tests_passed++;
}

TEST(double_bond_count) {
    ASSERT_EQ(get_desc_value("C=C", "DoubleBondCount"), 1);
    ASSERT_EQ(get_desc_value("C=CC=C", "DoubleBondCount"), 2);
    ASSERT_EQ(get_desc_value("CC(=O)O", "DoubleBondCount"), 1);  /* C=O in acid */
    ASSERT_EQ(get_desc_value("CC", "DoubleBondCount"), 0);

    tests_passed++;
}

TEST(triple_bond_count) {
    ASSERT_EQ(get_desc_value("C#C", "TripleBondCount"), 1);
    ASSERT_EQ(get_desc_value("CC#N", "TripleBondCount"), 1);
    ASSERT_EQ(get_desc_value("C#CC#C", "TripleBondCount"), 2);
    ASSERT_EQ(get_desc_value("CC", "TripleBondCount"), 0);

    tests_passed++;
}

TEST(aromatic_bond_count) {
    ASSERT_EQ(get_desc_value("c1ccccc1", "AromaticBondCount"), 6);
    ASSERT_EQ(get_desc_value("c1ccncc1", "AromaticBondCount"), 6);  /* Pyridine */
    ASSERT_EQ(get_desc_value("c1ccc2ccccc2c1", "AromaticBondCount"), 11);  /* Naphthalene */
    ASSERT_EQ(get_desc_value("CC", "AromaticBondCount"), 0);

    tests_passed++;
}

TEST(ring_bond_count) {
    ASSERT_EQ(get_desc_value("CC", "RingBondCount"), 0);
    ASSERT_EQ(get_desc_value("C1CC1", "RingBondCount"), 3);
    ASSERT_EQ(get_desc_value("c1ccccc1", "RingBondCount"), 6);
    ASSERT_EQ(get_desc_value("C1CCCCC1", "RingBondCount"), 6);

    tests_passed++;
}

TEST(rotatable_bond_count) {
    /* Simple cases */
    ASSERT_EQ(get_desc_value("C", "RotatableBondCount"), 0);
    ASSERT_EQ(get_desc_value("CC", "RotatableBondCount"), 0);  /* Terminal methyls */
    /* Chain with rotatable bonds */
    ASSERT_EQ(get_desc_value("CCCC", "RotatableBondCount"), 1);  /* Central C-C */
    ASSERT_EQ(get_desc_value("CCCCC", "RotatableBondCount"), 2);
    /* Ring bonds not rotatable */
    ASSERT_EQ(get_desc_value("c1ccccc1", "RotatableBondCount"), 0);

    tests_passed++;
}

/* ===========================================================================
 * Ring Count Tests
 * =========================================================================== */

TEST(ring_count) {
    ASSERT_EQ(get_desc_value("CC", "RingCount"), 0);
    ASSERT_EQ(get_desc_value("C1CC1", "RingCount"), 1);
    ASSERT_EQ(get_desc_value("c1ccccc1", "RingCount"), 1);
    ASSERT_EQ(get_desc_value("c1ccc2ccccc2c1", "RingCount"), 2);  /* Naphthalene */

    tests_passed++;
}

TEST(aromatic_ring_count) {
    ASSERT_EQ(get_desc_value("C1CC1", "AromaticRingCount"), 0);
    ASSERT_EQ(get_desc_value("c1ccccc1", "AromaticRingCount"), 1);
    ASSERT_EQ(get_desc_value("c1ccc2ccccc2c1", "AromaticRingCount"), 2);  /* Naphthalene */
    ASSERT_EQ(get_desc_value("C1CCCCC1", "AromaticRingCount"), 0);  /* Cyclohexane */

    tests_passed++;
}

TEST(aliphatic_ring_count) {
    ASSERT_EQ(get_desc_value("C1CC1", "AliphaticRingCount"), 1);
    ASSERT_EQ(get_desc_value("C1CCCCC1", "AliphaticRingCount"), 1);
    ASSERT_EQ(get_desc_value("c1ccccc1", "AliphaticRingCount"), 0);

    tests_passed++;
}

TEST(heterocycle_count) {
    ASSERT_EQ(get_desc_value("c1ccccc1", "HeterocycleCount"), 0);  /* Benzene: no heteroatoms */
    ASSERT_EQ(get_desc_value("c1ccncc1", "HeterocycleCount"), 1);  /* Pyridine */
    ASSERT_EQ(get_desc_value("c1ccoc1", "HeterocycleCount"), 1);   /* Furan */
    ASSERT_EQ(get_desc_value("c1ccsc1", "HeterocycleCount"), 1);   /* Thiophene */

    tests_passed++;
}

TEST(aromatic_heterocycle_count) {
    ASSERT_EQ(get_desc_value("c1ccncc1", "AromaticHeterocycleCount"), 1);  /* Pyridine */
    ASSERT_EQ(get_desc_value("c1ccoc1", "AromaticHeterocycleCount"), 1);   /* Furan */
    ASSERT_EQ(get_desc_value("c1ccccc1", "AromaticHeterocycleCount"), 0);  /* Benzene */

    tests_passed++;
}

/* ===========================================================================
 * Ring Size Tests
 * =========================================================================== */

TEST(ring_size_counts) {
    /* 3-membered ring */
    ASSERT_EQ(get_desc_value("C1CC1", "Ring3Count"), 1);
    ASSERT_EQ(get_desc_value("C1CC1", "Ring4Count"), 0);

    /* 4-membered ring */
    ASSERT_EQ(get_desc_value("C1CCC1", "Ring4Count"), 1);
    ASSERT_EQ(get_desc_value("C1CCC1", "Ring3Count"), 0);

    /* 5-membered ring */
    ASSERT_EQ(get_desc_value("C1CCCC1", "Ring5Count"), 1);
    ASSERT_EQ(get_desc_value("c1cccc1", "Ring5Count"), 1);  /* Cyclopentadienyl */

    /* 6-membered ring */
    ASSERT_EQ(get_desc_value("c1ccccc1", "Ring6Count"), 1);
    ASSERT_EQ(get_desc_value("C1CCCCC1", "Ring6Count"), 1);

    tests_passed++;
}

TEST(large_ring_count) {
    /* Normal rings */
    ASSERT_EQ(get_desc_value("c1ccccc1", "LargeRingCount"), 0);
    ASSERT_EQ(get_desc_value("C1CCCCC1", "LargeRingCount"), 0);

    /* 7-membered ring */
    ASSERT_EQ(get_desc_value("C1CCCCCC1", "LargeRingCount"), 1);

    tests_passed++;
}

/* ===========================================================================
 * Connectivity Tests
 * =========================================================================== */

TEST(ring_chain_atom_count) {
    /* All chain */
    ASSERT_EQ(get_desc_value("CCC", "ChainAtomCount"), 3);
    ASSERT_EQ(get_desc_value("CCC", "RingAtomCount"), 0);

    /* All ring */
    ASSERT_EQ(get_desc_value("C1CC1", "RingAtomCount"), 3);
    ASSERT_EQ(get_desc_value("C1CC1", "ChainAtomCount"), 0);

    /* Mixed */
    ASSERT_EQ(get_desc_value("Cc1ccccc1", "RingAtomCount"), 6);
    ASSERT_EQ(get_desc_value("Cc1ccccc1", "ChainAtomCount"), 1);

    tests_passed++;
}

TEST(terminal_atom_count) {
    /* Single atom: not terminal (no bond) */
    ASSERT_EQ(get_desc_value("C", "TerminalAtomCount"), 0);
    /* Ethane: 2 terminals */
    ASSERT_EQ(get_desc_value("CC", "TerminalAtomCount"), 2);
    /* Propane: 2 terminals */
    ASSERT_EQ(get_desc_value("CCC", "TerminalAtomCount"), 2);
    /* Branched: more terminals */
    ASSERT_EQ(get_desc_value("CC(C)C", "TerminalAtomCount"), 3);

    tests_passed++;
}

TEST(branch_point_count) {
    /* Linear: no branches */
    ASSERT_EQ(get_desc_value("CCCC", "BranchPointCount"), 0);
    /* Isobutane: 1 branch point */
    ASSERT_EQ(get_desc_value("CC(C)C", "BranchPointCount"), 1);
    /* Neopentane: 1 branch */
    ASSERT_EQ(get_desc_value("CC(C)(C)C", "BranchPointCount"), 1);
    /* Benzene: all atoms have degree 2 */
    ASSERT_EQ(get_desc_value("c1ccccc1", "BranchPointCount"), 0);

    tests_passed++;
}

TEST(quaternary_carbon_count) {
    /* No quaternary */
    ASSERT_EQ(get_desc_value("CCCC", "QuaternaryCarbonCount"), 0);
    ASSERT_EQ(get_desc_value("CC(C)C", "QuaternaryCarbonCount"), 0);  /* Only 3 heavy neighbors */
    /* Neopentane: 1 quaternary */
    ASSERT_EQ(get_desc_value("CC(C)(C)C", "QuaternaryCarbonCount"), 1);

    tests_passed++;
}

TEST(bridgehead_spiro_count) {
    /* Simple ring: no bridgehead */
    ASSERT_EQ(get_desc_value("c1ccccc1", "BridgeheadAtomCount"), 0);

    /* Naphthalene: 2 bridgeheads (shared atoms) */
    int64_t bridgeheads = get_desc_value("c1ccc2ccccc2c1", "BridgeheadAtomCount");
    ASSERT_EQ(bridgeheads, 2);

    tests_passed++;
}

/* ===========================================================================
 * H-Bond Tests
 * =========================================================================== */

TEST(hbond_donor_count) {
    ASSERT_EQ(get_desc_value("C", "HBondDonorCount"), 0);
    ASSERT_EQ(get_desc_value("CO", "HBondDonorCount"), 1);   /* Methanol OH */
    ASSERT_EQ(get_desc_value("CN", "HBondDonorCount"), 1);   /* Methylamine NH2 */
    ASSERT_EQ(get_desc_value("OCCO", "HBondDonorCount"), 2); /* Ethylene glycol: 2 OH */

    tests_passed++;
}

TEST(hbond_acceptor_count) {
    ASSERT_EQ(get_desc_value("C", "HBondAcceptorCount"), 0);
    ASSERT_EQ(get_desc_value("CO", "HBondAcceptorCount"), 1);  /* O is acceptor */
    ASSERT_EQ(get_desc_value("CC(=O)O", "HBondAcceptorCount"), 2);  /* 2 O's */
    ASSERT_EQ(get_desc_value("c1ccncc1", "HBondAcceptorCount"), 1);  /* Pyridine N */

    tests_passed++;
}

/* ===========================================================================
 * Charge Tests
 * =========================================================================== */

TEST(charge_counts) {
    /* Neutral molecules */
    ASSERT_EQ(get_desc_value("CC", "FormalChargeCount"), 0);
    ASSERT_EQ(get_desc_value("CC", "NetCharge"), 0);

    /* Sodium ion */
    ASSERT_EQ(get_desc_value("[Na+]", "FormalChargeCount"), 1);
    ASSERT_EQ(get_desc_value("[Na+]", "PositiveChargeCount"), 1);
    ASSERT_EQ(get_desc_value("[Na+]", "NegativeChargeCount"), 0);
    ASSERT_EQ(get_desc_value("[Na+]", "NetCharge"), 1);

    /* Chloride ion */
    ASSERT_EQ(get_desc_value("[Cl-]", "FormalChargeCount"), 1);
    ASSERT_EQ(get_desc_value("[Cl-]", "PositiveChargeCount"), 0);
    ASSERT_EQ(get_desc_value("[Cl-]", "NegativeChargeCount"), 1);
    ASSERT_EQ(get_desc_value("[Cl-]", "NetCharge"), -1);

    /* Salt: net charge 0 */
    ASSERT_EQ(get_desc_value("[Na+].[Cl-]", "NetCharge"), 0);

    /* Zwitterion */
    ASSERT_EQ(get_desc_value("[NH3+]CC([O-])=O", "PositiveChargeCount"), 1);
    ASSERT_EQ(get_desc_value("[NH3+]CC([O-])=O", "NegativeChargeCount"), 1);
    ASSERT_EQ(get_desc_value("[NH3+]CC([O-])=O", "NetCharge"), 0);

    tests_passed++;
}

/* ===========================================================================
 * Stereo Tests
 * =========================================================================== */

TEST(chiral_center_count) {
    ASSERT_EQ(get_desc_value("CC", "ChiralCenterCount"), 0);
    ASSERT_EQ(get_desc_value("C[C@H](O)F", "ChiralCenterCount"), 1);
    ASSERT_EQ(get_desc_value("C[C@@H](O)F", "ChiralCenterCount"), 1);

    tests_passed++;
}

TEST(stereo_double_bond_count) {
    /* Non-stereo double bond */
    ASSERT_EQ(get_desc_value("C=C", "StereoDoubleBondCount"), 0);

    tests_passed++;
}

/* ===========================================================================
 * Fragment Tests
 * =========================================================================== */

TEST(fragment_count) {
    ASSERT_EQ(get_desc_value("CC", "FragmentCount"), 1);
    ASSERT_EQ(get_desc_value("C.C", "FragmentCount"), 2);
    ASSERT_EQ(get_desc_value("C.C.C", "FragmentCount"), 3);
    ASSERT_EQ(get_desc_value("[Na+].[Cl-]", "FragmentCount"), 2);

    tests_passed++;
}

/* ===========================================================================
 * Unsaturation Tests
 * =========================================================================== */

TEST(unsaturation_count) {
    /* Saturated: DBE = 0 */
    ASSERT_EQ(get_desc_value("C", "UnsaturationCount"), 0);     /* CH4 */
    ASSERT_EQ(get_desc_value("CC", "UnsaturationCount"), 0);    /* C2H6 */

    /* One double bond: DBE = 1 */
    ASSERT_EQ(get_desc_value("C=C", "UnsaturationCount"), 1);   /* C2H4 */

    /* One triple bond: DBE = 2 */
    ASSERT_EQ(get_desc_value("C#C", "UnsaturationCount"), 2);   /* C2H2 */

    /* Benzene: DBE = 4 (3 double bonds + 1 ring) */
    ASSERT_EQ(get_desc_value("c1ccccc1", "UnsaturationCount"), 4);

    tests_passed++;
}

/* ===========================================================================
 * Complex Molecule Tests
 * =========================================================================== */

TEST(complex_molecules) {
    /* Caffeine: Cn1cnc2c1c(=O)n(c(=O)n2C)C */
    ASSERT_EQ(get_desc_value("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "CarbonCount"), 8);
    ASSERT_EQ(get_desc_value("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "NitrogenCount"), 4);
    ASSERT_EQ(get_desc_value("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "OxygenCount"), 2);
    ASSERT_EQ(get_desc_value("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "RingCount"), 2);

    /* Aspirin: CC(=O)Oc1ccccc1C(=O)O */
    ASSERT_EQ(get_desc_value("CC(=O)Oc1ccccc1C(=O)O", "CarbonCount"), 9);
    ASSERT_EQ(get_desc_value("CC(=O)Oc1ccccc1C(=O)O", "OxygenCount"), 4);
    ASSERT_EQ(get_desc_value("CC(=O)Oc1ccccc1C(=O)O", "AromaticRingCount"), 1);

    tests_passed++;
}

/* ===========================================================================
 * Batch Computation Test
 * =========================================================================== */

TEST(batch_computation) {
    molecule_t* mol = parse_smiles("c1ccccc1");
    ASSERT_NOT_NULL(mol);

    /* Compute multiple descriptors */
    const char* names[] = {"CarbonCount", "HydrogenCount", "AromaticRingCount"};
    descriptor_value_t values[3];

    cchem_status_t status = descriptors_compute_batch(mol, names, 3, values);
    ASSERT_EQ(status, CCHEM_OK);
    ASSERT_EQ(values[0].i, 6);  /* CarbonCount */
    ASSERT_EQ(values[1].i, 6);  /* HydrogenCount */
    ASSERT_EQ(values[2].i, 1);  /* AromaticRingCount */

    molecule_free(mol);
    tests_passed++;
}

TEST(compute_all_descriptors) {
    molecule_t* mol = parse_smiles("CCO");
    ASSERT_NOT_NULL(mol);

    const char* names[MAX_DESCRIPTORS];
    descriptor_value_t values[MAX_DESCRIPTORS];

    int count = descriptors_compute_all(mol, names, values, MAX_DESCRIPTORS);
    ASSERT_TRUE(count >= 50);

    /* Verify we got valid results for all */
    for (int i = 0; i < count; i++) {
        ASSERT_NOT_NULL(names[i]);
    }

    molecule_free(mol);
    tests_passed++;
}

/* ===========================================================================
 * Main
 * =========================================================================== */

int main(void) {
    printf("Running descriptor tests...\n\n");

    /* Initialize descriptors first */
    descriptors_init();

    /* Registration tests */
    printf("Registration:\n");
    RUN_TEST(descriptor_init);
    RUN_TEST(descriptor_lookup);
    RUN_TEST(descriptor_list);

    /* Element count tests */
    printf("\nElement counts:\n");
    RUN_TEST(carbon_count);
    RUN_TEST(hydrogen_count);
    RUN_TEST(oxygen_count);
    RUN_TEST(nitrogen_count);
    RUN_TEST(halogen_counts);
    RUN_TEST(sulfur_phosphorus_count);

    /* Heteroatom tests */
    printf("\nHeteroatom counts:\n");
    RUN_TEST(heteroatom_count);
    RUN_TEST(halogen_total_count);

    /* Atom count tests */
    printf("\nAtom counts:\n");
    RUN_TEST(atom_count);
    RUN_TEST(aromatic_atom_count);

    /* Carbon hybridization tests */
    printf("\nCarbon hybridization:\n");
    RUN_TEST(sp3_carbon_count);
    RUN_TEST(sp2_carbon_count);
    RUN_TEST(sp_carbon_count);

    /* Bond count tests */
    printf("\nBond counts:\n");
    RUN_TEST(total_bond_count);
    RUN_TEST(single_bond_count);
    RUN_TEST(double_bond_count);
    RUN_TEST(triple_bond_count);
    RUN_TEST(aromatic_bond_count);
    RUN_TEST(ring_bond_count);
    RUN_TEST(rotatable_bond_count);

    /* Ring count tests */
    printf("\nRing counts:\n");
    RUN_TEST(ring_count);
    RUN_TEST(aromatic_ring_count);
    RUN_TEST(aliphatic_ring_count);
    RUN_TEST(heterocycle_count);
    RUN_TEST(aromatic_heterocycle_count);
    RUN_TEST(ring_size_counts);
    RUN_TEST(large_ring_count);

    /* Connectivity tests */
    printf("\nConnectivity:\n");
    RUN_TEST(ring_chain_atom_count);
    RUN_TEST(terminal_atom_count);
    RUN_TEST(branch_point_count);
    RUN_TEST(quaternary_carbon_count);
    RUN_TEST(bridgehead_spiro_count);

    /* H-bond tests */
    printf("\nH-bonding:\n");
    RUN_TEST(hbond_donor_count);
    RUN_TEST(hbond_acceptor_count);

    /* Charge tests */
    printf("\nCharge:\n");
    RUN_TEST(charge_counts);

    /* Stereo tests */
    printf("\nStereochemistry:\n");
    RUN_TEST(chiral_center_count);
    RUN_TEST(stereo_double_bond_count);

    /* Fragment tests */
    printf("\nFragments:\n");
    RUN_TEST(fragment_count);

    /* Unsaturation tests */
    printf("\nUnsaturation:\n");
    RUN_TEST(unsaturation_count);

    /* Complex molecule tests */
    printf("\nComplex molecules:\n");
    RUN_TEST(complex_molecules);

    /* Batch computation tests */
    printf("\nBatch computation:\n");
    RUN_TEST(batch_computation);
    RUN_TEST(compute_all_descriptors);

    /* Summary */
    printf("\n");
    printf("==========================================\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);
    printf("==========================================\n");

    descriptors_cleanup();

    return tests_failed > 0 ? 1 : 0;
}
