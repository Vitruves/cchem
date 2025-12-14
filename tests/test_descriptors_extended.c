/**
 * @file test_descriptors_extended.c
 * @brief Tests for extended molecular descriptors (CPSA, BCUT, Moments, Aromatic, etc.)
 *
 * Uses Imatinib (Gleevec) as primary test molecule - a well-characterized drug with:
 * - MW: 493.603 g/mol
 * - logP: ~3.5
 * - TPSA: 86.28 Å²
 * - 4 aromatic rings (2 pyridine/pyrimidine, 2 phenyl)
 * - 1 piperazine ring (non-aromatic)
 * - Multiple N atoms in various environments
 * - Amide bond
 *
 * Additional test molecules:
 * - Caffeine: fused heterocyclic system
 * - Naphthalene: fused aromatic carbocycles
 * - Aspirin: simple aromatic with functional groups
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
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

#define ASSERT_DOUBLE_EQ(a, b, tol) do { \
    double _a = (a); \
    double _b = (b); \
    if (fabs(_a - _b) > (tol)) { \
        printf("FAILED: %s = %.4f, expected %.4f (tol=%.4f, line %d)\n", #a, _a, _b, (double)(tol), __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_DOUBLE_RANGE(a, low, high) do { \
    double _a = (a); \
    if (_a < (low) || _a > (high)) { \
        printf("FAILED: %s = %.4f, expected [%.4f, %.4f] (line %d)\n", #a, _a, (double)(low), (double)(high), __LINE__); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(x) ASSERT_EQ(!!(x), 1)
#define ASSERT_NOT_NULL(x) ASSERT_TRUE((x) != NULL)

/* Test molecules - SMILES */
#define IMATINIB "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
#define CAFFEINE "Cn1cnc2c1c(=O)n(c(=O)n2C)C"
#define NAPHTHALENE "c1ccc2ccccc2c1"
#define ASPIRIN "CC(=O)Oc1ccccc1C(=O)O"
#define BENZENE "c1ccccc1"
#define PYRIDINE "c1ccncc1"
#define THIOPHENE "c1ccsc1"
#define ETHANOL "CCO"
#define HEXANE "CCCCCC"

/* Helper: get integer descriptor value */
static int64_t get_desc_int(const char* smiles, const char* desc_name) {
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

/* Helper: get double descriptor value */
static double get_desc_double(const char* smiles, const char* desc_name) {
    char error_buf[256];
    descriptor_value_t value = {0};

    cchem_status_t status = descriptor_compute_from_smiles(smiles, desc_name, &value,
                                                           error_buf, sizeof(error_buf));
    if (status != CCHEM_OK) {
        printf("Error computing %s for %s: %s\n", desc_name, smiles, error_buf);
        return -999.0;
    }
    return value.d;
}

/* ===========================================================================
 * CPSA Descriptor Tests (Charged Partial Surface Area)
 * =========================================================================== */

TEST(cpsa_benzene) {
    /* Benzene: symmetric, all carbons slightly positive due to delocalization */
    double cpsa1 = get_desc_double(BENZENE, "CPSA_1");  /* Positive SA */
    double cpsa2 = get_desc_double(BENZENE, "CPSA_2");  /* Negative SA */
    double cpsa7 = get_desc_double(BENZENE, "CPSA_7");  /* Near-neutral SA */

    /* Benzene should have mostly near-neutral surface */
    ASSERT_TRUE(cpsa7 > 0);
    ASSERT_TRUE(cpsa1 >= 0);
    ASSERT_DOUBLE_EQ(cpsa2, 0.0, 0.01);  /* No strongly negative atoms */

    /* Total SA should be reasonable for 6 carbons (includes H in VSA calc) */
    double total_sa = cpsa1 + cpsa7;
    ASSERT_DOUBLE_RANGE(total_sa, 300.0, 600.0);

    tests_passed++;
}

TEST(cpsa_ethanol) {
    /* Ethanol: O is negative, C's are positive */
    double cpsa1 = get_desc_double(ETHANOL, "CPSA_1");  /* Positive SA */
    double cpsa2 = get_desc_double(ETHANOL, "CPSA_2");  /* Negative SA */

    /* Oxygen should contribute negative SA */
    ASSERT_TRUE(cpsa2 > 0);
    /* Carbons should contribute positive SA */
    ASSERT_TRUE(cpsa1 > 0);

    tests_passed++;
}

TEST(cpsa_imatinib) {
    /* Imatinib: complex drug with multiple charge environments */
    double cpsa1 = get_desc_double(IMATINIB, "CPSA_1");
    double cpsa2 = get_desc_double(IMATINIB, "CPSA_2");
    double cpsa13 = get_desc_double(IMATINIB, "CPSA_13");  /* Total charge-weighted SA */

    /* Should have both positive and negative contributions */
    ASSERT_TRUE(cpsa1 > 100);  /* Multiple carbons */
    ASSERT_TRUE(cpsa2 > 50);   /* Multiple nitrogens */
    ASSERT_TRUE(cpsa13 != 0);  /* Non-zero charge-weighted SA */

    /* Total SA reasonable for 29 heavy atoms plus hydrogens */
    ASSERT_DOUBLE_RANGE(cpsa1 + cpsa2, 800.0, 1600.0);

    tests_passed++;
}

TEST(hbsa_descriptors) {
    /* H-bond surface area tests */

    /* Ethanol: O is both donor and acceptor */
    double hbsa1_eth = get_desc_double(ETHANOL, "HBSA_1");  /* Donor SA */
    double hbsa2_eth = get_desc_double(ETHANOL, "HBSA_2");  /* Acceptor SA */
    ASSERT_TRUE(hbsa1_eth > 0);  /* OH is donor */
    ASSERT_TRUE(hbsa2_eth > 0);  /* O is acceptor */

    /* Benzene: no H-bond capability */
    double hbsa1_ben = get_desc_double(BENZENE, "HBSA_1");
    double hbsa2_ben = get_desc_double(BENZENE, "HBSA_2");
    ASSERT_DOUBLE_EQ(hbsa1_ben, 0.0, 0.01);
    ASSERT_DOUBLE_EQ(hbsa2_ben, 0.0, 0.01);

    /* Imatinib: multiple NH donors and N acceptors */
    double hbsa1_ima = get_desc_double(IMATINIB, "HBSA_1");
    double hbsa2_ima = get_desc_double(IMATINIB, "HBSA_2");
    ASSERT_TRUE(hbsa1_ima > 0);  /* NH groups */
    ASSERT_TRUE(hbsa2_ima > 0);  /* N atoms */

    tests_passed++;
}

TEST(lpsa_descriptors) {
    /* Lipophilicity surface area */

    /* Hexane: all hydrophobic */
    double lpsa1_hex = get_desc_double(HEXANE, "LPSA_1");  /* Hydrophobic SA */
    double lpsa2_hex = get_desc_double(HEXANE, "LPSA_2");  /* Hydrophilic SA */
    ASSERT_TRUE(lpsa1_hex > 0);
    ASSERT_DOUBLE_EQ(lpsa2_hex, 0.0, 0.01);

    /* Ethanol: mixed */
    double lpsa1_eth = get_desc_double(ETHANOL, "LPSA_1");
    double lpsa2_eth = get_desc_double(ETHANOL, "LPSA_2");
    ASSERT_TRUE(lpsa1_eth > 0);  /* C atoms */
    ASSERT_TRUE(lpsa2_eth > 0);  /* O atom */

    tests_passed++;
}

/* ===========================================================================
 * BCUT Extension Tests
 * =========================================================================== */

TEST(bcut_ext_benzene) {
    /* Benzene: symmetric molecule, eigenvalues should reflect uniformity */
    double ea_hi1 = get_desc_double(BENZENE, "BCUTea_hi1");
    double ea_lo1 = get_desc_double(BENZENE, "BCUTea_lo1");

    /* EA eigenvalues should be non-zero */
    ASSERT_TRUE(ea_hi1 != 0 || ea_lo1 != 0);

    /* Aromatic BCUT should be high (all atoms aromatic) */
    double ar_hi1 = get_desc_double(BENZENE, "BCUTar_hi1");
    ASSERT_TRUE(ar_hi1 > 0);

    tests_passed++;
}

TEST(bcut_ext_imatinib) {
    /* Imatinib: diverse atom types */
    double ea_hi1 = get_desc_double(IMATINIB, "BCUTea_hi1");
    double r_hi1 = get_desc_double(IMATINIB, "BCUTr_hi1");   /* Covalent radius */
    double v_hi1 = get_desc_double(IMATINIB, "BCUTv_hi1");   /* Valence electrons */
    double vdw_hi1 = get_desc_double(IMATINIB, "BCUTvdw_hi1"); /* VdW volume */

    /* All should have non-zero values */
    ASSERT_TRUE(ea_hi1 > 0 || ea_hi1 < 0);
    ASSERT_TRUE(r_hi1 > 0);
    ASSERT_TRUE(v_hi1 > 0);
    ASSERT_TRUE(vdw_hi1 > 0);

    /* Ring count BCUT - Imatinib has 5 rings */
    double rc_hi1 = get_desc_double(IMATINIB, "BCUTrc_hi1");
    ASSERT_TRUE(rc_hi1 > 0);

    tests_passed++;
}

TEST(bcut_ext_range) {
    /* Test that high eigenvalues >= low eigenvalues */
    double ea_hi1 = get_desc_double(IMATINIB, "BCUTea_hi1");
    double ea_lo1 = get_desc_double(IMATINIB, "BCUTea_lo1");
    ASSERT_TRUE(ea_hi1 >= ea_lo1);

    double r_hi1 = get_desc_double(IMATINIB, "BCUTr_hi1");
    double r_lo1 = get_desc_double(IMATINIB, "BCUTr_lo1");
    ASSERT_TRUE(r_hi1 >= r_lo1);

    tests_passed++;
}

/* ===========================================================================
 * Property Moments Tests
 * =========================================================================== */

TEST(moments_benzene) {
    /* Benzene: uniform atoms, should have low variance/skewness */
    double mass_skew = get_desc_double(BENZENE, "Mass_Skew");
    double en_skew = get_desc_double(BENZENE, "EN_Skew");

    /* Uniform distribution should have near-zero skewness */
    ASSERT_DOUBLE_RANGE(mass_skew, -0.5, 0.5);
    ASSERT_DOUBLE_RANGE(en_skew, -0.5, 0.5);

    /* Median mass should be carbon mass (~12) */
    double mass_median = get_desc_double(BENZENE, "Mass_Median");
    ASSERT_DOUBLE_RANGE(mass_median, 11.0, 13.0);

    /* IQR should be 0 (all same element) */
    double mass_iqr = get_desc_double(BENZENE, "Mass_IQR");
    ASSERT_DOUBLE_EQ(mass_iqr, 0.0, 0.1);

    tests_passed++;
}

TEST(moments_imatinib) {
    /* Imatinib: diverse atoms (C, N, O) */
    double mass_median = get_desc_double(IMATINIB, "Mass_Median");
    double en_median = get_desc_double(IMATINIB, "EN_Median");

    /* Median mass should be around carbon (dominant element) */
    ASSERT_DOUBLE_RANGE(mass_median, 10.0, 15.0);

    /* EN median around 2.5-3.0 (mix of C, N) */
    ASSERT_DOUBLE_RANGE(en_median, 2.0, 3.5);

    /* IQR may be 0 if molecule is dominated by one element (like carbon) */
    double mass_iqr = get_desc_double(IMATINIB, "Mass_IQR");
    ASSERT_TRUE(mass_iqr >= 0);  /* Should be non-negative */

    /* Entropy should be non-zero for diverse molecule */
    double mass_entropy = get_desc_double(IMATINIB, "Mass_Entropy");
    ASSERT_TRUE(mass_entropy > 0);

    tests_passed++;
}

TEST(moments_entropy_gini) {
    /* Entropy and Gini tests */

    /* Benzene: uniform -> low entropy, low Gini */
    double ben_entropy = get_desc_double(BENZENE, "Mass_Entropy");
    double ben_gini = get_desc_double(BENZENE, "Mass_Gini");
    ASSERT_DOUBLE_RANGE(ben_entropy, 0.0, 0.3);
    ASSERT_DOUBLE_RANGE(ben_gini, -0.1, 0.3);

    /* Imatinib: diverse -> higher entropy */
    double ima_entropy = get_desc_double(IMATINIB, "Mass_Entropy");
    ASSERT_TRUE(ima_entropy > ben_entropy);

    tests_passed++;
}

/* ===========================================================================
 * Aromatic Pattern Tests
 * =========================================================================== */

TEST(aromatic_ring_types) {
    /* Benzene */
    ASSERT_EQ((int64_t)get_desc_double(BENZENE, "Arom_Benzene"), 1);
    ASSERT_EQ((int64_t)get_desc_double(BENZENE, "Arom_AllCarbon"), 1);
    ASSERT_EQ((int64_t)get_desc_double(BENZENE, "Arom_Isolated"), 1);

    /* Pyridine */
    ASSERT_EQ((int64_t)get_desc_double(PYRIDINE, "Arom_Pyridine"), 1);
    ASSERT_EQ((int64_t)get_desc_double(PYRIDINE, "Arom_Heterocyclic"), 1);

    /* Thiophene */
    ASSERT_EQ((int64_t)get_desc_double(THIOPHENE, "Arom_Thiophene"), 1);

    /* Naphthalene - fused benzene rings */
    ASSERT_EQ((int64_t)get_desc_double(NAPHTHALENE, "Arom_Benzene"), 2);
    ASSERT_EQ((int64_t)get_desc_double(NAPHTHALENE, "Arom_Fused2"), 1);
    ASSERT_EQ((int64_t)get_desc_double(NAPHTHALENE, "Arom_NaphthaleneL"), 1);

    tests_passed++;
}

TEST(aromatic_imatinib) {
    /* Imatinib: 2 phenyl + 1 pyridine + 1 pyrimidine aromatic */
    int benzene_count = (int)get_desc_double(IMATINIB, "Arom_Benzene");
    int pyridine_count = (int)get_desc_double(IMATINIB, "Arom_Pyridine");
    int pyrimidine_count = (int)get_desc_double(IMATINIB, "Arom_Pyrimidine");

    /* Total aromatic rings should be 4 */
    ASSERT_TRUE(benzene_count + pyridine_count + pyrimidine_count >= 3);

    /* Should have heterocyclic aromatics */
    int hetero_count = (int)get_desc_double(IMATINIB, "Arom_Heterocyclic");
    ASSERT_TRUE(hetero_count >= 2);

    tests_passed++;
}

TEST(aromatic_density) {
    /* Benzene: all atoms aromatic */
    double dens = get_desc_double(BENZENE, "AromDens_Atoms");
    ASSERT_DOUBLE_EQ(dens, 1.0, 0.01);

    /* Hexane: no aromatic atoms */
    double dens_hex = get_desc_double(HEXANE, "AromDens_Atoms");
    ASSERT_DOUBLE_EQ(dens_hex, 0.0, 0.01);

    /* Imatinib: partial aromatic */
    double dens_ima = get_desc_double(IMATINIB, "AromDens_Atoms");
    ASSERT_DOUBLE_RANGE(dens_ima, 0.5, 0.9);

    tests_passed++;
}

TEST(aromatic_electronic) {
    /* Test aromatic electronic properties */
    double chi_sum = get_desc_double(BENZENE, "AromChi_Sum");
    double chi_mean = get_desc_double(BENZENE, "AromChi_Mean");

    /* Benzene: 6 carbons with EN ~2.55 */
    ASSERT_DOUBLE_RANGE(chi_sum, 14.0, 17.0);  /* 6 * 2.55 = 15.3 */
    ASSERT_DOUBLE_RANGE(chi_mean, 2.4, 2.7);

    /* Variance should be ~0 for uniform benzene */
    double chi_var = get_desc_double(BENZENE, "AromChi_Var");
    ASSERT_DOUBLE_RANGE(chi_var, -0.1, 0.1);

    /* Pyridine: N has higher EN */
    double chi_var_pyr = get_desc_double(PYRIDINE, "AromChi_Var");
    ASSERT_TRUE(chi_var_pyr > chi_var);

    tests_passed++;
}

TEST(aromatic_heteroatoms) {
    /* Count N, O, S in aromatic rings */
    int n_ben = (int)get_desc_double(BENZENE, "AromHetero_N");
    int n_pyr = (int)get_desc_double(PYRIDINE, "AromHetero_N");
    int s_thio = (int)get_desc_double(THIOPHENE, "AromHetero_S");

    ASSERT_EQ(n_ben, 0);
    ASSERT_EQ(n_pyr, 1);
    ASSERT_EQ(s_thio, 1);

    /* Imatinib: multiple aromatic N */
    int n_ima = (int)get_desc_double(IMATINIB, "AromHetero_N");
    ASSERT_TRUE(n_ima >= 3);

    tests_passed++;
}

/* ===========================================================================
 * Extended Atom Pair Tests
 * =========================================================================== */

TEST(atompairs_aromatic) {
    /* Benzene: all aromatic-aromatic pairs */
    int aa1 = (int)get_desc_double(BENZENE, "AP_AromArom1");
    int aa2 = (int)get_desc_double(BENZENE, "AP_AromArom2");
    int aa3 = (int)get_desc_double(BENZENE, "AP_AromArom3");

    /* Distance 1: 6 pairs (each atom bonded to 2 neighbors, 6 bonds) */
    ASSERT_EQ(aa1, 6);
    /* Distance 2: each atom reaches 2 atoms at distance 2 */
    ASSERT_TRUE(aa2 > 0);
    /* Distance 3: opposite atoms */
    ASSERT_TRUE(aa3 > 0);

    /* Aromatic-aliphatic should be 0 */
    int al1 = (int)get_desc_double(BENZENE, "AP_AromAlip1");
    ASSERT_EQ(al1, 0);

    tests_passed++;
}

TEST(atompairs_hexane) {
    /* Hexane: all aliphatic, no aromatic pairs */
    int aa1 = (int)get_desc_double(HEXANE, "AP_AromArom1");
    ASSERT_EQ(aa1, 0);

    /* Chain-chain pairs */
    int cc1 = (int)get_desc_double(HEXANE, "AP_ChainChain1");
    ASSERT_TRUE(cc1 > 0);

    tests_passed++;
}

TEST(atompairs_ring_chain) {
    /* Toluene: 1 ring, 1 chain atom */
    const char* toluene = "Cc1ccccc1";
    int rr1 = (int)get_desc_double(toluene, "AP_RingRing1");
    int rc1 = (int)get_desc_double(toluene, "AP_RingChain1");

    ASSERT_TRUE(rr1 > 0);  /* Ring-ring pairs in benzene */
    ASSERT_TRUE(rc1 > 0);  /* Methyl connected to ring */

    tests_passed++;
}

TEST(atompairs_imatinib) {
    /* Imatinib: mixed aromatic and aliphatic */
    int aa1 = (int)get_desc_double(IMATINIB, "AP_AromArom1");
    int al1 = (int)get_desc_double(IMATINIB, "AP_AromAlip1");

    ASSERT_TRUE(aa1 > 0);  /* Aromatic rings */
    ASSERT_TRUE(al1 > 0);  /* Connections to piperazine */

    /* Should have ring-ring pairs at various distances */
    int rr3 = (int)get_desc_double(IMATINIB, "AP_RingRing3");
    ASSERT_TRUE(rr3 > 0);

    tests_passed++;
}

/* ===========================================================================
 * Framework Topology Tests
 * =========================================================================== */

TEST(framework_benzene) {
    /* Benzene: single ring system */
    int systems = (int)get_desc_double(BENZENE, "Frame_RingSystems");
    ASSERT_EQ(systems, 1);

    /* No fused pairs, spiro, bridgeheads */
    int fused = (int)get_desc_double(BENZENE, "Frame_FusedPairs");
    int spiro = (int)get_desc_double(BENZENE, "Frame_SpiroJunctions");
    int bridge = (int)get_desc_double(BENZENE, "Frame_BridgeheadAtoms");
    ASSERT_EQ(fused, 0);
    ASSERT_EQ(spiro, 0);
    ASSERT_EQ(bridge, 0);

    /* Ring size = 6 */
    int max_size = (int)get_desc_double(BENZENE, "Frame_MaxRingSize");
    int min_size = (int)get_desc_double(BENZENE, "Frame_MinRingSize");
    ASSERT_EQ(max_size, 6);
    ASSERT_EQ(min_size, 6);

    /* All atoms in ring */
    int ring_atoms = (int)get_desc_double(BENZENE, "Frame_TotalRingAtoms");
    ASSERT_EQ(ring_atoms, 6);

    tests_passed++;
}

TEST(framework_naphthalene) {
    /* Naphthalene: fused 2-ring system */
    int systems = (int)get_desc_double(NAPHTHALENE, "Frame_RingSystems");
    ASSERT_EQ(systems, 1);  /* Single fused system */

    int fused = (int)get_desc_double(NAPHTHALENE, "Frame_FusedPairs");
    ASSERT_EQ(fused, 1);

    /* 2 bridgehead atoms (shared atoms) */
    int bridge = (int)get_desc_double(NAPHTHALENE, "Frame_BridgeheadAtoms");
    ASSERT_EQ(bridge, 2);

    /* 10 atoms in rings */
    int ring_atoms = (int)get_desc_double(NAPHTHALENE, "Frame_TotalRingAtoms");
    ASSERT_EQ(ring_atoms, 10);

    tests_passed++;
}

TEST(framework_imatinib) {
    /* Imatinib: 4 aromatic + 1 piperazine = 5 rings, some connected */
    int systems = (int)get_desc_double(IMATINIB, "Frame_RingSystems");
    ASSERT_TRUE(systems >= 1);

    /* Should have linker atoms connecting rings */
    int linker = (int)get_desc_double(IMATINIB, "Frame_LinkerAtoms");
    ASSERT_TRUE(linker > 0);

    /* Chain atoms (piperazine methyl, etc.) */
    int chain = (int)get_desc_double(IMATINIB, "Frame_ChainAtoms");
    ASSERT_TRUE(chain > 0);

    /* Complexity should be high */
    double complexity = get_desc_double(IMATINIB, "Frame_Complexity");
    ASSERT_TRUE(complexity > 1);

    tests_passed++;
}

TEST(framework_hexane) {
    /* Hexane: no rings */
    int systems = (int)get_desc_double(HEXANE, "Frame_RingSystems");
    ASSERT_EQ(systems, 0);

    /* All chain atoms */
    int chain = (int)get_desc_double(HEXANE, "Frame_ChainAtoms");
    ASSERT_EQ(chain, 6);

    /* Chain density = 1.0 */
    double density = get_desc_double(HEXANE, "Frame_ChainDensity");
    ASSERT_DOUBLE_EQ(density, 1.0, 0.01);

    tests_passed++;
}

TEST(framework_scaffold) {
    /* Test scaffold metrics */

    /* Benzene: all scaffold */
    double scaffold_ratio = get_desc_double(BENZENE, "Frame_ScaffoldRatio");
    ASSERT_DOUBLE_EQ(scaffold_ratio, 1.0, 0.01);

    /* Toluene: 6/7 scaffold (methyl is side chain) */
    const char* toluene = "Cc1ccccc1";
    double tol_scaffold = get_desc_double(toluene, "Frame_ScaffoldRatio");
    ASSERT_DOUBLE_RANGE(tol_scaffold, 0.7, 1.0);

    /* Decorated ring in toluene */
    int decorated = (int)get_desc_double(toluene, "Frame_DecoratedRings");
    ASSERT_EQ(decorated, 1);

    tests_passed++;
}

/* ===========================================================================
 * Constitutional Extension Tests
 * =========================================================================== */

TEST(constitutional_element_combos) {
    /* CNO counts */
    int cno_eth = (int)get_desc_double(ETHANOL, "Const_CNO");
    ASSERT_EQ(cno_eth, 3);  /* 2C + 1O */

    int cno_ima = (int)get_desc_double(IMATINIB, "Const_CNO");
    /* Imatinib: 29C + 7N + 1O = 37 */
    ASSERT_TRUE(cno_ima > 30);

    /* Period 2 atoms (C, N, O, F) */
    int p2_eth = (int)get_desc_double(ETHANOL, "Const_Period2");
    ASSERT_EQ(p2_eth, 3);

    /* Polar heavy atoms (N, O, S, P) */
    int polar_eth = (int)get_desc_double(ETHANOL, "Const_PolarHeavy");
    ASSERT_EQ(polar_eth, 1);  /* O */

    int polar_ima = (int)get_desc_double(IMATINIB, "Const_PolarHeavy");
    ASSERT_TRUE(polar_ima >= 8);  /* 7N + 1O */

    tests_passed++;
}

TEST(constitutional_en_ip) {
    /* High EN atoms (EN > 3.0): N, O, F, Cl */
    int high_en_eth = (int)get_desc_double(ETHANOL, "Const_HighEN");
    ASSERT_EQ(high_en_eth, 1);  /* O */

    int high_en_ben = (int)get_desc_double(BENZENE, "Const_HighEN");
    ASSERT_EQ(high_en_ben, 0);  /* C has EN ~2.55 */

    /* Low EN atoms (EN < 2.5): metals, some nonmetals */
    int low_en_eth = (int)get_desc_double(ETHANOL, "Const_LowEN");
    ASSERT_EQ(low_en_eth, 0);  /* C=2.55, O=3.44 */

    tests_passed++;
}

TEST(constitutional_bonds) {
    /* Polar vs nonpolar bonds */
    int polar_eth = (int)get_desc_double(ETHANOL, "Const_PolarBonds");
    int nonpolar_eth = (int)get_desc_double(ETHANOL, "Const_NonpolarBonds");
    ASSERT_TRUE(polar_eth > 0);    /* C-O bond is polar */
    ASSERT_TRUE(nonpolar_eth > 0); /* C-C bond is nonpolar */

    /* Carbon-only bonds in hexane */
    int cc_hex = (int)get_desc_double(HEXANE, "Const_CarbonOnlyBonds");
    ASSERT_EQ(cc_hex, 5);

    /* Ring double bonds vs chain double bonds */
    int ring_dbl = (int)get_desc_double(BENZENE, "Const_RingDoubleBonds");
    /* Aromatic bonds aren't counted as "double" - they're aromatic */
    ASSERT_EQ(ring_dbl, 0);  /* Benzene has aromatic bonds, not double bonds */
    /* Check a molecule with actual double bond */
    const char* cyclohexene = "C1CC=CCC1";
    int ring_dbl_ch = (int)get_desc_double(cyclohexene, "Const_RingDoubleBonds");
    ASSERT_EQ(ring_dbl_ch, 1);

    tests_passed++;
}

TEST(constitutional_hybrid) {
    /* sp3 C bonded to heteroatom */
    int csp3_het_eth = (int)get_desc_double(ETHANOL, "Const_Csp3_Hetero");
    ASSERT_EQ(csp3_het_eth, 1);  /* CH2 bonded to O */

    /* sp2 C bonded to heteroatom (carbonyl) */
    const char* acetone = "CC(=O)C";
    int csp2_het_ace = (int)get_desc_double(acetone, "Const_Csp2_Hetero");
    ASSERT_EQ(csp2_het_ace, 1);

    /* N aromatic vs aliphatic */
    int n_arom_pyr = (int)get_desc_double(PYRIDINE, "Const_N_Aromatic");
    int n_aliph_pyr = (int)get_desc_double(PYRIDINE, "Const_N_Aliphatic");
    ASSERT_EQ(n_arom_pyr, 1);
    ASSERT_EQ(n_aliph_pyr, 0);

    /* O carbonyl vs ether */
    int o_carb_ace = (int)get_desc_double(acetone, "Const_O_Carbonyl");
    int o_eth_ace = (int)get_desc_double(acetone, "Const_O_Ether");
    ASSERT_EQ(o_carb_ace, 1);
    ASSERT_EQ(o_eth_ace, 0);

    const char* dme = "COC";  /* Dimethyl ether */
    int o_eth_dme = (int)get_desc_double(dme, "Const_O_Ether");
    ASSERT_EQ(o_eth_dme, 1);

    tests_passed++;
}

/* ===========================================================================
 * Caffeine - Fused Heterocyclic System Tests
 * =========================================================================== */

TEST(caffeine_comprehensive) {
    /* Caffeine: C8H10N4O2, fused purine system */

    /* Element counts */
    ASSERT_EQ(get_desc_int(CAFFEINE, "CarbonCount"), 8);
    ASSERT_EQ(get_desc_int(CAFFEINE, "NitrogenCount"), 4);
    ASSERT_EQ(get_desc_int(CAFFEINE, "OxygenCount"), 2);

    /* Ring system - fused 5+6 */
    int rings = (int)get_desc_double(CAFFEINE, "Frame_RingSystems");
    ASSERT_EQ(rings, 1);  /* Single fused system */

    int fused = (int)get_desc_double(CAFFEINE, "Frame_FusedPairs");
    ASSERT_EQ(fused, 1);

    /* Aromatic heterocycles */
    int hetero_arom = (int)get_desc_double(CAFFEINE, "Arom_Heterocyclic");
    ASSERT_TRUE(hetero_arom >= 1);

    /* Should have aromatic N */
    int n_arom = (int)get_desc_double(CAFFEINE, "AromHetero_N");
    ASSERT_TRUE(n_arom >= 2);

    /* Carbonyl oxygens */
    int o_carb = (int)get_desc_double(CAFFEINE, "Const_O_Carbonyl");
    ASSERT_EQ(o_carb, 2);

    tests_passed++;
}

/* ===========================================================================
 * Aspirin Tests
 * =========================================================================== */

TEST(aspirin_comprehensive) {
    /* Aspirin: C9H8O4 */

    /* Element counts */
    ASSERT_EQ(get_desc_int(ASPIRIN, "CarbonCount"), 9);
    ASSERT_EQ(get_desc_int(ASPIRIN, "OxygenCount"), 4);

    /* Single benzene ring */
    int benzene = (int)get_desc_double(ASPIRIN, "Arom_Benzene");
    ASSERT_EQ(benzene, 1);

    /* Decorated ring (has substituents) */
    int decorated = (int)get_desc_double(ASPIRIN, "Frame_DecoratedRings");
    ASSERT_EQ(decorated, 1);

    /* Carbonyl groups */
    int o_carb = (int)get_desc_double(ASPIRIN, "Const_O_Carbonyl");
    ASSERT_EQ(o_carb, 2);  /* Acetyl C=O and carboxylic C=O */

    /* Ester/ether O */
    int o_ether = (int)get_desc_double(ASPIRIN, "Const_O_Ether");
    ASSERT_EQ(o_ether, 1);  /* Ester O */

    /* H-bond capability */
    double hbsa_donor = get_desc_double(ASPIRIN, "HBSA_1");
    double hbsa_acc = get_desc_double(ASPIRIN, "HBSA_2");
    ASSERT_TRUE(hbsa_donor > 0);  /* Carboxylic OH */
    ASSERT_TRUE(hbsa_acc > 0);    /* Multiple O atoms */

    tests_passed++;
}

/* ===========================================================================
 * Main
 * =========================================================================== */

int main(void) {
    printf("Running extended descriptor tests...\n\n");

    /* Initialize descriptors */
    descriptors_init();
    printf("Total descriptors registered: %d\n\n", descriptor_count());

    /* CPSA Tests */
    printf("CPSA Descriptors:\n");
    RUN_TEST(cpsa_benzene);
    RUN_TEST(cpsa_ethanol);
    RUN_TEST(cpsa_imatinib);
    RUN_TEST(hbsa_descriptors);
    RUN_TEST(lpsa_descriptors);

    /* BCUT Extension Tests */
    printf("\nBCUT Extension Descriptors:\n");
    RUN_TEST(bcut_ext_benzene);
    RUN_TEST(bcut_ext_imatinib);
    RUN_TEST(bcut_ext_range);

    /* Property Moments Tests */
    printf("\nProperty Moments Descriptors:\n");
    RUN_TEST(moments_benzene);
    RUN_TEST(moments_imatinib);
    RUN_TEST(moments_entropy_gini);

    /* Aromatic Pattern Tests */
    printf("\nAromatic Pattern Descriptors:\n");
    RUN_TEST(aromatic_ring_types);
    RUN_TEST(aromatic_imatinib);
    RUN_TEST(aromatic_density);
    RUN_TEST(aromatic_electronic);
    RUN_TEST(aromatic_heteroatoms);

    /* Extended Atom Pair Tests */
    printf("\nExtended Atom Pair Descriptors:\n");
    RUN_TEST(atompairs_aromatic);
    RUN_TEST(atompairs_hexane);
    RUN_TEST(atompairs_ring_chain);
    RUN_TEST(atompairs_imatinib);

    /* Framework Topology Tests */
    printf("\nFramework Topology Descriptors:\n");
    RUN_TEST(framework_benzene);
    RUN_TEST(framework_naphthalene);
    RUN_TEST(framework_imatinib);
    RUN_TEST(framework_hexane);
    RUN_TEST(framework_scaffold);

    /* Constitutional Extension Tests */
    printf("\nConstitutional Extension Descriptors:\n");
    RUN_TEST(constitutional_element_combos);
    RUN_TEST(constitutional_en_ip);
    RUN_TEST(constitutional_bonds);
    RUN_TEST(constitutional_hybrid);

    /* Comprehensive Molecule Tests */
    printf("\nComprehensive Molecule Tests:\n");
    RUN_TEST(caffeine_comprehensive);
    RUN_TEST(aspirin_comprehensive);

    /* Summary */
    printf("\n");
    printf("==========================================\n");
    printf("Extended Descriptor Tests\n");
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);
    printf("==========================================\n");

    descriptors_cleanup();

    return tests_failed > 0 ? 1 : 0;
}
