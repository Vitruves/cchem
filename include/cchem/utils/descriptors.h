/**
 * @file descriptors.h
 * @brief Modular molecular descriptor computation system
 *
 * This module provides a framework for computing molecular descriptors
 * from canonical SMILES. Descriptors are organized into categories
 * (counts, ratios, etc.) and can be computed individually or in batches.
 *
 * All descriptors are computed from parsed molecule structures for speed.
 */

#ifndef CCHEM_DESCRIPTORS_H
#define CCHEM_DESCRIPTORS_H

#include <stddef.h>
#include <stdbool.h>
#include "cchem/canonicalizer/types.h"
#include "cchem/canonicalizer/molecule.h"

/* Maximum number of registered descriptors */
#define MAX_DESCRIPTORS 2048

/* Maximum descriptor name length */
#define MAX_DESCRIPTOR_NAME 64

/* Hash table size for fast descriptor lookup (must be power of 2) */
#define DESC_HASH_TABLE_SIZE 4096
#define DESC_HASH_TABLE_MASK (DESC_HASH_TABLE_SIZE - 1)

/* Descriptor categories */
typedef enum {
    DESC_CATEGORY_COUNTS = 0,
    DESC_CATEGORY_RATIOS,
    DESC_CATEGORY_PROPERTIES,
    DESC_CATEGORY_ELECTRONIC,
    DESC_CATEGORY_STERIC,
    DESC_CATEGORY_ENERGETIC,
    DESC_CATEGORY_BITSTRING,
    DESC_CATEGORY_FINGERPRINT,
    DESC_CATEGORY_CUSTOM,
    DESC_CATEGORY_CONSTITUTIONAL,
    DESC_CATEGORY_COUNT
} descriptor_category_t;

/* Descriptor value (can be int or double) */
typedef union {
    int64_t i;
    double d;
} descriptor_value_t;

/* Descriptor value type */
typedef enum {
    DESC_VALUE_INT = 0,
    DESC_VALUE_DOUBLE
} descriptor_value_type_t;

/* Forward declaration */
typedef struct descriptor_def descriptor_def_t;

/**
 * Descriptor compute function type
 * @param mol Parsed molecule structure
 * @param value Output value
 * @return CCHEM_OK on success, error code otherwise
 */
typedef cchem_status_t (*descriptor_compute_fn)(const molecule_t* mol, descriptor_value_t* value);

/**
 * Batch descriptor compute function type (optimized for computing multiple descriptors)
 * @param mol Parsed molecule structure
 * @param descriptors Array of descriptor definitions to compute
 * @param num_descriptors Number of descriptors
 * @param values Output values array (same order as descriptors)
 * @return CCHEM_OK on success, error code otherwise
 */
typedef cchem_status_t (*descriptor_batch_compute_fn)(const molecule_t* mol,
                                                       const descriptor_def_t** descriptors,
                                                       int num_descriptors,
                                                       descriptor_value_t* values);

/* Descriptor definition */
struct descriptor_def {
    char name[MAX_DESCRIPTOR_NAME];           /* CamelCase name (e.g., "CarbonCount") */
    char name_lower[MAX_DESCRIPTOR_NAME];     /* Pre-computed lowercase for fast comparison */
    char description[128];                     /* Human-readable description */
    descriptor_category_t category;            /* Category */
    descriptor_value_type_t value_type;        /* Return type */
    descriptor_compute_fn compute;             /* Compute function */
    void* user_data;                           /* Optional user data */
    bool registered;                           /* Is registered in registry */
};

/* Hash table entry for fast descriptor lookup */
typedef struct {
    int16_t index;      /* Index into descriptors array, -1 if empty */
    int16_t next;       /* Next index in chain for collision resolution, -1 if none */
} desc_hash_entry_t;

/* Descriptor registry */
typedef struct {
    descriptor_def_t descriptors[MAX_DESCRIPTORS];
    int num_descriptors;
    bool initialized;

    /* Category-specific batch compute functions for optimization */
    descriptor_batch_compute_fn batch_compute[DESC_CATEGORY_COUNT];

    /* Hash table for O(1) descriptor lookup by name */
    desc_hash_entry_t hash_table[DESC_HASH_TABLE_SIZE];
} descriptor_registry_t;

/* Global registry */
extern descriptor_registry_t g_descriptor_registry;

/* ============================================================================
 * Registry Management
 * ============================================================================ */

/**
 * Initialize the descriptor registry and register all built-in descriptors
 */
void descriptors_init(void);

/**
 * Cleanup the descriptor registry
 */
void descriptors_cleanup(void);

/**
 * Register a new descriptor
 * @param def Descriptor definition (copied into registry)
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t descriptor_register(const descriptor_def_t* def);

/**
 * Get descriptor by name (case-insensitive)
 * @param name Descriptor name
 * @return Pointer to descriptor definition, or NULL if not found
 */
const descriptor_def_t* descriptor_get(const char* name);

/**
 * Get all descriptors in a category
 * @param category Category to filter by
 * @param out_descriptors Output array of descriptor pointers
 * @param max_descriptors Maximum number of descriptors to return
 * @return Number of descriptors found
 */
int descriptor_get_by_category(descriptor_category_t category,
                               const descriptor_def_t** out_descriptors,
                               int max_descriptors);

/**
 * Get all registered descriptors
 * @param out_descriptors Output array of descriptor pointers
 * @param max_descriptors Maximum number of descriptors to return
 * @return Number of descriptors found
 */
int descriptor_get_all(const descriptor_def_t** out_descriptors, int max_descriptors);

/**
 * Get number of registered descriptors
 */
int descriptor_count(void);

/**
 * List all descriptor names (for help/documentation)
 */
void descriptor_list_all(void);

/* ============================================================================
 * Cache Management
 * ============================================================================ */

/**
 * Get current cache generation counter
 * Used by descriptor modules to detect when caches should be invalidated
 * @return Current cache generation number
 */
uint64_t descriptor_cache_generation(void);

/* ============================================================================
 * Descriptor Computation
 * ============================================================================ */

/**
 * Compute a single descriptor from a molecule
 * @param mol Parsed molecule
 * @param name Descriptor name
 * @param value Output value
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t descriptor_compute(const molecule_t* mol,
                                  const char* name,
                                  descriptor_value_t* value);

/**
 * Compute a single descriptor from SMILES (parses internally)
 * @param smiles SMILES string (will be canonicalized first)
 * @param name Descriptor name
 * @param value Output value
 * @param error_buf Buffer for error message
 * @param error_buf_size Size of error buffer
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t descriptor_compute_from_smiles(const char* smiles,
                                              const char* name,
                                              descriptor_value_t* value,
                                              char* error_buf,
                                              size_t error_buf_size);

/**
 * Compute multiple descriptors from a molecule (optimized batch)
 * @param mol Parsed molecule
 * @param names Array of descriptor names
 * @param num_descriptors Number of descriptors
 * @param values Output values array
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t descriptors_compute_batch(const molecule_t* mol,
                                         const char** names,
                                         int num_descriptors,
                                         descriptor_value_t* values);

/**
 * Compute all descriptors from a molecule
 * @param mol Parsed molecule
 * @param out_names Output array for descriptor names (can be NULL)
 * @param out_values Output array for values
 * @param max_descriptors Maximum number of descriptors
 * @return Number of descriptors computed, or -1 on error
 */
int descriptors_compute_all(const molecule_t* mol,
                            const char** out_names,
                            descriptor_value_t* out_values,
                            int max_descriptors);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Format descriptor value as string
 * @param def Descriptor definition
 * @param value Descriptor value
 * @param buf Output buffer
 * @param buf_size Buffer size
 * @return Pointer to buf
 */
char* descriptor_format_value(const descriptor_def_t* def,
                              const descriptor_value_t* value,
                              char* buf,
                              size_t buf_size);

/**
 * Parse descriptor names from comma-separated string
 * @param input Comma-separated descriptor names or "all"
 * @param out_names Output array of names
 * @param max_names Maximum number of names
 * @return Number of names parsed, or -1 on error
 */
int descriptor_parse_names(const char* input,
                           char** out_names,
                           int max_names);

/**
 * Free parsed names array
 */
void descriptor_free_names(char** names, int num_names);

/**
 * Get category name string
 */
const char* descriptor_category_name(descriptor_category_t category);

/* ============================================================================
 * Built-in Descriptor Registration Functions
 * Called by descriptors_init() to register all built-in descriptors
 * ============================================================================ */

/**
 * Register count-based descriptors (from counts.c)
 * - Element counts (CarbonCount, HydrogenCount, etc.)
 * - AtomCount, BondCount
 * - SingleBondCount, DoubleBondCount, TripleBondCount
 * - RingCount, AromaticRingCount
 */
void descriptors_register_counts(void);

/**
 * Register ratio-based descriptors (from ratios.c)
 * - Elemental stoichiometry ratios
 * - Electronic hybridization ratios
 * - Electronegativity & hardness ratios
 * - Bond dynamics ratios
 * - Charge & ionization ratios
 * - Topological electronic ratios
 */
void descriptors_register_ratios(void);

/**
 * Register LogP descriptor (from logp.c)
 * - WCLogP: Wildman-Crippen LogP
 */
void descriptors_register_logp(void);

/**
 * Register electronic descriptors (from electronic.c)
 * - Gasteiger-Marsili partial charges (PEOE)
 * - Electronegativity-weighted topology
 * - Electrotopological states (E-States)
 * - Polarizability and refractivity
 * - Electronic connectivity descriptors
 */
void descriptors_register_electronic(void);

/**
 * Register steric descriptors (from steric.c)
 * - Volume descriptors (VdW, McGowan)
 * - Surface area descriptors (TPSA)
 * - Volume partitioning
 * - Shape and compactness topology
 */
void descriptors_register_steric(void);

/**
 * Register energetic descriptors (from energetics.c)
 * - Born solvation proxies
 * - Electronic hardness & reactivity (FMO theory)
 * - Hydrophobic effect & cavity formation
 * - Hansen solubility parameters
 * - Linear moments of distribution
 * - Rotational & vibrational entropy
 */
void descriptors_register_energetic(void);

/**
 * Register bitstring descriptors (from bitstring.c)
 * - Electronic case bits (aromaticity)
 * - Vertical bit integration (columnar)
 * - Horizontal bit dynamics (dipoles)
 * - Atomic density (popcounts)
 * - Structural bit-masks
 * - Information entropy
 * - Electronegativity proxies
 */
void descriptors_register_bitstring(void);

/**
 * Register acid-base descriptors (from acidbase.c)
 * - Acidic group counts by pKa range (super, strong, moderate, weak)
 * - Basic group counts by conjugate acid pKa range
 * - Total acidic/basic function counts
 * - Mean acidic/basic potential (weighted pKa)
 */
void descriptors_register_acidbase(void);

/**
 * Register solubility descriptors (from solubility.c)
 * - CLogS: Calculated aqueous solubility (log mol/L)
 */
void descriptors_register_solubility(void);

/**
 * Register topological descriptors (from topology.c)
 * - Kier-Hall connectivity indices (Chi0, Chi1, Chi2Path, Chi3Path, Chi0v, Chi1v)
 * - Zagreb indices (Zagreb1, Zagreb2, HyperZagreb)
 * - Path-based indices (WienerIndex, MeanPathLength, MolecularDiameter, BalabanJ)
 * - Neighborhood descriptors (MeanNeighborDegree, TwoHopReachability)
 * - Ionization topology (IonizableDensity, PolarPeripheryRatio, HydrophobicCoreIndex)
 * - Ring topology (RingSystemCount, FusedRingAtomRatio, BridgeBondCount)
 * - LogD-relevant (ChargeablePathSum, LipophilicChainLength, PolarClusterSize)
 * - Additional indices (NormZagreb1, RandicIndex)
 */
void descriptors_register_topology(void);

/**
 * Register fractional descriptors (from fractional.c)
 * - MW fractions (FcC, FcN, FcO, FcS, FcF, FcCl, FcBr, FcI, FcHalo, FcHetero)
 * - EN fractions (FcPolar, FcApolar, FcENAboveAvg, FcENBelowAvg, FcENHigh)
 * - Bond fractions (FcCSp3, FcCSp2, FcUnpol, FcPol, FcBondC/N/O/S/P)
 * - Structural fractions (FcAromaticAtoms, FcRingAtoms, FcChargedAtoms, etc.)
 * - Physical property fractions (FcSmallR, FcLargeR, FcHighPolz, FcSmallVdW, etc.)
 */
void descriptors_register_fractional(void);

/**
 * Batch compute all fractional descriptors (60 descriptors)
 * Single pass through atoms and bonds.
 */
int descriptors_compute_fractional_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Register hash-based descriptors (from hash.c)
 * - Full SMILES hashes (HashFull0, HashFull1, HashFull2)
 * - N-gram hashes (HashBigram, HashTrigram, HashQuadgram)
 * - Pattern hashes (HashAromatic, HashAliphatic, HashHetero, HashRing, HashBranch)
 * - MinHash signatures (MinHash0-7)
 * - Window statistics (HashWindowMin, HashWindowMax, HashWindowMean)
 * - Order-sensitive (HashPosWeight, HashXorChain, HashTransition)
 * - Ionizable patterns (HashAcidic, HashBasic, HashIonizable)
 */
void descriptors_register_hash(void);

/**
 * Batch compute all hash descriptors (30 descriptors)
 * Single pass through SMILES string using xxHash.
 */
int descriptors_compute_hash_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Register graph-based descriptors (from graph.c)
 * - Connectivity: density, degree statistics
 * - Path-based: diameter, Wiener index, eccentricity
 * - Centrality: betweenness, closeness, eigenvector
 * - Clustering: transitivity, local clustering
 * - Spectral: spectral radius, graph energy
 * - ADME-specific: polar centrality, branching index
 */
void descriptors_register_graph(void);

/**
 * Batch compute all graph descriptors (30 descriptors)
 * Single igraph construction, all metrics computed in one pass.
 */
int descriptors_compute_graph_all(const molecule_t* mol, descriptor_value_t* values);

/* ============================================================================
 * Optimized Batch Computation Functions (thread-safe)
 * These compute all descriptors in a category with minimal redundant work.
 * Each function uses stack-allocated caches for thread safety.
 * ============================================================================ */

/**
 * Batch compute all bitstring descriptors (30 descriptors)
 * Generates SMILES only once, computes all descriptors in single pass.
 */
int descriptors_compute_bitstring_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Batch compute all electronic descriptors (30 descriptors)
 * Computes Gasteiger charges and E-States only once.
 */
int descriptors_compute_electronic_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Batch compute all ratio descriptors (30 descriptors)
 * Collects molecular statistics only once.
 */
int descriptors_compute_ratios_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Batch compute all energetic descriptors (30 descriptors)
 * Collects energetic statistics only once.
 */
int descriptors_compute_energetic_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Batch compute all steric descriptors (20 descriptors)
 * Collects steric statistics only once.
 */
int descriptors_compute_steric_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Batch compute all count descriptors (83 descriptors)
 * Single pass through atoms and bonds.
 */
int descriptors_compute_counts_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Register autocorrelation descriptors (from autocorrelations.c)
 * - Broto-Moreau 2D autocorrelations
 * - ATSm: Mass-weighted (lags 0-8)
 * - ATSv: Volume-weighted (lags 0-8)
 * - ATSe: Electronegativity-weighted (lags 0-8)
 * - ATSp: Polarizability-weighted (lags 0-8)
 * - ATSi: Ionization potential-weighted (lags 0-8)
 * - ATSc: Charge-weighted (lags 0-8)
 */
void descriptors_register_autocorrelations(void);

/**
 * Batch compute all autocorrelation descriptors (54 descriptors)
 * Computes distance matrix once, then all ATS values.
 */
int descriptors_compute_autocorr_all(const molecule_t* mol, descriptor_value_t* values);

/**
 * Register LogD 7.4 descriptor (from logd74.c)
 * - ccLogD74: Distribution coefficient at pH 7.4 (neural network model)
 * Trained on 26,036 molecules: R² = 0.9427, MAE = 0.2261, RMSE = 0.3194
 */
void descriptors_register_logd74(void);

/**
 * Register functional group descriptors (from functional.c)
 * - Carbonyl variants (aldehyde, ketone, acyl halide, anhydride)
 * - Nitrogen groups (imine, oxime, hydrazine, azo, azide, guanidine, etc.)
 * - Sulfur groups (sulfide, disulfide, sulfoxide, sulfonamide, etc.)
 * - Oxygen groups (epoxide, peroxide, acetal, hemiacetal, enol)
 * - Mixed groups (urea, thiourea, carbamate, isocyanate, phosphate)
 * - Heterocyclic scaffolds (pyridine, pyrimidine, pyrrole, furan, thiophene, imidazole, etc.)
 * - Aliphatic rings (cyclopropane, piperidine, piperazine, morpholine, pyrrolidine)
 * - Drug features (lactam, lactone, beta-lactam, vinyl, allyl, benzyl, phenyl)
 */
void descriptors_register_functional(void);

/**
 * Register pharmacophore and molecular complexity descriptors (from pharmacophore.c)
 * - Pharmacophore points (lipophilic, basic N, acidic O, halogen bond donors)
 * - Pharmacophore density and balance metrics
 * - Molecular complexity (Bertz, flexibility, stereocenters)
 * - Drug-likeness (Ro5, lead-likeness, CNS MPO)
 * - Structural diversity metrics
 */
void descriptors_register_pharmacophore(void);

/**
 * Register MQN (Molecular Quantum Numbers) descriptors (from mqn.c)
 * - 42 integer descriptors for chemical space navigation
 * - Atom counts by type (C, F, Cl, Br, I, S, P, acyclic/cyclic N/O)
 * - Bond counts (acyclic/cyclic single/double/triple)
 * - Polarity (HBA, HBD, charges, valence nodes)
 * - Topology (rings by size, aromatic atoms/bonds, rotatable bonds)
 */
void descriptors_register_mqn(void);

/**
 * Register E-State atom type count descriptors (from estate.c)
 * - Counts atoms by electrotopological state type
 * - Carbon types (ssCH3, ssCH2, sssCH, ssssC, dCH2, dsCH, aaCH, etc.)
 * - Nitrogen types (sNH2, ssNH, sssN, aaNH, aaN, tN)
 * - Oxygen types (sOH, ssO, dO, aaO)
 * - Sulfur types (sSH, ssS, dS, aaS, dssS, ddssS)
 * - Halogen types (sF, sCl, sBr, sI)
 */
void descriptors_register_estate(void);

/**
 * Register VSA (Van der Waals Surface Area) descriptors (from vsa.c)
 * - SlogP_VSA (12 bins): Surface area by LogP contribution
 * - SMR_VSA (10 bins): Surface area by molar refractivity
 * - PEOE_VSA (14 bins): Surface area by partial charge
 * - EState_VSA (11 bins): Surface area by E-state
 */
void descriptors_register_vsa(void);

/**
 * Register extended Chi connectivity indices (from chi.c)
 * - Simple Chi (Chi0-Chi4): Path connectivity
 * - Valence Chi (Chi0v-Chi4v): Valence-weighted connectivity
 * - Cluster Chi (Chi3c, Chi4c, Chi4pc): Branching patterns
 * - Kappa shape indices (Kappa1-3)
 */
void descriptors_register_chi(void);

/**
 * Register topological atom pair descriptors (from atompairs.c)
 * - Pairs of atom types (C, N, O, S, Hal) at distances 1-7
 * - Summary descriptors (total pairs, heteroatom pairs)
 */
void descriptors_register_atompairs(void);

/**
 * Register BCUT eigenvalue descriptors (from bcut.c)
 * - Burden matrix eigenvalues weighted by mass, charge, EN, polarizability, etc.
 * - 48 descriptors (8 properties × 6 eigenvalues)
 */
void descriptors_register_bcut(void);

/**
 * Register Zagreb and related topological indices (from zagreb.c)
 * - Zagreb M1, M2 indices and variants
 * - Randic, Harmonic, ABC, GA connectivity indices
 * - Balaban J, Wiener, and other distance-based indices
 * - 24 descriptors
 */
void descriptors_register_zagreb(void);

/**
 * Register Information Content descriptors (from infocontent.c)
 * - Shannon entropy based molecular complexity
 * - IC, SIC, CIC, TIC at various orders
 * - Bonding information content, Bertz complexity
 * - 24 descriptors
 */
void descriptors_register_infocontent(void);

/**
 * Register Walk Count descriptors (from walkcounts.c)
 * - Molecular walk counts (adjacency matrix powers)
 * - Path counts, self-returning walks
 * - 36 descriptors
 */
void descriptors_register_walkcounts(void);

/**
 * Register E-State Sum descriptors (from estate_sums.c)
 * - Sum of E-State indices by atom type
 * - E-State statistics (min, max, mean, range)
 * - 32 descriptors
 */
void descriptors_register_estate_sums(void);

/**
 * Register Extended Topological Atom (ETA) descriptors (from eta.c)
 * - Shape, core, epsilon indices
 * - Composite ETA indices
 * - 24 descriptors
 */
void descriptors_register_eta(void);

/**
 * Register Ring Complexity descriptors (from ringcomplexity.c)
 * - Ring size distribution
 * - Fusion patterns (fused, spiro, bridgehead)
 * - Ring complexity indices
 * - 18 descriptors
 */
void descriptors_register_ringcomplexity(void);

/**
 * Register CPSA descriptors (from cpsa.c)
 * - Charged partial surface area
 * - H-bond surface area
 * - Lipophilicity surface area
 * - Polarizability/EN surface area
 * - IP/EA surface area
 * - 70 descriptors
 */
void descriptors_register_cpsa(void);

/**
 * Register Extended BCUT descriptors (from bcut_ext.c)
 * - Electron affinity, covalent radius, valence electrons
 * - VdW volume, oxidation state, lone pairs
 * - Aromatic flag, ring count
 * - 48 descriptors
 */
void descriptors_register_bcut_ext(void);

/**
 * Register Property Moment descriptors (from moments.c)
 * - Statistical moments for 7 atomic properties
 * - Skewness, kurtosis, median, IQR, entropy, Gini
 * - 42 descriptors
 */
void descriptors_register_moments(void);

/**
 * Register Aromatic Pattern descriptors (from aromatic.c)
 * - Ring classification (benzene, pyridine, fused systems, etc.)
 * - Aromatic density and distribution
 * - Aromatic electronic properties
 * - Aromatic topology
 * - 64 descriptors
 */
void descriptors_register_aromatic(void);

/**
 * Register Extended Atom Pair descriptors (from atompairs_ext.c)
 * - Aromatic-aliphatic pairs at various distances
 * - Degree-based pairs (terminal, branch, core)
 * - Charge-based pairs
 * - Ring-based pairs
 * - 56 descriptors
 */
void descriptors_register_atompairs_ext(void);

/**
 * Register Framework Topology descriptors (from framework.c)
 * - Ring system topology
 * - Chain/linker analysis
 * - Scaffold metrics
 * - 40 descriptors
 */
void descriptors_register_framework(void);

/**
 * Register Constitutional Extension descriptors (from constitutional.c)
 * - Element combination counts
 * - Bond environment counts
 * - Hybrid counts
 * - 34 descriptors
 */
void descriptors_register_constitutional(void);

/**
 * Batch compute all constitutional extension descriptors (34 descriptors)
 * Single pass through atoms and bonds.
 */
int descriptors_compute_constitutional_all(const molecule_t* mol, descriptor_value_t* values);

#endif /* CCHEM_DESCRIPTORS_H */
