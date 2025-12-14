/**
 * @file sanitize.h
 * @brief Molecular structure sanitization and normalization
 *
 * Provides comprehensive sanitization operations for molecular structures:
 * - Salt/fragment removal (unsalt)
 * - Aromaticity perception and conversion
 * - Kekulization
 * - Charge neutralization
 * - Structure normalization
 * - Tautomer enumeration
 */

#ifndef CCHEM_CANONICALIZER_SANITIZE_H
#define CCHEM_CANONICALIZER_SANITIZE_H

#include "types.h"
#include "molecule.h"

/* Sanitization operation flags */
typedef enum {
    SANITIZE_NONE           = 0,
    SANITIZE_UNSALT         = (1 << 0),  /* Remove salts, keep largest fragment */
    SANITIZE_AROMATIZE      = (1 << 1),  /* Perceive and apply aromaticity */
    SANITIZE_KEKULIZE       = (1 << 2),  /* Convert aromatic to Kekule form */
    SANITIZE_NEUTRALIZE     = (1 << 3),  /* Neutralize charges */
    SANITIZE_REMOVE_STEREO  = (1 << 4),  /* Remove stereochemistry */
    SANITIZE_REMOVE_ISOTOPES = (1 << 5), /* Remove isotope labels */
    SANITIZE_REMOVE_H       = (1 << 6),  /* Remove explicit hydrogens */
    SANITIZE_NORMALIZE      = (1 << 7),  /* Normalize functional groups */
    SANITIZE_VALIDATE       = (1 << 8),  /* Validate structure */
    SANITIZE_CLEANUP        = (1 << 9),  /* General cleanup */
    SANITIZE_ALL            = 0xFFFF     /* All operations */
} sanitize_flags_t;

/* Predefined sanitization presets */
#define SANITIZE_COMPLETE   (SANITIZE_UNSALT | SANITIZE_AROMATIZE | \
                             SANITIZE_NEUTRALIZE | SANITIZE_NORMALIZE | \
                             SANITIZE_VALIDATE)
#define SANITIZE_STANDARD   (SANITIZE_AROMATIZE | SANITIZE_VALIDATE)
#define SANITIZE_MINIMAL    (SANITIZE_VALIDATE)

/* Sanitization options */
typedef struct {
    sanitize_flags_t flags;         /* Operations to perform */

    /* Unsalt options */
    bool keep_all_organic;          /* Keep all organic fragments, not just largest */
    int min_fragment_atoms;         /* Minimum atoms to keep fragment (default: 3) */
    bool remove_water;              /* Remove water molecules (HOH) */

    /* Neutralization options */
    bool neutralize_acids;          /* Deprotonate carboxylic acids, etc. */
    bool neutralize_bases;          /* Protonate amines, etc. */
    bool preserve_quaternary_n;     /* Keep quaternary nitrogen charges */

    /* Normalization options */
    bool normalize_nitro;           /* Normalize nitro groups */
    bool normalize_sulfoxide;       /* Normalize sulfoxides */
    bool normalize_phosphate;       /* Normalize phosphates */

    /* Hydrogen options */
    bool remove_all_h;              /* Remove all explicit H, not just non-polar */
    bool add_explicit_h;            /* Add explicit hydrogens */

    /* Validation options */
    bool check_valence;             /* Check valence rules */
    bool check_charges;             /* Check charge balance */
} sanitize_options_t;

/* Default sanitization options */
extern const sanitize_options_t SANITIZE_OPTIONS_DEFAULT;

/* Strict sanitization options (more aggressive cleanup) */
extern const sanitize_options_t SANITIZE_OPTIONS_STRICT;

/* Minimal sanitization options */
extern const sanitize_options_t SANITIZE_OPTIONS_MINIMAL;

/**
 * @brief Main sanitization function
 *
 * Applies sanitization operations to a molecule in place.
 * Operations are applied in the following order:
 * 1. Validate (if requested)
 * 2. Unsalt (remove fragments)
 * 3. Normalize functional groups
 * 4. Neutralize charges
 * 5. Handle aromaticity (aromatize or kekulize)
 * 6. Remove/add hydrogens
 * 7. Remove stereo/isotopes
 *
 * @param mol Molecule to sanitize (modified in place)
 * @param options Sanitization options (NULL for defaults)
 * @param error_buf Buffer for error message
 * @param error_buf_size Size of error buffer
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t molecule_sanitize(molecule_t* mol,
                                 const sanitize_options_t* options,
                                 char* error_buf, size_t error_buf_size);

/**
 * @brief Remove salts and keep largest organic fragment
 *
 * Identifies disconnected fragments and keeps only the largest one
 * (by heavy atom count). Common counter-ions are recognized:
 * Na+, K+, Ca2+, Mg2+, Cl-, Br-, I-, etc.
 *
 * @param mol Molecule to modify in place
 * @param options Sanitization options (for unsalt-specific settings)
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_unsalt(molecule_t* mol, const sanitize_options_t* options);

/**
 * @brief Perceive aromaticity and convert bonds to aromatic type
 *
 * Uses Hückel's rule (4n+2 pi electrons) to detect aromatic systems
 * and converts bonds to aromatic type. Also marks atoms as aromatic.
 *
 * @param mol Molecule to modify in place
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_aromatize(molecule_t* mol);

/**
 * @brief Convert aromatic bonds to alternating single/double (Kekule form)
 *
 * Uses backtracking algorithm to assign single/double bonds to aromatic
 * systems while satisfying valence constraints.
 *
 * @param mol Molecule to modify in place
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_kekulize(molecule_t* mol);

/**
 * @brief Neutralize charges in molecule
 *
 * Removes protonation states by:
 * - Deprotonating positively charged groups (R-NH3+ -> R-NH2)
 * - Protonating negatively charged groups (R-O- -> R-OH)
 *
 * @param mol Molecule to modify in place
 * @param options Sanitization options (for neutralize-specific settings)
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_neutralize(molecule_t* mol, const sanitize_options_t* options);

/**
 * @brief Normalize functional group representations
 *
 * Converts functional groups to canonical forms:
 * - Nitro groups: [N+](=O)[O-] or N(=O)=O
 * - Sulfoxides: S(=O)
 * - N-oxides: [N+][O-] -> N=O where appropriate
 *
 * @param mol Molecule to modify in place
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_normalize(molecule_t* mol);

/**
 * @brief Remove explicit hydrogens
 *
 * Converts explicit hydrogens to implicit where chemically valid.
 * Preserves hydrogens that carry stereochemistry or isotope labels
 * unless options specify otherwise.
 *
 * @param mol Molecule to modify in place
 * @param options Sanitization options
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_remove_explicit_h(molecule_t* mol,
                                          const sanitize_options_t* options);

/**
 * @brief Add explicit hydrogens
 *
 * Converts implicit hydrogens to explicit atoms.
 *
 * @param mol Molecule to modify in place
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_add_explicit_h(molecule_t* mol);

/**
 * @brief Remove stereochemistry information
 *
 * Clears all stereochemistry:
 * - Tetrahedral chirality (@, @@)
 * - Double bond E/Z (/, \)
 *
 * @param mol Molecule to modify in place
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_remove_stereo(molecule_t* mol);

/**
 * @brief Remove isotope labels
 *
 * Clears all isotope mass labels, converting to natural abundance.
 *
 * @param mol Molecule to modify in place
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_remove_isotopes(molecule_t* mol);

/**
 * @brief Get fragment count information
 *
 * @param mol Molecule to analyze
 * @param num_fragments Output: number of fragments
 * @param largest_fragment_idx Output: index of largest fragment (by heavy atoms)
 * @param largest_fragment_atoms Output: number of heavy atoms in largest fragment
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_get_fragment_info(molecule_t* mol,
                                          int* num_fragments,
                                          int* largest_fragment_idx,
                                          int* largest_fragment_atoms);

/**
 * @brief Keep only specified fragment
 *
 * Removes all atoms not belonging to the specified fragment.
 *
 * @param mol Molecule to modify in place
 * @param fragment_idx Fragment index to keep (from fragment_ids)
 * @return CCHEM_OK on success
 */
cchem_status_t molecule_keep_fragment(molecule_t* mol, int fragment_idx);

/**
 * @brief Check if fragment is a common salt/counter-ion
 *
 * Recognizes common inorganic counter-ions:
 * Na+, K+, Li+, Ca2+, Mg2+, Zn2+, Fe2+, Fe3+,
 * Cl-, Br-, I-, F-, NO3-, SO4-, PO4-, etc.
 *
 * @param mol Molecule containing the fragment
 * @param fragment_idx Fragment index to check
 * @return true if fragment is a recognized salt
 */
bool molecule_fragment_is_salt(const molecule_t* mol, int fragment_idx);

/**
 * @brief Check if fragment is organic
 *
 * A fragment is organic if it contains carbon.
 *
 * @param mol Molecule containing the fragment
 * @param fragment_idx Fragment index to check
 * @return true if fragment contains carbon
 */
bool molecule_fragment_is_organic(const molecule_t* mol, int fragment_idx);

/* =========================================================================
 * Tautomer Enumeration
 * ========================================================================= */

/* Maximum tautomers to enumerate */
#define MAX_TAUTOMERS 100

/* Tautomer enumeration options */
typedef struct {
    int max_tautomers;              /* Maximum tautomers to generate */
    bool include_input;             /* Include input structure in results */
    bool canonical_only;            /* Only return canonical tautomer */
    bool consider_rings;            /* Consider ring-chain tautomerism */
    bool keto_enol;                 /* Enumerate keto-enol tautomers */
    bool amide_imidic;              /* Enumerate amide-imidic acid tautomers */
    bool lactam_lactim;             /* Enumerate lactam-lactim tautomers */
    bool nitroso_oxime;             /* Enumerate nitroso-oxime tautomers */
    bool phosphate;                 /* Enumerate phosphate tautomers */
} tautomer_options_t;

/* Default tautomer options */
extern const tautomer_options_t TAUTOMER_OPTIONS_DEFAULT;

/* Tautomer result */
typedef struct {
    char** smiles;                  /* Array of tautomer SMILES */
    double* scores;                 /* Stability scores (lower = more stable) */
    int num_tautomers;              /* Number of tautomers generated */
    int canonical_idx;              /* Index of canonical tautomer */
} tautomer_result_t;

/**
 * @brief Enumerate tautomers for a molecule
 *
 * Generates tautomeric forms by considering:
 * - 1,3 hydrogen shifts (keto-enol, amide-imidic)
 * - 1,5 hydrogen shifts (ring systems)
 * - Prototropic tautomerism
 *
 * @param mol Input molecule
 * @param options Tautomer enumeration options (NULL for defaults)
 * @param result Output: tautomer result (caller must free with tautomer_result_free)
 * @return CCHEM_OK on success
 */
cchem_status_t tautomer_enumerate(const molecule_t* mol,
                                  const tautomer_options_t* options,
                                  tautomer_result_t* result);

/**
 * @brief Get canonical tautomer
 *
 * Returns the most stable/canonical tautomeric form.
 *
 * @param smiles Input SMILES
 * @param options Tautomer options (NULL for defaults)
 * @param error_buf Error buffer
 * @param error_buf_size Size of error buffer
 * @return Canonical tautomer SMILES (caller must free), or NULL on error
 */
char* tautomer_canonical(const char* smiles,
                         const tautomer_options_t* options,
                         char* error_buf, size_t error_buf_size);

/**
 * @brief Free tautomer result
 */
void tautomer_result_free(tautomer_result_t* result);

/* =========================================================================
 * High-level Convenience Functions
 * ========================================================================= */

/**
 * @brief Sanitize SMILES string
 *
 * Parses SMILES, applies sanitization, returns sanitized canonical SMILES.
 *
 * @param smiles Input SMILES
 * @param flags Sanitization flags
 * @param error_buf Error buffer
 * @param error_buf_size Size of error buffer
 * @return Sanitized SMILES (caller must free), or NULL on error
 */
char* smiles_sanitize(const char* smiles, sanitize_flags_t flags,
                      char* error_buf, size_t error_buf_size);

/**
 * @brief Parse sanitization flags from string
 *
 * Parses comma-separated list of sanitization operation names:
 * "complete", "unsalt", "aromatize", "kekulize", "neutralize",
 * "remove-stereo", "remove-isotopes", "remove-h", "normalize", "validate"
 *
 * @param str Input string
 * @param flags Output: parsed flags
 * @return CCHEM_OK on success, CCHEM_ERROR_INVALID_INPUT if unknown flag
 */
cchem_status_t sanitize_parse_flags(const char* str, sanitize_flags_t* flags);

/**
 * @brief Get string representation of sanitization flags
 *
 * @param flags Sanitization flags
 * @param buf Output buffer
 * @param buf_size Size of buffer
 */
void sanitize_flags_to_string(sanitize_flags_t flags, char* buf, size_t buf_size);

#endif /* CCHEM_CANONICALIZER_SANITIZE_H */
