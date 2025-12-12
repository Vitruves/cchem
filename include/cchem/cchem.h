/**
 * @file cchem.h
 * @brief Main header for cchem chemoinformatics library
 */

#ifndef CCHEM_H
#define CCHEM_H

/* Core canonicalizer components */
#include "canonicalizer/types.h"
#include "canonicalizer/element.h"
#include "canonicalizer/atom.h"
#include "canonicalizer/bond.h"
#include "canonicalizer/molecule.h"
#include "canonicalizer/lexer.h"
#include "canonicalizer/parser.h"
#include "canonicalizer/invariant.h"
#include "canonicalizer/ring_finder.h"
#include "canonicalizer/stereo.h"
#include "canonicalizer/canon.h"
#include "canonicalizer/smiles_writer.h"

/* Utilities */
#include "csv.h"
#include "progress.h"
#include "parallel.h"

/* Version information */
#define CCHEM_VERSION_MAJOR 1
#define CCHEM_VERSION_MINOR 0
#define CCHEM_VERSION_PATCH 0
#define CCHEM_VERSION_STRING "1.0.0"

/* Get version string */
const char* cchem_version(void);

/* Initialize cchem library (call once at startup) */
cchem_status_t cchem_init(void);

/* Cleanup cchem library (call once at shutdown) */
void cchem_cleanup(void);

/**
 * @brief Canonicalize a SMILES string
 *
 * @param smiles Input SMILES string
 * @param error_buf Buffer for error message (can be NULL)
 * @param error_buf_size Size of error buffer
 * @return Canonical SMILES string (caller must free) or NULL on error
 */
char* cchem_canonicalize(const char* smiles, char* error_buf, size_t error_buf_size);

/**
 * @brief Check if two SMILES represent the same molecule
 *
 * @param smiles1 First SMILES string
 * @param smiles2 Second SMILES string
 * @return true if equivalent, false otherwise
 */
bool cchem_smiles_equal(const char* smiles1, const char* smiles2);

/**
 * @brief Validate a SMILES string
 *
 * @param smiles SMILES string to validate
 * @param error_buf Buffer for error message (can be NULL)
 * @param error_buf_size Size of error buffer
 * @return CCHEM_OK if valid, error code otherwise
 */
cchem_status_t cchem_validate_smiles(const char* smiles,
                                     char* error_buf, size_t error_buf_size);

/**
 * @brief Process CSV file with parallel canonicalization
 *
 * @param input_file Input CSV file path
 * @param output_file Output CSV file path
 * @param smiles_column Name of column containing SMILES
 * @param output_column Name of column for canonical SMILES
 * @param num_threads Number of threads to use (0 = auto)
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t cchem_process_csv(const char* input_file,
                                 const char* output_file,
                                 const char* smiles_column,
                                 const char* output_column,
                                 int num_threads);

#endif /* CCHEM_H */
