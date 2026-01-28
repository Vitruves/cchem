/**
 * @file cchem.h
 * @brief Main header for cchem chemoinformatics library
 *
 * This is the umbrella header that includes all cchem functionality:
 * - Canonicalizer: SMILES parsing, canonicalization, and writing
 * - Descriptors: Molecular descriptor computation
 * - Depictor: 2D/3D molecular structure visualization
 * - Splitter: Dataset splitting for ML workflows
 * - Utilities: CSV processing, threading, progress bars
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
#include "canonicalizer/sanitize.h"

/* Molecular descriptors */
#include "utils/descriptors.h"

/* Molecular depiction (2D/3D rendering) */
#include "depictor/depictor.h"

/* Dataset splitting */
#include "splitter/splitter.h"

/* Utilities */
#include "utils/csv.h"
#include "utils/progress.h"
#include "utils/threading.h"
#include "utils/memory.h"

/* ============================================================================
 * Version Information
 * ============================================================================ */

#define CCHEM_VERSION_MAJOR 1
#define CCHEM_VERSION_MINOR 1
#define CCHEM_VERSION_PATCH 0
#define CCHEM_VERSION_STRING "1.1.0"

/** @brief Get version string */
const char* cchem_version(void);

/** @brief Initialize cchem library (call once at startup) */
cchem_status_t cchem_init(void);

/** @brief Cleanup cchem library (call once at shutdown) */
void cchem_cleanup(void);

/* ============================================================================
 * Quick Reference - Main Functions by Module
 * ============================================================================
 *
 * CANONICALIZER (canonicalizer/canon.h, canonicalizer/parser.h):
 *   smiles_canonicalize()     - Canonicalize SMILES string
 *   smiles_to_molecule()      - Parse SMILES to molecule
 *   smiles_validate()         - Validate SMILES syntax
 *   smiles_are_equivalent()   - Check if two SMILES are equivalent
 *   molecule_to_smiles()      - Convert molecule to SMILES
 *
 * DESCRIPTORS (utils/descriptors.h):
 *   descriptors_init()        - Initialize descriptor registry
 *   descriptor_get()          - Get descriptor by name
 *   descriptors_compute_all() - Compute all descriptors for molecule
 *
 * DEPICTOR (depictor/depictor.h):
 *   depict_smiles()           - Render SMILES to image file
 *   depict_molecule()         - Render molecule to image file
 *
 * SPLITTER (splitter/splitter.h):
 *   splitter_split_csv()      - Split CSV dataset
 *
 * CSV PROCESSING (utils/csv.h):
 *   csv_batch_context_create() - Create batch processing context
 *   csv_batch_canonicalize()   - Batch canonicalize CSV file
 */

#endif /* CCHEM_H */
