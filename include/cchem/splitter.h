/**
 * @file splitter.h
 * @brief Dataset splitting utilities for train/test/validation splits
 *
 * Supports both random and scaffold-based splitting for machine learning
 * dataset preparation. Scaffold splitting ensures that molecules with the
 * same Murcko scaffold are grouped together in splits.
 */

#ifndef CCHEM_SPLITTER_H
#define CCHEM_SPLITTER_H

#include "canonicalizer/types.h"
#include "canonicalizer/molecule.h"
#include <stdbool.h>
#include <stddef.h>

/* Maximum number of output splits */
#define SPLITTER_MAX_SPLITS 10

/* Maximum number of molecules to split */
#define SPLITTER_MAX_MOLECULES 10000000

/* Splitting method */
typedef enum {
    SPLIT_METHOD_RANDOM,    /* Random shuffling */
    SPLIT_METHOD_SCAFFOLD   /* Murcko scaffold-based grouping */
} split_method_t;

/* Split options */
typedef struct {
    split_method_t method;           /* Splitting method */
    int num_splits;                  /* Number of output splits */
    double ratios[SPLITTER_MAX_SPLITS]; /* Split ratios (should sum to 1.0) */
    unsigned int seed;               /* Random seed (0 = use time) */
    int num_threads;                 /* Number of threads (0 = auto) */
    bool verbose;                    /* Print progress info */
    bool stratified;                 /* For scaffold: distribute each scaffold across all splits */
} split_options_t;

/* Default split options */
#define SPLIT_OPTIONS_DEFAULT { \
    .method = SPLIT_METHOD_RANDOM, \
    .num_splits = 2, \
    .ratios = {0.8, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, \
    .seed = 0, \
    .num_threads = 0, \
    .verbose = false, \
    .stratified = false \
}

/* Split result for a single molecule */
typedef struct {
    int split_idx;     /* Which split this row belongs to (0-indexed) */
    int original_idx;  /* Original row index in input */
} split_assignment_t;

/* Split result structure */
typedef struct {
    split_assignment_t* assignments;  /* Array of split assignments */
    int num_assignments;              /* Number of assignments */
    int* split_counts;                /* Count per split */
    int num_splits;                   /* Number of splits */
    int num_scaffolds;                /* Number of unique scaffolds (for scaffold split) */
} split_result_t;

/* Scaffold information for a molecule */
typedef struct {
    char* scaffold_smiles;  /* Canonical SMILES of Murcko scaffold */
    int group_id;           /* Scaffold group ID */
} scaffold_info_t;

/**
 * @brief Parse split ratios from string
 *
 * Parses a comma-separated string of ratios (e.g., "80,10,10" or "0.8,0.1,0.1")
 * and normalizes them to sum to 1.0.
 *
 * @param ratio_str Input ratio string
 * @param ratios Output array of ratios (must be at least SPLITTER_MAX_SPLITS)
 * @param num_ratios Output: number of ratios parsed
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t splitter_parse_ratios(const char* ratio_str, double* ratios, int* num_ratios);

/**
 * @brief Extract Murcko scaffold from molecule
 *
 * Computes the Murcko scaffold (ring systems + linkers) of a molecule.
 * The scaffold is returned as a canonical SMILES string.
 *
 * @param mol Input molecule
 * @param scaffold_buf Buffer for scaffold SMILES
 * @param scaffold_buf_size Size of buffer
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t molecule_get_murcko_scaffold(const molecule_t* mol,
                                            char* scaffold_buf,
                                            size_t scaffold_buf_size);

/**
 * @brief Extract Murcko scaffold from SMILES string
 *
 * @param smiles Input SMILES string
 * @param scaffold_buf Buffer for scaffold SMILES
 * @param scaffold_buf_size Size of buffer
 * @param error_buf Error message buffer (can be NULL)
 * @param error_buf_size Size of error buffer
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t smiles_get_murcko_scaffold(const char* smiles,
                                          char* scaffold_buf,
                                          size_t scaffold_buf_size,
                                          char* error_buf,
                                          size_t error_buf_size);

/**
 * @brief Create split result structure
 *
 * @param num_rows Number of rows to split
 * @param num_splits Number of output splits
 * @return Split result or NULL on error
 */
split_result_t* split_result_create(int num_rows, int num_splits);

/**
 * @brief Free split result
 *
 * @param result Split result to free
 */
void split_result_free(split_result_t* result);

/**
 * @brief Perform random split
 *
 * Randomly shuffles indices and partitions according to ratios.
 *
 * @param num_rows Total number of rows
 * @param options Split options
 * @return Split result or NULL on error
 */
split_result_t* splitter_random_split(int num_rows, const split_options_t* options);

/**
 * @brief Perform scaffold-based split
 *
 * Groups molecules by Murcko scaffold and partitions scaffold groups
 * to maintain structural diversity across splits.
 *
 * @param smiles Array of SMILES strings
 * @param num_molecules Number of molecules
 * @param options Split options
 * @param error_buf Error buffer (can be NULL)
 * @param error_buf_size Error buffer size
 * @return Split result or NULL on error
 */
split_result_t* splitter_scaffold_split(const char** smiles,
                                        int num_molecules,
                                        const split_options_t* options,
                                        char* error_buf,
                                        size_t error_buf_size);

/**
 * @brief Split a CSV file
 *
 * Main entry point for splitting CSV files. Reads input, computes splits,
 * and writes output files.
 *
 * @param input_file Input CSV file path
 * @param output_files Array of output file paths
 * @param num_outputs Number of output files
 * @param smiles_col Name of SMILES column
 * @param options Split options
 * @param error_buf Error buffer
 * @param error_buf_size Error buffer size
 * @return CCHEM_OK on success, error code otherwise
 */
cchem_status_t splitter_split_csv(const char* input_file,
                                  const char** output_files,
                                  int num_outputs,
                                  const char* smiles_col,
                                  const split_options_t* options,
                                  char* error_buf,
                                  size_t error_buf_size);

#endif /* CCHEM_SPLITTER_H */
