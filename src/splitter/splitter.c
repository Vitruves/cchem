/**
 * @file splitter.c
 * @brief Dataset splitting implementation for train/test/validation splits
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#include "cchem/splitter/splitter.h"
#include "cchem/cchem.h"
#include "cchem/utils/csv.h"
#include "cchem/utils/threading.h"
#include "cchem/utils/progress.h"

/* ========================================================================
 * Ratio parsing
 * ======================================================================== */

cchem_status_t splitter_parse_ratios(const char* ratio_str, double* ratios, int* num_ratios) {
    if (!ratio_str || !ratios || !num_ratios) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Parse comma-separated values */
    char* str_copy = strdup(ratio_str);
    if (!str_copy) return CCHEM_ERROR_MEMORY;

    int count = 0;
    double sum = 0.0;

    char* token = strtok(str_copy, ",");
    while (token && count < SPLITTER_MAX_SPLITS) {
        double val = atof(token);
        if (val <= 0) {
            free(str_copy);
            return CCHEM_ERROR_INVALID_INPUT;
        }
        ratios[count++] = val;
        sum += val;
        token = strtok(NULL, ",");
    }

    free(str_copy);

    if (count < 2) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Normalize to sum to 1.0 */
    for (int i = 0; i < count; i++) {
        ratios[i] /= sum;
    }

    /* Zero out remaining slots */
    for (int i = count; i < SPLITTER_MAX_SPLITS; i++) {
        ratios[i] = 0.0;
    }

    *num_ratios = count;
    return CCHEM_OK;
}

/* ========================================================================
 * Murcko scaffold extraction
 * ======================================================================== */

/**
 * @brief Find atoms that are part of ring systems
 *
 * @param mol Input molecule (rings must be computed)
 * @param in_ring Output: boolean array indicating ring membership
 * @return Number of ring atoms
 */
static int find_ring_atoms(const molecule_t* mol, bool* in_ring) {
    int count = 0;
    memset(in_ring, 0, mol->num_atoms * sizeof(bool));

    for (int i = 0; i < mol->num_rings; i++) {
        const ring_t* ring = &mol->rings[i];
        for (int j = 0; j < ring->size; j++) {
            int atom_idx = ring->atoms[j];
            if (!in_ring[atom_idx]) {
                in_ring[atom_idx] = true;
                count++;
            }
        }
    }
    return count;
}

/**
 * @brief BFS to find shortest path between two ring atoms
 *
 * @param mol Molecule
 * @param start Start atom index
 * @param end End atom index
 * @param in_ring Ring membership array
 * @param path Output: path atoms (excluding start and end)
 * @param max_path Maximum path length
 * @return Path length or -1 if no path
 */
static int find_linker_path(const molecule_t* mol, int start, int end,
                           const bool* in_ring, int* path, int max_path) {
    if (start == end) return 0;

    bool* visited = calloc(mol->num_atoms, sizeof(bool));
    int* parent = malloc(mol->num_atoms * sizeof(int));
    int* queue = malloc(mol->num_atoms * sizeof(int));

    if (!visited || !parent || !queue) {
        free(visited);
        free(parent);
        free(queue);
        return -1;
    }

    for (int i = 0; i < mol->num_atoms; i++) {
        parent[i] = -1;
    }

    int q_start = 0, q_end = 0;
    queue[q_end++] = start;
    visited[start] = true;

    bool found = false;

    while (q_start < q_end && !found) {
        int curr = queue[q_start++];
        const atom_t* atom = molecule_get_atom_const(mol, curr);

        for (int i = 0; i < atom->num_neighbors; i++) {
            int neighbor = atom->neighbors[i];

            if (!visited[neighbor]) {
                visited[neighbor] = true;
                parent[neighbor] = curr;

                if (neighbor == end) {
                    found = true;
                    break;
                }

                /* Only continue through non-ring atoms */
                if (!in_ring[neighbor]) {
                    queue[q_end++] = neighbor;
                }
            }
        }
    }

    int path_len = 0;
    if (found) {
        /* Reconstruct path (excluding start and end) */
        int curr = parent[end];
        while (curr != start && curr != -1 && path_len < max_path) {
            path[path_len++] = curr;
            curr = parent[curr];
        }
    }

    free(visited);
    free(parent);
    free(queue);

    return found ? path_len : -1;
}

/**
 * @brief Find all linker atoms connecting ring systems
 *
 * @param mol Input molecule
 * @param in_ring Ring membership array
 * @param is_linker Output: boolean array for linker atoms
 * @return Number of linker atoms
 */
static int find_linker_atoms(const molecule_t* mol, const bool* in_ring, bool* is_linker) {
    memset(is_linker, 0, mol->num_atoms * sizeof(bool));
    int linker_count = 0;

    int* path = malloc(mol->num_atoms * sizeof(int));
    if (!path) return 0;

    /* For each ring atom, check if any non-ring neighbor leads to another ring */
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!in_ring[i]) continue;

        const atom_t* atom = molecule_get_atom_const(mol, i);

        for (int n = 0; n < atom->num_neighbors; n++) {
            int neighbor = atom->neighbors[n];

            /* Skip ring atoms - we want paths through non-ring atoms */
            if (in_ring[neighbor]) continue;
            if (is_linker[neighbor]) continue;

            /* Check if this neighbor connects to another ring atom */
            for (int j = i + 1; j < mol->num_atoms; j++) {
                if (!in_ring[j]) continue;

                int path_len = find_linker_path(mol, i, j, in_ring, path, mol->num_atoms);
                if (path_len > 0) {
                    for (int p = 0; p < path_len; p++) {
                        if (!is_linker[path[p]]) {
                            is_linker[path[p]] = true;
                            linker_count++;
                        }
                    }
                }
            }
        }
    }

    free(path);
    return linker_count;
}

/**
 * @brief Create scaffold molecule from original molecule
 *
 * @param mol Original molecule
 * @param in_scaffold Boolean array indicating scaffold membership
 * @return New molecule containing only scaffold atoms, or NULL on error
 */
static molecule_t* create_scaffold_molecule(const molecule_t* mol, const bool* in_scaffold) {
    molecule_t* scaffold = molecule_create();
    if (!scaffold) return NULL;

    /* Map old atom indices to new indices */
    int* idx_map = malloc(mol->num_atoms * sizeof(int));
    if (!idx_map) {
        molecule_free(scaffold);
        return NULL;
    }

    for (int i = 0; i < mol->num_atoms; i++) {
        idx_map[i] = -1;
    }

    /* Add scaffold atoms */
    int new_idx = 0;
    for (int i = 0; i < mol->num_atoms; i++) {
        if (!in_scaffold[i]) continue;

        const atom_t* old_atom = molecule_get_atom_const(mol, i);

        atom_t new_atom = *old_atom;
        new_atom.num_neighbors = 0;

        /* Reset implicit H - will be recalculated */
        new_atom.implicit_h_count = 0;

        int added_idx = molecule_add_atom_full(scaffold, &new_atom);
        if (added_idx < 0) {
            free(idx_map);
            molecule_free(scaffold);
            return NULL;
        }
        idx_map[i] = new_idx++;
    }

    /* Add bonds between scaffold atoms */
    for (int i = 0; i < mol->num_bonds; i++) {
        const bond_t* bond = molecule_get_bond_const(mol, i);
        int new_a1 = idx_map[bond->atom1];
        int new_a2 = idx_map[bond->atom2];

        if (new_a1 >= 0 && new_a2 >= 0) {
            if (molecule_add_bond(scaffold, new_a1, new_a2, bond->type) < 0) {
                free(idx_map);
                molecule_free(scaffold);
                return NULL;
            }
        }
    }

    free(idx_map);

    /* Recalculate implicit hydrogens and perceive aromaticity */
    molecule_calc_implicit_h(scaffold);
    molecule_find_rings(scaffold);
    molecule_perceive_aromaticity(scaffold);

    return scaffold;
}

cchem_status_t molecule_get_murcko_scaffold(const molecule_t* mol,
                                            char* scaffold_buf,
                                            size_t scaffold_buf_size) {
    if (!mol || !scaffold_buf || scaffold_buf_size < 2) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* For molecules without rings, scaffold is empty */
    if (mol->num_rings == 0) {
        scaffold_buf[0] = '\0';
        return CCHEM_OK;
    }

    /* Allocate working arrays */
    bool* in_ring = calloc(mol->num_atoms, sizeof(bool));
    bool* is_linker = calloc(mol->num_atoms, sizeof(bool));
    bool* in_scaffold = calloc(mol->num_atoms, sizeof(bool));

    if (!in_ring || !is_linker || !in_scaffold) {
        free(in_ring);
        free(is_linker);
        free(in_scaffold);
        return CCHEM_ERROR_MEMORY;
    }

    /* Find ring atoms */
    find_ring_atoms(mol, in_ring);

    /* Find linker atoms */
    find_linker_atoms(mol, in_ring, is_linker);

    /* Combine: scaffold = ring atoms + linker atoms */
    for (int i = 0; i < mol->num_atoms; i++) {
        in_scaffold[i] = in_ring[i] || is_linker[i];
    }

    /* Create scaffold molecule */
    molecule_t* scaffold_mol = create_scaffold_molecule(mol, in_scaffold);

    free(in_ring);
    free(is_linker);
    free(in_scaffold);

    if (!scaffold_mol) {
        return CCHEM_ERROR_MEMORY;
    }

    /* Generate canonical SMILES for scaffold */
    char* smiles = molecule_to_canonical_smiles(scaffold_mol, NULL);
    molecule_free(scaffold_mol);

    if (!smiles) {
        scaffold_buf[0] = '\0';
        return CCHEM_OK;  /* Empty scaffold is valid */
    }

    /* Copy to output buffer */
    size_t len = strlen(smiles);
    if (len >= scaffold_buf_size) {
        len = scaffold_buf_size - 1;
    }
    memcpy(scaffold_buf, smiles, len);
    scaffold_buf[len] = '\0';

    free(smiles);
    return CCHEM_OK;
}

cchem_status_t smiles_get_murcko_scaffold(const char* smiles,
                                          char* scaffold_buf,
                                          size_t scaffold_buf_size,
                                          char* error_buf,
                                          size_t error_buf_size) {
    if (!smiles || !scaffold_buf) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Invalid input");
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Parse SMILES */
    molecule_t* mol = smiles_to_molecule(smiles, error_buf, error_buf_size);
    if (!mol) {
        return CCHEM_ERROR_PARSE;
    }

    /* Find rings if not already done */
    if (!mol->rings_computed) {
        molecule_find_rings(mol);
    }

    cchem_status_t status = molecule_get_murcko_scaffold(mol, scaffold_buf, scaffold_buf_size);

    molecule_free(mol);
    return status;
}

/* ========================================================================
 * Split result management
 * ======================================================================== */

split_result_t* split_result_create(int num_rows, int num_splits) {
    split_result_t* result = calloc(1, sizeof(split_result_t));
    if (!result) return NULL;

    result->assignments = calloc(num_rows, sizeof(split_assignment_t));
    result->split_counts = calloc(num_splits, sizeof(int));

    if (!result->assignments || !result->split_counts) {
        split_result_free(result);
        return NULL;
    }

    result->num_assignments = num_rows;
    result->num_splits = num_splits;
    return result;
}

void split_result_free(split_result_t* result) {
    if (result) {
        free(result->assignments);
        free(result->split_counts);
        free(result);
    }
}

/* ========================================================================
 * Random splitting
 * ======================================================================== */

/* Fisher-Yates shuffle */
static void shuffle_indices(int* indices, int n, unsigned int seed) {
    if (seed == 0) {
        seed = (unsigned int)time(NULL);
    }
    srand(seed);

    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
    }
}

split_result_t* splitter_random_split(int num_rows, const split_options_t* options) {
    if (num_rows <= 0 || !options || options->num_splits < 2) {
        return NULL;
    }

    split_result_t* result = split_result_create(num_rows, options->num_splits);
    if (!result) return NULL;

    /* Create and shuffle indices */
    int* indices = malloc(num_rows * sizeof(int));
    if (!indices) {
        split_result_free(result);
        return NULL;
    }

    for (int i = 0; i < num_rows; i++) {
        indices[i] = i;
    }

    shuffle_indices(indices, num_rows, options->seed);

    /* Partition according to ratios */
    int assigned = 0;
    for (int split = 0; split < options->num_splits; split++) {
        int split_size;
        if (split == options->num_splits - 1) {
            /* Last split gets all remaining */
            split_size = num_rows - assigned;
        } else {
            split_size = (int)round(options->ratios[split] * num_rows);
            /* Ensure we don't exceed remaining */
            if (assigned + split_size > num_rows) {
                split_size = num_rows - assigned;
            }
        }

        for (int i = 0; i < split_size; i++) {
            int idx = indices[assigned + i];
            result->assignments[idx].split_idx = split;
            result->assignments[idx].original_idx = idx;
        }

        result->split_counts[split] = split_size;
        assigned += split_size;
    }

    free(indices);
    return result;
}

/* ========================================================================
 * Scaffold-based splitting
 * ======================================================================== */

/* Scaffold computation task argument */
typedef struct {
    int row_idx;
    const char* smiles;
    char scaffold[512];
    bool success;
} scaffold_task_arg_t;

/* Worker function for scaffold computation */
static void* scaffold_worker(void* arg) {
    scaffold_task_arg_t* task = (scaffold_task_arg_t*)arg;
    task->scaffold[0] = '\0';

    if (!task->smiles || task->smiles[0] == '\0') {
        task->success = false;
        return task;
    }

    char error_buf[256];
    cchem_status_t status = smiles_get_murcko_scaffold(
        task->smiles, task->scaffold, sizeof(task->scaffold),
        error_buf, sizeof(error_buf));

    task->success = (status == CCHEM_OK);
    return task;
}

/* Scaffold group for grouping molecules */
typedef struct {
    char* scaffold;      /* Canonical scaffold SMILES */
    int* indices;        /* Molecule indices in this group */
    int count;           /* Number of molecules */
    int capacity;        /* Allocated capacity */
} scaffold_group_t;

/* Add index to scaffold group */
static int scaffold_group_add(scaffold_group_t* group, int idx) {
    if (group->count >= group->capacity) {
        int new_cap = group->capacity == 0 ? 16 : group->capacity * 2;
        int* new_indices = realloc(group->indices, new_cap * sizeof(int));
        if (!new_indices) return -1;
        group->indices = new_indices;
        group->capacity = new_cap;
    }
    group->indices[group->count++] = idx;
    return 0;
}

split_result_t* splitter_scaffold_split(const char** smiles,
                                        int num_molecules,
                                        const split_options_t* options,
                                        char* error_buf,
                                        size_t error_buf_size) {
    if (!smiles || num_molecules <= 0 || !options || options->num_splits < 2) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Invalid input parameters");
        return NULL;
    }

    int num_threads = options->num_threads > 0 ? options->num_threads : parallel_get_num_cores();

    /* Create thread pool for scaffold computation */
    thread_pool_t* pool = thread_pool_create(num_threads);
    if (!pool) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Failed to create thread pool");
        return NULL;
    }

    thread_pool_ensure_capacity(pool, num_molecules);

    /* Allocate task arguments */
    scaffold_task_arg_t* tasks = calloc(num_molecules, sizeof(scaffold_task_arg_t));
    if (!tasks) {
        thread_pool_free(pool);
        if (error_buf) snprintf(error_buf, error_buf_size, "Memory allocation failed");
        return NULL;
    }

    /* Progress bar for scaffold computation */
    progress_t* progress = NULL;
    if (options->verbose) {
        progress_config_t prog_config = PROGRESS_CONFIG_DEFAULT;
        prog_config.prefix = "Scaffolds";
        progress = progress_create(num_molecules, &prog_config);
    }

    /* Submit scaffold computation tasks */
    for (int i = 0; i < num_molecules; i++) {
        tasks[i].row_idx = i;
        tasks[i].smiles = smiles[i];
        thread_pool_submit(pool, scaffold_worker, &tasks[i]);
    }

    /* Wait for completion with progress updates */
    while (thread_pool_num_completed(pool) < num_molecules) {
        if (progress) {
            progress_update(progress, thread_pool_num_completed(pool));
        }
        usleep(50000);
    }

    if (progress) {
        progress_finish(progress);
        progress_free(progress);
    }

    /* Group molecules by scaffold */
    scaffold_group_t* groups = NULL;
    int num_groups = 0;
    int groups_capacity = 0;

    for (int i = 0; i < num_molecules; i++) {
        scaffold_task_arg_t* task = &tasks[i];

        /* Find or create group for this scaffold */
        int group_idx = -1;
        for (int g = 0; g < num_groups; g++) {
            if (strcmp(groups[g].scaffold, task->scaffold) == 0) {
                group_idx = g;
                break;
            }
        }

        if (group_idx < 0) {
            /* Create new group */
            if (num_groups >= groups_capacity) {
                int new_cap = groups_capacity == 0 ? 256 : groups_capacity * 2;
                scaffold_group_t* new_groups = realloc(groups, new_cap * sizeof(scaffold_group_t));
                if (!new_groups) {
                    goto cleanup_error;
                }
                groups = new_groups;
                /* Initialize new entries */
                for (int j = groups_capacity; j < new_cap; j++) {
                    memset(&groups[j], 0, sizeof(scaffold_group_t));
                }
                groups_capacity = new_cap;
            }

            group_idx = num_groups++;
            groups[group_idx].scaffold = strdup(task->scaffold);
            if (!groups[group_idx].scaffold) {
                goto cleanup_error;
            }
        }

        if (scaffold_group_add(&groups[group_idx], i) < 0) {
            goto cleanup_error;
        }
    }

    thread_pool_free(pool);
    pool = NULL;

    if (options->verbose) {
        printf("Found %d unique scaffolds\n", num_groups);
    }

    /* Sort groups by size (largest first) for better distribution */
    for (int i = 0; i < num_groups - 1; i++) {
        for (int j = i + 1; j < num_groups; j++) {
            if (groups[j].count > groups[i].count) {
                scaffold_group_t tmp = groups[i];
                groups[i] = groups[j];
                groups[j] = tmp;
            }
        }
    }

    /* Shuffle groups for randomness */
    unsigned int seed = options->seed ? options->seed : (unsigned int)time(NULL);
    srand(seed);

    for (int i = num_groups - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        scaffold_group_t tmp = groups[i];
        groups[i] = groups[j];
        groups[j] = tmp;
    }

    /* Create result */
    split_result_t* result = split_result_create(num_molecules, options->num_splits);
    if (!result) {
        goto cleanup_error;
    }
    result->num_scaffolds = num_groups;

    /* Assign scaffold groups to splits */
    int* split_sizes = calloc(options->num_splits, sizeof(int));
    int* target_sizes = malloc(options->num_splits * sizeof(int));

    if (!split_sizes || !target_sizes) {
        free(split_sizes);
        free(target_sizes);
        split_result_free(result);
        result = NULL;
        goto cleanup_error;
    }

    /* Calculate target sizes */
    int remaining = num_molecules;
    for (int s = 0; s < options->num_splits; s++) {
        if (s == options->num_splits - 1) {
            target_sizes[s] = remaining;
        } else {
            target_sizes[s] = (int)round(options->ratios[s] * num_molecules);
            remaining -= target_sizes[s];
        }
    }

    if (options->stratified) {
        /* Stratified split: distribute each scaffold group across all splits proportionally */
        for (int g = 0; g < num_groups; g++) {
            /* Shuffle indices within this scaffold group */
            for (int i = groups[g].count - 1; i > 0; i--) {
                int j = rand() % (i + 1);
                int tmp = groups[g].indices[i];
                groups[g].indices[i] = groups[g].indices[j];
                groups[g].indices[j] = tmp;
            }

            /* Distribute molecules from this group across splits according to ratios */
            int assigned = 0;
            for (int s = 0; s < options->num_splits; s++) {
                int group_split_size;
                if (s == options->num_splits - 1) {
                    /* Last split gets remaining molecules from this group */
                    group_split_size = groups[g].count - assigned;
                } else {
                    group_split_size = (int)round(options->ratios[s] * groups[g].count);
                    /* Ensure at least 1 if group is large enough and we haven't assigned all */
                    if (group_split_size == 0 && groups[g].count > options->num_splits &&
                        assigned < groups[g].count) {
                        group_split_size = 1;
                    }
                    /* Don't exceed remaining */
                    if (assigned + group_split_size > groups[g].count) {
                        group_split_size = groups[g].count - assigned;
                    }
                }

                for (int i = 0; i < group_split_size; i++) {
                    int mol_idx = groups[g].indices[assigned + i];
                    result->assignments[mol_idx].split_idx = s;
                    result->assignments[mol_idx].original_idx = mol_idx;
                    split_sizes[s]++;
                }
                assigned += group_split_size;
            }
        }
    } else {
        /* Non-stratified: assign entire scaffold groups to one split (greedy algorithm) */
        for (int g = 0; g < num_groups; g++) {
            /* Find split with largest deficit */
            int best_split = 0;
            int best_deficit = target_sizes[0] - split_sizes[0];

            for (int s = 1; s < options->num_splits; s++) {
                int deficit = target_sizes[s] - split_sizes[s];
                if (deficit > best_deficit) {
                    best_deficit = deficit;
                    best_split = s;
                }
            }

            /* Assign all molecules in this group to the chosen split */
            for (int i = 0; i < groups[g].count; i++) {
                int mol_idx = groups[g].indices[i];
                result->assignments[mol_idx].split_idx = best_split;
                result->assignments[mol_idx].original_idx = mol_idx;
            }

            split_sizes[best_split] += groups[g].count;
        }
    }

    /* Copy final counts */
    for (int s = 0; s < options->num_splits; s++) {
        result->split_counts[s] = split_sizes[s];
    }

    free(split_sizes);
    free(target_sizes);

    /* Cleanup */
    for (int g = 0; g < num_groups; g++) {
        free(groups[g].scaffold);
        free(groups[g].indices);
    }
    free(groups);
    free(tasks);

    return result;

cleanup_error:
    if (pool) thread_pool_free(pool);
    for (int g = 0; g < num_groups; g++) {
        free(groups[g].scaffold);
        free(groups[g].indices);
    }
    free(groups);
    free(tasks);

    if (error_buf) snprintf(error_buf, error_buf_size, "Memory allocation failed");
    return NULL;
}

/* ========================================================================
 * Main CSV splitting function
 * ======================================================================== */

cchem_status_t splitter_split_csv(const char* input_file,
                                  const char** output_files,
                                  int num_outputs,
                                  const char* smiles_col,
                                  const split_options_t* options,
                                  char* error_buf,
                                  size_t error_buf_size) {
    if (!input_file || !output_files || !smiles_col || !options) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Invalid input parameters");
        return CCHEM_ERROR_INVALID_INPUT;
    }

    if (num_outputs != options->num_splits) {
        if (error_buf) snprintf(error_buf, error_buf_size,
                               "Number of output files (%d) must match number of splits (%d)",
                               num_outputs, options->num_splits);
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Read input CSV */
    csv_document_t* doc = csv_document_create();
    if (!doc) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Memory allocation failed");
        return CCHEM_ERROR_MEMORY;
    }

    if (csv_document_read(doc, input_file, true) != CSV_OK) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Failed to read input file: %s", input_file);
        csv_document_free(doc);
        return CCHEM_ERROR_FILE_IO;
    }

    /* Find SMILES column */
    int smiles_col_idx = csv_find_column(&doc->header, smiles_col);
    if (smiles_col_idx < 0) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Column '%s' not found", smiles_col);
        csv_document_free(doc);
        return CCHEM_ERROR_INVALID_INPUT;
    }

    int num_rows = doc->num_rows;

    if (options->verbose) {
        printf("Input file: %s\n", input_file);
        printf("Rows: %d\n", num_rows);
        printf("SMILES column: %s (index %d)\n", smiles_col, smiles_col_idx);
        printf("Split method: %s\n",
               options->method == SPLIT_METHOD_RANDOM ? "random" : "scaffold");
        printf("Split ratios:");
        for (int i = 0; i < options->num_splits; i++) {
            printf(" %.1f%%", options->ratios[i] * 100.0);
        }
        printf("\n\n");
    }

    /* Perform split */
    split_result_t* result = NULL;

    if (options->method == SPLIT_METHOD_RANDOM) {
        result = splitter_random_split(num_rows, options);
    } else {
        /* Extract SMILES strings for scaffold splitting */
        const char** smiles_arr = malloc(num_rows * sizeof(char*));
        if (!smiles_arr) {
            if (error_buf) snprintf(error_buf, error_buf_size, "Memory allocation failed");
            csv_document_free(doc);
            return CCHEM_ERROR_MEMORY;
        }

        for (int i = 0; i < num_rows; i++) {
            csv_row_t* row = &doc->rows[i];
            smiles_arr[i] = (smiles_col_idx < row->num_fields) ?
                            row->fields[smiles_col_idx] : "";
        }

        result = splitter_scaffold_split(smiles_arr, num_rows, options, error_buf, error_buf_size);
        free(smiles_arr);
    }

    if (!result) {
        if (error_buf && error_buf[0] == '\0') {
            snprintf(error_buf, error_buf_size, "Split operation failed");
        }
        csv_document_free(doc);
        return CCHEM_ERROR_MEMORY;
    }

    /* Create output writers */
    csv_writer_t** writers = calloc(num_outputs, sizeof(csv_writer_t*));
    if (!writers) {
        if (error_buf) snprintf(error_buf, error_buf_size, "Memory allocation failed");
        split_result_free(result);
        csv_document_free(doc);
        return CCHEM_ERROR_MEMORY;
    }

    for (int i = 0; i < num_outputs; i++) {
        writers[i] = csv_writer_create(output_files[i]);
        if (!writers[i]) {
            if (error_buf) snprintf(error_buf, error_buf_size,
                                   "Failed to create output file: %s", output_files[i]);
            for (int j = 0; j < i; j++) {
                csv_writer_free(writers[j]);
            }
            free(writers);
            split_result_free(result);
            csv_document_free(doc);
            return CCHEM_ERROR_FILE_IO;
        }

        /* Write header to each output */
        csv_writer_write_fields(writers[i],
                               (const char**)doc->header.fields,
                               doc->header.num_fields);
    }

    /* Write rows to appropriate output files */
    for (int i = 0; i < num_rows; i++) {
        int split_idx = result->assignments[i].split_idx;
        csv_row_t* row = &doc->rows[i];

        csv_writer_write_fields(writers[split_idx],
                               (const char**)row->fields,
                               row->num_fields);
    }

    /* Cleanup writers */
    for (int i = 0; i < num_outputs; i++) {
        csv_writer_free(writers[i]);
    }
    free(writers);

    /* Print summary */
    if (options->verbose) {
        printf("\nSplit summary:\n");
        if (options->method == SPLIT_METHOD_SCAFFOLD) {
            printf("  Unique scaffolds: %d\n", result->num_scaffolds);
        }
    }

    for (int i = 0; i < num_outputs; i++) {
        double pct = num_rows > 0 ? (double)result->split_counts[i] / num_rows * 100.0 : 0.0;
        printf("  %s: %d rows (%.1f%%)\n", output_files[i], result->split_counts[i], pct);
    }

    split_result_free(result);
    csv_document_free(doc);

    return CCHEM_OK;
}
