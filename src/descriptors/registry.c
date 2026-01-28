/**
 * @file registry.c
 * @brief Descriptor registration and lookup system
 */

#include "cchem/compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include "cchem/utils/descriptors.h"
#include "cchem/cchem.h"

/* Global registry instance */
descriptor_registry_t g_descriptor_registry = {0};

/* Global cache generation counter - incremented each time a new molecule is parsed */
static _Thread_local uint64_t g_cache_generation = 0;

/* Get current cache generation */
uint64_t descriptor_cache_generation(void) {
    return g_cache_generation;
}

/* Increment cache generation (called before computing descriptors on new molecule) */
static void descriptor_cache_invalidate(void) {
    g_cache_generation++;
}

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

/* Convert string to lowercase (for pre-computation at registration time) */
static void str_to_lower(char* dst, const char* src, size_t max_len) {
    size_t i = 0;
    while (src[i] && i < max_len - 1) {
        dst[i] = (char)tolower((unsigned char)src[i]);
        i++;
    }
    dst[i] = '\0';
}

/* FNV-1a hash for already-lowercase string (fast path) */
static uint32_t hash_name_lower(const char* name_lower) {
    uint32_t hash = 2166136261u;  /* FNV offset basis */
    while (*name_lower) {
        hash ^= (uint32_t)(unsigned char)*name_lower;
        hash *= 16777619u;  /* FNV prime */
        name_lower++;
    }
    return hash & DESC_HASH_TABLE_MASK;
}

/* FNV-1a hash for case-insensitive string hashing (used at registration time only) */
static uint32_t hash_name(const char* name) {
    uint32_t hash = 2166136261u;  /* FNV offset basis */
    while (*name) {
        hash ^= (uint32_t)tolower((unsigned char)*name);
        hash *= 16777619u;  /* FNV prime */
        name++;
    }
    return hash & DESC_HASH_TABLE_MASK;
}

/* Initialize hash table */
static void init_hash_table(void) {
    for (int i = 0; i < DESC_HASH_TABLE_SIZE; i++) {
        g_descriptor_registry.hash_table[i].index = -1;
        g_descriptor_registry.hash_table[i].next = -1;
    }
}

/* Add descriptor to hash table */
static void hash_table_insert(int desc_index) {
    const char* name = g_descriptor_registry.descriptors[desc_index].name;
    uint32_t bucket = hash_name(name);

    if (g_descriptor_registry.hash_table[bucket].index == -1) {
        /* Empty bucket - direct insert */
        g_descriptor_registry.hash_table[bucket].index = (int16_t)desc_index;
    } else {
        /* Collision - chain at end */
        int16_t curr = (int16_t)bucket;
        while (g_descriptor_registry.hash_table[curr].next != -1) {
            curr = g_descriptor_registry.hash_table[curr].next;
        }
        /* Find empty slot for chaining (use descriptor slots after hash table) */
        /* For simplicity, use linear probing to find empty bucket */
        uint32_t probe = (bucket + 1) & DESC_HASH_TABLE_MASK;
        while (g_descriptor_registry.hash_table[probe].index != -1) {
            probe = (probe + 1) & DESC_HASH_TABLE_MASK;
        }
        g_descriptor_registry.hash_table[curr].next = (int16_t)probe;
        g_descriptor_registry.hash_table[probe].index = (int16_t)desc_index;
    }
}

/* ============================================================================
 * Registry Management
 * ============================================================================ */

void descriptors_init(void) {
    if (g_descriptor_registry.initialized) return;

    memset(&g_descriptor_registry, 0, sizeof(g_descriptor_registry));
    init_hash_table();
    g_descriptor_registry.initialized = true;

    /* Register all built-in descriptors */
    descriptors_register_counts();
    descriptors_register_logp();
    descriptors_register_ratios();
    descriptors_register_electronic();
    descriptors_register_steric();
    descriptors_register_energetic();
    descriptors_register_bitstring();
    descriptors_register_acidbase();
    descriptors_register_solubility();
    descriptors_register_topology();
    descriptors_register_fractional();
    descriptors_register_hash();
    descriptors_register_graph();
    descriptors_register_autocorrelations();
    descriptors_register_logd74();
    descriptors_register_functional();
    descriptors_register_pharmacophore();
    descriptors_register_mqn();
    descriptors_register_estate();
    descriptors_register_vsa();
    descriptors_register_chi();
    descriptors_register_atompairs();
    descriptors_register_bcut();
    descriptors_register_zagreb();
    descriptors_register_infocontent();
    descriptors_register_walkcounts();
    descriptors_register_estate_sums();
    descriptors_register_eta();
    descriptors_register_ringcomplexity();
    descriptors_register_cpsa();
    descriptors_register_bcut_ext();
    descriptors_register_moments();
    descriptors_register_aromatic();
    descriptors_register_atompairs_ext();
    descriptors_register_framework();
    descriptors_register_constitutional();
}

void descriptors_cleanup(void) {
    memset(&g_descriptor_registry, 0, sizeof(g_descriptor_registry));
}

cchem_status_t descriptor_register(const descriptor_def_t* def) {
    if (!def || !def->name[0] || !def->compute) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    if (!g_descriptor_registry.initialized) {
        descriptors_init();
    }

    if (g_descriptor_registry.num_descriptors >= MAX_DESCRIPTORS) {
        return CCHEM_ERROR_MEMORY;
    }

    /* Check for duplicate name */
    if (descriptor_get(def->name) != NULL) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Copy descriptor definition */
    int new_index = g_descriptor_registry.num_descriptors;
    descriptor_def_t* reg = &g_descriptor_registry.descriptors[new_index];
    memcpy(reg, def, sizeof(descriptor_def_t));
    reg->registered = true;

    /* Pre-compute lowercase name for fast lookup */
    str_to_lower(reg->name_lower, reg->name, MAX_DESCRIPTOR_NAME);

    g_descriptor_registry.num_descriptors++;

    /* Add to hash table for fast lookup */
    hash_table_insert(new_index);

    return CCHEM_OK;
}

const descriptor_def_t* descriptor_get(const char* name) {
    if (!name || !name[0]) return NULL;

    /* Pre-compute lowercase query once (on stack, no allocation) */
    char name_lower[MAX_DESCRIPTOR_NAME];
    str_to_lower(name_lower, name, MAX_DESCRIPTOR_NAME);

    /* Use hash table for O(1) lookup - hash the lowercase string directly */
    uint32_t bucket = hash_name_lower(name_lower);
    int16_t curr = (int16_t)bucket;

    while (curr != -1 && g_descriptor_registry.hash_table[curr].index != -1) {
        int idx = g_descriptor_registry.hash_table[curr].index;
        if (idx >= 0 && idx < g_descriptor_registry.num_descriptors) {
            /* Fast strcmp on pre-computed lowercase names (no per-char tolower) */
            if (strcmp(g_descriptor_registry.descriptors[idx].name_lower, name_lower) == 0) {
                return &g_descriptor_registry.descriptors[idx];
            }
        }
        curr = g_descriptor_registry.hash_table[curr].next;
    }
    return NULL;
}

int descriptor_get_by_category(descriptor_category_t category,
                               const descriptor_def_t** out_descriptors,
                               int max_descriptors) {
    if (!out_descriptors || max_descriptors <= 0) return 0;

    int count = 0;
    for (int i = 0; i < g_descriptor_registry.num_descriptors && count < max_descriptors; i++) {
        if (g_descriptor_registry.descriptors[i].category == category) {
            out_descriptors[count++] = &g_descriptor_registry.descriptors[i];
        }
    }
    return count;
}

int descriptor_get_all(const descriptor_def_t** out_descriptors, int max_descriptors) {
    if (!out_descriptors || max_descriptors <= 0) return 0;

    int count = 0;
    for (int i = 0; i < g_descriptor_registry.num_descriptors && count < max_descriptors; i++) {
        out_descriptors[count++] = &g_descriptor_registry.descriptors[i];
    }
    return count;
}

int descriptor_count(void) {
    return g_descriptor_registry.num_descriptors;
}

void descriptor_list_all(void) {
    printf("Available descriptors (%d):\n", g_descriptor_registry.num_descriptors);
    printf("%-25s %-15s %s\n", "Name", "Category", "Description");
    printf("%-25s %-15s %s\n", "----", "--------", "-----------");

    for (int i = 0; i < g_descriptor_registry.num_descriptors; i++) {
        const descriptor_def_t* d = &g_descriptor_registry.descriptors[i];
        printf("%-25s %-15s %s\n",
               d->name,
               descriptor_category_name(d->category),
               d->description);
    }
}

/* ============================================================================
 * Descriptor Computation
 * ============================================================================ */

cchem_status_t descriptor_compute(const molecule_t* mol,
                                  const char* name,
                                  descriptor_value_t* value) {
    if (!mol || !name || !value) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    const descriptor_def_t* def = descriptor_get(name);
    if (!def) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    return def->compute(mol, value);
}

cchem_status_t descriptor_compute_from_smiles(const char* smiles,
                                              const char* name,
                                              descriptor_value_t* value,
                                              char* error_buf,
                                              size_t error_buf_size) {
    if (!smiles || !name || !value) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "Invalid input parameters");
        }
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Canonicalize SMILES first */
    char* canonical = smiles_canonicalize(smiles, NULL, error_buf, error_buf_size);
    if (!canonical) {
        return CCHEM_ERROR_INVALID_SMILES;
    }

    /* Parse canonical SMILES into molecule */
    molecule_t* mol = smiles_to_molecule(canonical, error_buf, error_buf_size);
    free(canonical);

    if (!mol) {
        return CCHEM_ERROR_PARSE;
    }

    /* Invalidate descriptor caches for new molecule */
    descriptor_cache_invalidate();

    /* Compute descriptor */
    cchem_status_t status = descriptor_compute(mol, name, value);
    molecule_free(mol);

    return status;
}

cchem_status_t descriptors_compute_batch(const molecule_t* mol,
                                         const char** names,
                                         int num_descriptors,
                                         descriptor_value_t* values) {
    if (!mol || !names || !values || num_descriptors <= 0) {
        return CCHEM_ERROR_INVALID_INPUT;
    }

    /* Group descriptors by category for optimized batch computation */
    const descriptor_def_t* descs[MAX_DESCRIPTORS];

    for (int i = 0; i < num_descriptors; i++) {
        descs[i] = descriptor_get(names[i]);
        if (!descs[i]) {
            return CCHEM_ERROR_INVALID_INPUT;
        }
    }

    /* Check if we have a batch compute function for the category */
    /* For now, fall back to individual computation */
    /* TODO: Implement optimized batch compute for each category */

    for (int i = 0; i < num_descriptors; i++) {
        cchem_status_t status = descs[i]->compute(mol, &values[i]);
        if (status != CCHEM_OK) {
            return status;
        }
    }

    return CCHEM_OK;
}

int descriptors_compute_all(const molecule_t* mol,
                            const char** out_names,
                            descriptor_value_t* out_values,
                            int max_descriptors) {
    if (!mol || !out_values || max_descriptors <= 0) {
        return -1;
    }

    int count = 0;

    /* Use optimized batch computation for categories that support it.
     * Categories are registered in order: counts, logp, ratios, electronic, steric, energetic, bitstring.
     * We iterate through and use batch functions where available. */

    for (int i = 0; i < g_descriptor_registry.num_descriptors && count < max_descriptors; i++) {
        const descriptor_def_t* def = &g_descriptor_registry.descriptors[i];
        descriptor_category_t cat = def->category;

        /* Check if this is the first descriptor of a category with batch support */
        bool is_first_in_cat = (i == 0) || (g_descriptor_registry.descriptors[i-1].category != cat);

        if (is_first_in_cat) {
            /* Try batch computation for supported categories */
            int batch_count = 0;
            descriptor_value_t batch_values[128];  /* Enough for largest category (counts: 83) */

            switch (cat) {
                case DESC_CATEGORY_COUNTS:
                    batch_count = descriptors_compute_counts_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_ELECTRONIC:
                    batch_count = descriptors_compute_electronic_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_RATIOS:
                    batch_count = descriptors_compute_ratios_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_STERIC:
                    batch_count = descriptors_compute_steric_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_ENERGETIC:
                    batch_count = descriptors_compute_energetic_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_BITSTRING:
                    batch_count = descriptors_compute_bitstring_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_FINGERPRINT:
                    batch_count = descriptors_compute_hash_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_PROPERTIES:
                    batch_count = descriptors_compute_fractional_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_CUSTOM:
                    batch_count = descriptors_compute_graph_all(mol, batch_values);
                    break;
                case DESC_CATEGORY_CONSTITUTIONAL:
                    batch_count = descriptors_compute_constitutional_all(mol, batch_values);
                    break;
                default:
                    batch_count = 0;
                    break;
            }

            if (batch_count > 0) {
                /* Copy batch results */
                int j = i;
                int b = 0;
                while (j < g_descriptor_registry.num_descriptors &&
                       g_descriptor_registry.descriptors[j].category == cat &&
                       count < max_descriptors &&
                       b < batch_count) {
                    out_values[count] = batch_values[b];
                    if (out_names) {
                        out_names[count] = g_descriptor_registry.descriptors[j].name;
                    }
                    count++;
                    j++;
                    b++;
                }
                /* Skip to end of this category */
                i = j - 1;
                continue;
            }
        }

        /* Fallback to individual computation */
        cchem_status_t status = def->compute(mol, &out_values[count]);
        if (status == CCHEM_OK) {
            if (out_names) {
                out_names[count] = def->name;
            }
            count++;
        }
    }

    return count;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

char* descriptor_format_value(const descriptor_def_t* def,
                              const descriptor_value_t* value,
                              char* buf,
                              size_t buf_size) {
    if (!def || !value || !buf || buf_size == 0) {
        return NULL;
    }

    if (def->value_type == DESC_VALUE_INT) {
        snprintf(buf, buf_size, "%lld", (long long)value->i);
    } else {
        snprintf(buf, buf_size, "%.6g", value->d);
    }

    return buf;
}

int descriptor_parse_names(const char* input,
                           char** out_names,
                           int max_names) {
    if (!input || !out_names || max_names <= 0) {
        return -1;
    }

    /* Handle "all" keyword (case-insensitive) */
    char input_lower[8];
    str_to_lower(input_lower, input, sizeof(input_lower));
    if (strcmp(input_lower, "all") == 0) {
        int count = 0;
        for (int i = 0; i < g_descriptor_registry.num_descriptors && count < max_names; i++) {
            out_names[count] = strdup(g_descriptor_registry.descriptors[i].name);
            if (!out_names[count]) {
                descriptor_free_names(out_names, count);
                return -1;
            }
            count++;
        }
        return count;
    }

    /* Parse comma-separated names */
    char* copy = strdup(input);
    if (!copy) return -1;

    int count = 0;
    char* token = strtok(copy, ",");

    while (token && count < max_names) {
        /* Skip leading whitespace */
        while (*token && isspace((unsigned char)*token)) token++;

        /* Skip trailing whitespace */
        char* end = token + strlen(token) - 1;
        while (end > token && isspace((unsigned char)*end)) *end-- = '\0';

        if (*token) {
            /* Verify descriptor exists */
            if (!descriptor_get(token)) {
                free(copy);
                descriptor_free_names(out_names, count);
                return -1;
            }

            out_names[count] = strdup(token);
            if (!out_names[count]) {
                free(copy);
                descriptor_free_names(out_names, count);
                return -1;
            }
            count++;
        }

        token = strtok(NULL, ",");
    }

    free(copy);
    return count;
}

void descriptor_free_names(char** names, int num_names) {
    if (!names) return;
    for (int i = 0; i < num_names; i++) {
        if (names[i]) free(names[i]);
    }
}

const char* descriptor_category_name(descriptor_category_t category) {
    switch (category) {
        case DESC_CATEGORY_COUNTS:        return "Counts";
        case DESC_CATEGORY_RATIOS:        return "Ratios";
        case DESC_CATEGORY_PROPERTIES:    return "Properties";
        case DESC_CATEGORY_ELECTRONIC:    return "Electronic";
        case DESC_CATEGORY_STERIC:        return "Steric";
        case DESC_CATEGORY_ENERGETIC:     return "Energetic";
        case DESC_CATEGORY_BITSTRING:     return "Bitstring";
        case DESC_CATEGORY_FINGERPRINT:   return "Fingerprint";
        case DESC_CATEGORY_CUSTOM:        return "Custom";
        case DESC_CATEGORY_CONSTITUTIONAL: return "Constitutional";
        default:                          return "Unknown";
    }
}
