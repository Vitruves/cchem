/**
 * @file canon.c
 * @brief Main canonicalization algorithm
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cchem/canonicalizer/canon.h"
#include "cchem/canonicalizer/invariant.h"
#include "cchem/canonicalizer/ring_finder.h"
#include "cchem/canonicalizer/smiles_writer.h"

/* Default options */
const canon_options_t CANON_OPTIONS_DEFAULT = {
    .preserve_stereo = true,
    .kekulize = false,
    .remove_explicit_h = false,
    .remove_atom_maps = true,
    .canonical_stereo = true,
    .use_isomeric = true,
    .use_isotopes = true,
    .use_charges = true
};

const canon_options_t CANON_OPTIONS_STRICT = {
    .preserve_stereo = true,
    .kekulize = false,
    .remove_explicit_h = false,
    .remove_atom_maps = false,
    .canonical_stereo = true,
    .use_isomeric = true,
    .use_isotopes = true,
    .use_charges = true
};

const canon_options_t CANON_OPTIONS_MINIMAL = {
    .preserve_stereo = false,
    .kekulize = false,
    .remove_explicit_h = true,
    .remove_atom_maps = true,
    .canonical_stereo = false,
    .use_isomeric = false,
    .use_isotopes = false,
    .use_charges = true
};

canon_context_t* canon_context_create(molecule_t* mol, const canon_options_t* options) {
    if (!mol) return NULL;

    canon_context_t* ctx = (canon_context_t*)calloc(1, sizeof(canon_context_t));
    if (!ctx) return NULL;

    ctx->mol = mol;
    ctx->options = options ? *options : CANON_OPTIONS_DEFAULT;
    ctx->num_atoms = mol->num_atoms;
    ctx->computed = false;
    ctx->iterations = 0;
    ctx->error_msg[0] = '\0';

    /* Allocate working arrays - use calloc for deterministic initialization */
    ctx->ranks = (int*)calloc(mol->num_atoms, sizeof(int));
    ctx->invariants = (uint64_t*)calloc(mol->num_atoms, sizeof(uint64_t));
    ctx->order = (int*)calloc(mol->num_atoms, sizeof(int));

    if (!ctx->ranks || !ctx->invariants || !ctx->order) {
        canon_context_free(ctx);
        return NULL;
    }

    /* Initialize order to identity */
    for (int i = 0; i < mol->num_atoms; i++) {
        ctx->order[i] = i;
        ctx->ranks[i] = 0;
    }

    /* Create stereo info if needed */
    if (ctx->options.preserve_stereo) {
        ctx->stereo = stereo_info_create();
    }

    return ctx;
}

void canon_context_free(canon_context_t* ctx) {
    if (!ctx) return;

    if (ctx->ranks) free(ctx->ranks);
    if (ctx->invariants) free(ctx->invariants);
    if (ctx->order) free(ctx->order);
    if (ctx->stereo) stereo_info_free(ctx->stereo);

    free(ctx);
}

cchem_status_t canon_calc_invariants(canon_context_t* ctx) {
    if (!ctx || !ctx->mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Ensure rings are computed */
    if (!ctx->mol->rings_computed) {
        cchem_status_t status = molecule_find_rings(ctx->mol);
        if (status != CCHEM_OK) {
            snprintf(ctx->error_msg, sizeof(ctx->error_msg), "Failed to find rings");
            return status;
        }
    }

    /* Calculate initial invariants */
    cchem_status_t status = invariant_calc_initial(ctx->mol);
    if (status != CCHEM_OK) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg), "Failed to calculate invariants");
        return status;
    }

    /* Copy to context */
    for (int i = 0; i < ctx->num_atoms; i++) {
        ctx->invariants[i] = ctx->mol->atoms[i].invariant;
    }

    return CCHEM_OK;
}

cchem_status_t canon_refine_invariants(canon_context_t* ctx) {
    if (!ctx || !ctx->mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Refine using Morgan algorithm */
    int max_iter = ctx->num_atoms;
    cchem_status_t status = invariant_refine(ctx->mol, max_iter);
    if (status != CCHEM_OK) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg), "Failed to refine invariants");
        return status;
    }

    /* Update context */
    for (int i = 0; i < ctx->num_atoms; i++) {
        ctx->invariants[i] = ctx->mol->atoms[i].invariant;
    }

    return CCHEM_OK;
}

cchem_status_t canon_break_ties(canon_context_t* ctx) {
    if (!ctx || !ctx->mol) return CCHEM_ERROR_INVALID_INPUT;

    cchem_status_t status = invariant_break_ties(ctx->mol);
    if (status != CCHEM_OK) {
        snprintf(ctx->error_msg, sizeof(ctx->error_msg), "Failed to break ties");
        return status;
    }

    /* Update context ranks */
    for (int i = 0; i < ctx->num_atoms; i++) {
        ctx->ranks[i] = ctx->mol->atoms[i].canon_rank;
    }

    return CCHEM_OK;
}

/* Comparison for sorting by rank */
typedef struct {
    int index;
    int rank;
} rank_pair_t;

static int rank_pair_compare(const void* a, const void* b) {
    const rank_pair_t* pa = (const rank_pair_t*)a;
    const rank_pair_t* pb = (const rank_pair_t*)b;
    /* Primary sort by canonical rank */
    if (pa->rank != pb->rank) {
        return pa->rank - pb->rank;
    }
    /* Secondary sort by atom index for deterministic ordering when ranks are equal */
    return pa->index - pb->index;
}

cchem_status_t canon_generate_order(canon_context_t* ctx) {
    if (!ctx || !ctx->mol) return CCHEM_ERROR_INVALID_INPUT;

    /* Create pairs of (index, rank) and sort */
    rank_pair_t* pairs = (rank_pair_t*)calloc(ctx->num_atoms, sizeof(rank_pair_t));
    if (!pairs) return CCHEM_ERROR_MEMORY;

    for (int i = 0; i < ctx->num_atoms; i++) {
        pairs[i].index = i;
        pairs[i].rank = ctx->mol->atoms[i].canon_rank;
    }

    qsort(pairs, ctx->num_atoms, sizeof(rank_pair_t), rank_pair_compare);

    /* Extract order */
    for (int i = 0; i < ctx->num_atoms; i++) {
        ctx->order[i] = pairs[i].index;
    }

    free(pairs);

    /* Store in molecule */
    if (ctx->mol->canon_order) free(ctx->mol->canon_order);
    ctx->mol->canon_order = (int*)calloc(ctx->num_atoms, sizeof(int));
    if (ctx->mol->canon_order) {
        memcpy(ctx->mol->canon_order, ctx->order, ctx->num_atoms * sizeof(int));
    }

    return CCHEM_OK;
}

cchem_status_t canon_stereo(canon_context_t* ctx) {
    if (!ctx || !ctx->mol) return CCHEM_ERROR_INVALID_INPUT;

    if (!ctx->options.preserve_stereo || !ctx->stereo) {
        return CCHEM_OK;
    }

    /* Detect stereocenters */
    cchem_status_t status = stereo_detect_centers(ctx->mol, ctx->stereo);
    if (status != CCHEM_OK) return status;

    /* Canonicalize stereochemistry */
    status = stereo_canonicalize(ctx->mol, ctx->stereo, ctx->order);
    if (status != CCHEM_OK) return status;

    return CCHEM_OK;
}

cchem_status_t canon_canonicalize(canon_context_t* ctx) {
    if (!ctx || !ctx->mol) return CCHEM_ERROR_INVALID_INPUT;

    cchem_status_t status;

    /* Step 1: Calculate initial invariants */
    status = canon_calc_invariants(ctx);
    if (status != CCHEM_OK) return status;

    /* Step 2: Refine invariants */
    status = canon_refine_invariants(ctx);
    if (status != CCHEM_OK) return status;

    /* Step 3: Break ties */
    status = canon_break_ties(ctx);
    if (status != CCHEM_OK) return status;

    /* Step 4: Generate canonical order */
    status = canon_generate_order(ctx);
    if (status != CCHEM_OK) return status;

    /* Step 5: Canonicalize stereochemistry */
    status = canon_stereo(ctx);
    if (status != CCHEM_OK) return status;

    ctx->computed = true;
    ctx->mol->is_canonical = true;

    return CCHEM_OK;
}

const int* canon_get_order(const canon_context_t* ctx) {
    if (!ctx || !ctx->computed) return NULL;
    return ctx->order;
}

int canon_get_rank(const canon_context_t* ctx, int atom_idx) {
    if (!ctx || !ctx->computed) return -1;
    if (atom_idx < 0 || atom_idx >= ctx->num_atoms) return -1;
    return ctx->ranks[atom_idx];
}

cchem_status_t molecule_canonicalize(molecule_t* mol, const canon_options_t* options) {
    if (!mol) return CCHEM_ERROR_INVALID_INPUT;

    canon_context_t* ctx = canon_context_create(mol, options);
    if (!ctx) return CCHEM_ERROR_MEMORY;

    cchem_status_t status = canon_canonicalize(ctx);

    canon_context_free(ctx);
    return status;
}

char* molecule_to_canonical_smiles(const molecule_t* mol, const canon_options_t* options) {
    if (!mol) return NULL;

    /* Clone molecule to avoid modifying original */
    molecule_t* mol_copy = molecule_clone(mol);
    if (!mol_copy) return NULL;

    /* Canonicalize */
    cchem_status_t status = molecule_canonicalize(mol_copy, options);
    if (status != CCHEM_OK) {
        molecule_free(mol_copy);
        return NULL;
    }

    /* Generate SMILES */
    smiles_output_options_t output_opts = SMILES_OUTPUT_CANONICAL;
    if (options) {
        output_opts.show_stereo = options->preserve_stereo;
        output_opts.show_isotopes = options->use_isotopes;
    }

    char* smiles = molecule_to_smiles(mol_copy, &output_opts);

    molecule_free(mol_copy);
    return smiles;
}

char* smiles_canonicalize(const char* smiles, const canon_options_t* options,
                          char* error_buf, size_t error_buf_size) {
    if (!smiles) {
        if (error_buf && error_buf_size > 0) {
            snprintf(error_buf, error_buf_size, "NULL SMILES input");
        }
        return NULL;
    }

    /* Parse SMILES */
    molecule_t* mol = smiles_to_molecule(smiles, error_buf, error_buf_size);
    if (!mol) {
        return NULL;
    }

    /* Generate canonical SMILES */
    char* canonical = molecule_to_canonical_smiles(mol, options);
    if (!canonical && error_buf && error_buf_size > 0) {
        snprintf(error_buf, error_buf_size, "Failed to generate canonical SMILES");
    }

    molecule_free(mol);
    return canonical;
}

bool smiles_are_equivalent(const char* smiles1, const char* smiles2) {
    if (!smiles1 || !smiles2) return false;

    char error_buf[256];

    char* canon1 = smiles_canonicalize(smiles1, NULL, error_buf, sizeof(error_buf));
    if (!canon1) return false;

    char* canon2 = smiles_canonicalize(smiles2, NULL, error_buf, sizeof(error_buf));
    if (!canon2) {
        free(canon1);
        return false;
    }

    bool equal = (strcmp(canon1, canon2) == 0);

    free(canon1);
    free(canon2);

    return equal;
}
