/**
 * @file main.c
 * @brief cchem CLI - SMILES canonicalizer command line interface
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cchem/cchem.h"
#include "cchem/utils/csv.h"
#include "cchem/utils/commands.h"

static void print_version(void) {
    printf("cchem version %s\n", VERSION);
    printf("SMILES canonicalization tool\n");
}

static void print_usage(const char* prog_name) {
    printf("Usage: %s <command> [options]\n\n", prog_name);
    printf("Commands:\n");
    printf("  canonicalize   Canonicalize SMILES strings\n");
    printf("  compute        Compute molecular descriptors\n");
    printf("  depict         Generate 2D/3D molecular structure images\n");
    printf("  split          Split dataset into train/test/validation sets\n");
    printf("  validate       Validate SMILES syntax\n");
    printf("  version        Print version information\n");
    printf("  help           Print this help message\n");
    printf("\n");
    printf("For command-specific help:\n");
    printf("  %s <command> --help\n", prog_name);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char* command = argv[1];

    /* Shift arguments for subcommand parsing */
    argc--;
    argv++;

    if (strcmp(command, "canonicalize") == 0) {
        return cmd_canonicalize(argc, argv);
    } else if (strcmp(command, "compute") == 0) {
        return cmd_compute(argc, argv);
    } else if (strcmp(command, "depict") == 0) {
        return cmd_depict(argc, argv);
    } else if (strcmp(command, "split") == 0) {
        return cmd_split(argc, argv);
    } else if (strcmp(command, "validate") == 0) {
        return cmd_validate(argc, argv);
    } else if (strcmp(command, "version") == 0 || strcmp(command, "--version") == 0) {
        print_version();
        return 0;
    } else if (strcmp(command, "help") == 0 || strcmp(command, "--help") == 0 || strcmp(command, "-h") == 0) {
        print_usage(argv[0]);
        return 0;
    } else {
        fprintf(stderr, "Unknown command: %s\n", command);
        print_usage(argv[0]);
        return 1;
    }
}

/* ============================================================================
 * Library API Functions
 * ============================================================================ */

const char* cchem_version(void) {
    return CCHEM_VERSION_STRING;
}

cchem_status_t cchem_init(void) {
    return CCHEM_OK;
}

void cchem_cleanup(void) {
}

char* cchem_canonicalize(const char* smiles, char* error_buf, size_t error_buf_size) {
    return smiles_canonicalize(smiles, NULL, error_buf, error_buf_size);
}

bool cchem_smiles_equal(const char* smiles1, const char* smiles2) {
    return smiles_are_equivalent(smiles1, smiles2);
}

cchem_status_t cchem_validate_smiles(const char* smiles, char* error_buf, size_t error_buf_size) {
    return smiles_validate(smiles, error_buf, error_buf_size);
}

cchem_status_t cchem_process_csv(const char* input_file, const char* output_file,
                                 const char* smiles_column, const char* output_column,
                                 int num_threads) {
    csv_batch_context_t* ctx = csv_batch_context_create(
        input_file, output_file, smiles_column, output_column, num_threads);

    if (!ctx) return CCHEM_ERROR_MEMORY;

    cchem_status_t status = csv_batch_canonicalize(ctx);
    csv_batch_context_free(ctx);

    return status;
}
