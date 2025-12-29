/**
 * @file cmd_validate.c
 * @brief Validate command implementation
 */

#include "cchem/utils/commands.h"

void print_validate_usage(const char* prog_name) {
    printf("Usage: %s validate [options]\n\n", prog_name);
    printf("Validate SMILES syntax\n\n");
    printf("Options:\n");
    printf("  -S, --smiles <SMILES>   SMILES string to validate\n");
    printf("  -v, --verbose           Verbose output\n");
    printf("  -h, --help              Print this help message\n");
}

int cmd_validate(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"smiles",  required_argument, 0, 'S'},
        {"verbose", no_argument,       0, 'v'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* smiles = NULL;
    bool verbose = false;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "S:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'S':
                smiles = optarg;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_validate_usage(argv[0]);
                return 0;
            default:
                print_validate_usage(argv[0]);
                return 1;
        }
    }

    if (!smiles) {
        fprintf(stderr, "Error: SMILES string required (-S/--smiles)\n");
        return 1;
    }

    char error_buf[256];
    cchem_status_t status = smiles_validate(smiles, error_buf, sizeof(error_buf));

    if (status == CCHEM_OK) {
        if (verbose) {
            printf("Valid SMILES: %s\n", smiles);
        } else {
            printf("valid\n");
        }
        return 0;
    } else {
        if (verbose) {
            printf("Invalid SMILES: %s\n", smiles);
            printf("Error: %s\n", error_buf);
        } else {
            printf("invalid: %s\n", error_buf);
        }
        return 1;
    }
}
