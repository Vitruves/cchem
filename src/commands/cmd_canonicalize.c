/**
 * @file cmd_canonicalize.c
 * @brief Canonicalize command implementation
 */

#include "cchem/utils/commands.h"

void print_canonicalize_usage(const char* prog_name) {
    printf("Usage: %s canonicalize [options]\n\n", prog_name);
    printf("Canonicalize SMILES strings (aromatized by default)\n\n");
    printf("Options:\n");
    printf("  -S, --smiles <SMILES>   Single SMILES string to canonicalize\n");
    printf("  -f, --file <file.csv>   Input CSV file\n");
    printf("  -s, --col <name>        Input SMILES column name (required with -f)\n");
    printf("  -c, --outcol <name>     Output column name for canonical SMILES (default: canonical_smiles)\n");
    printf("  -o, --output <file.csv> Output CSV file (required with -f)\n");
    printf("  -n, --ncpu <N>          Number of CPU cores for parallel processing (default: auto)\n");
    printf("  --sanitize <opts>       Apply sanitization before canonicalization\n");
    printf("                          Values: \"complete\" or comma-separated list of:\n");
    printf("                            unsalt          - Remove salts, keep largest fragment\n");
    printf("                            aromatize       - Perceive and apply aromaticity (default)\n");
    printf("                            kekulize        - Convert aromatic to Kekule form\n");
    printf("                            neutralize      - Neutralize charges\n");
    printf("                            normalize       - Normalize functional groups (nitro, etc.)\n");
    printf("                            remove-stereo   - Remove stereochemistry\n");
    printf("                            remove-isotopes - Remove isotope labels\n");
    printf("                            remove-h        - Remove explicit hydrogens\n");
    printf("                            validate        - Validate structure\n");
    printf("  --list-tautomers        List all tautomeric forms instead of canonicalizing\n");
    printf("  -v, --verbose           Verbose output\n");
    printf("  -h, --help              Print this help message\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s canonicalize -S \"c1ccccc1\"\n", prog_name);
    printf("  %s canonicalize -S \"[Na+].CC(=O)[O-]\" --sanitize unsalt\n", prog_name);
    printf("  %s canonicalize -S \"C1=CC=CC=C1\" --sanitize kekulize\n", prog_name);
    printf("  %s canonicalize -S \"CC(=O)NC\" --list-tautomers\n", prog_name);
    printf("  %s canonicalize -f input.csv -s smiles -c canonical -o output.csv -n 4\n", prog_name);
}

int cmd_canonicalize(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"smiles",        required_argument, 0, 'S'},
        {"file",          required_argument, 0, 'f'},
        {"col",           required_argument, 0, 's'},
        {"outcol",        required_argument, 0, 'c'},
        {"output",        required_argument, 0, 'o'},
        {"ncpu",          required_argument, 0, 'n'},
        {"sanitize",      required_argument, 0, 1000},
        {"list-tautomers", no_argument,      0, 1001},
        {"verbose",       no_argument,       0, 'v'},
        {"help",          no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* smiles = NULL;
    char* input_file = NULL;
    char* input_col = NULL;
    char* output_col = "canonical_smiles";
    char* output_file = NULL;
    int ncpu = 0;  /* 0 = auto */
    bool verbose = false;
    char* sanitize_str = NULL;
    bool list_tautomers = false;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "S:f:s:c:o:n:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'S':
                smiles = optarg;
                break;
            case 'f':
                input_file = optarg;
                break;
            case 's':
                input_col = optarg;
                break;
            case 'c':
                output_col = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'n':
                ncpu = atoi(optarg);
                break;
            case 1000:
                sanitize_str = optarg;
                break;
            case 1001:
                list_tautomers = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_canonicalize_usage(argv[0]);
                return 0;
            default:
                print_canonicalize_usage(argv[0]);
                return 1;
        }
    }

    /* Single SMILES mode */
    if (smiles) {
        char error_buf[256];

        /* Parse sanitization flags if provided */
        sanitize_flags_t sanitize_flags = SANITIZE_NONE;
        if (sanitize_str) {
            if (sanitize_parse_flags(sanitize_str, &sanitize_flags) != CCHEM_OK) {
                fprintf(stderr, "Error: Invalid sanitize option: %s\n", sanitize_str);
                fprintf(stderr, "Valid options: complete, unsalt, aromatize, kekulize, neutralize,\n");
                fprintf(stderr, "               normalize, remove-stereo, remove-isotopes, remove-h, validate\n");
                return 1;
            }
        }

        /* Handle list-tautomers mode */
        if (list_tautomers) {
            molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
            if (!mol) {
                fprintf(stderr, "Error: %s\n", error_buf);
                return 1;
            }

            /* Apply sanitization if requested */
            if (sanitize_flags != SANITIZE_NONE) {
                sanitize_options_t san_opts = SANITIZE_OPTIONS_DEFAULT;
                san_opts.flags = sanitize_flags;
                if (molecule_sanitize(mol, &san_opts, error_buf, sizeof(error_buf)) != CCHEM_OK) {
                    fprintf(stderr, "Error: %s\n", error_buf);
                    molecule_free(mol);
                    return 1;
                }
            }

            tautomer_result_t result;
            tautomer_options_t taut_opts = TAUTOMER_OPTIONS_DEFAULT;

            if (tautomer_enumerate(mol, &taut_opts, &result) != CCHEM_OK) {
                fprintf(stderr, "Error: Tautomer enumeration failed\n");
                molecule_free(mol);
                return 1;
            }

            if (verbose) {
                printf("Found %d tautomer(s):\n", result.num_tautomers);
            }

            for (int i = 0; i < result.num_tautomers; i++) {
                if (verbose) {
                    printf("%d. %s%s\n", i + 1, result.smiles[i],
                           i == result.canonical_idx ? " (canonical)" : "");
                } else {
                    printf("%s\n", result.smiles[i]);
                }
            }

            tautomer_result_free(&result);
            molecule_free(mol);
            return 0;
        }

        /* Standard canonicalization with optional sanitization */
        if (sanitize_flags != SANITIZE_NONE) {
            /* Parse, sanitize, then canonicalize */
            molecule_t* mol = smiles_to_molecule(smiles, error_buf, sizeof(error_buf));
            if (!mol) {
                fprintf(stderr, "Error: %s\n", error_buf);
                return 1;
            }

            /* Apply sanitization */
            sanitize_options_t san_opts = SANITIZE_OPTIONS_DEFAULT;
            san_opts.flags = sanitize_flags;

            if (molecule_sanitize(mol, &san_opts, error_buf, sizeof(error_buf)) != CCHEM_OK) {
                fprintf(stderr, "Error: %s\n", error_buf);
                molecule_free(mol);
                return 1;
            }

            /* Generate canonical SMILES */
            char* canonical = molecule_to_canonical_smiles(mol, NULL);
            molecule_free(mol);

            if (canonical) {
                printf("%s\n", canonical);
                free(canonical);
                return 0;
            } else {
                fprintf(stderr, "Error: Failed to generate canonical SMILES\n");
                return 1;
            }
        }

        /* Standard canonicalization (no sanitization) */
        char* canonical = smiles_canonicalize(smiles, NULL, error_buf, sizeof(error_buf));

        if (canonical) {
            printf("%s\n", canonical);
            free(canonical);
            return 0;
        } else {
            fprintf(stderr, "Error: %s\n", error_buf);
            return 1;
        }
    }

    /* CSV batch mode */
    if (input_file) {
        if (!input_col) {
            fprintf(stderr, "Error: Input column name required (-s/--col)\n");
            return 1;
        }
        if (!output_file) {
            fprintf(stderr, "Error: Output file required (-o/--output)\n");
            return 1;
        }

        if (verbose) {
            printf("Processing: %s\n", input_file);
            printf("Input column: %s\n", input_col);
            printf("Output column: %s\n", output_col);
            printf("Output file: %s\n", output_file);
            printf("CPU cores: %s\n", ncpu > 0 ? "specified" : "auto");
            printf("\n");
        }

        csv_batch_context_t* ctx = csv_batch_context_create(
            input_file, output_file, input_col, output_col, ncpu);

        if (!ctx) {
            fprintf(stderr, "Error: Failed to create batch context\n");
            return 1;
        }

        cchem_status_t status = csv_batch_canonicalize(ctx);

        if (status == CCHEM_OK) {
            size_t total, success, errors;
            csv_batch_get_stats(ctx, &total, &success, &errors);

            if (verbose) {
                printf("\n");
            }
            printf("Processed: %zu/%zu (%.1f%% success)\n",
                   success, total, total > 0 ? (double)success / total * 100.0 : 0.0);

            if (errors > 0) {
                printf("Errors: %zu\n", errors);
            }

            printf("Output written to: %s\n", output_file);

            csv_batch_context_free(ctx);
            return 0;
        } else {
            fprintf(stderr, "Error: Batch processing failed\n");
            csv_batch_context_free(ctx);
            return 1;
        }
    }

    /* No input specified */
    fprintf(stderr, "Error: No input specified. Use -S for single SMILES or -f for CSV file.\n");
    print_canonicalize_usage(argv[0]);
    return 1;
}
