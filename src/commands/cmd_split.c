/**
 * @file cmd_split.c
 * @brief Split command implementation
 */

#include "cchem/utils/commands.h"

void print_split_usage(const char* prog_name) {
    printf("Usage: %s split [options]\n\n", prog_name);
    printf("Split dataset into train/test/validation sets\n\n");
    printf("Options:\n");
    printf("  -f, --file <file.csv>         Input CSV file (required)\n");
    printf("  -o, --output <f1,f2,...>      Comma-separated output file names (required)\n");
    printf("  -s, --col <name>              SMILES column name (required)\n");
    printf("  --split-ratios <ratios>       Split ratios, e.g., \"80,10,10\" or \"0.8,0.1,0.1\"\n");
    printf("                                (default: 80,20 for 2 outputs, 80,10,10 for 3)\n");
    printf("  --splitting-method <method>   Splitting method: 'random' or 'scaffold'\n");
    printf("                                (default: random)\n");
    printf("  --stratified                  For scaffold split: distribute each scaffold\n");
    printf("                                proportionally across all outputs (default: off)\n");
    printf("  -n, --ncpu <N>                Number of CPU cores (default: all)\n");
    printf("  --seed <N>                    Random seed for reproducibility\n");
    printf("  -v, --verbose                 Verbose output\n");
    printf("  -h, --help                    Print this help message\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,val.csv,test.csv --split-ratios 80,10,10\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv --splitting-method scaffold\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv --splitting-method scaffold --stratified\n", prog_name);
    printf("  %s split -f data.csv -s smiles -o train.csv,test.csv --split-ratios 8,2 --seed 42\n", prog_name);
}

int cmd_split(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"file",             required_argument, 0, 'f'},
        {"output",           required_argument, 0, 'o'},
        {"col",              required_argument, 0, 's'},
        {"split-ratios",     required_argument, 0, 1000},
        {"splitting-method", required_argument, 0, 1001},
        {"stratified",       no_argument,       0, 1003},
        {"ncpu",             required_argument, 0, 'n'},
        {"seed",             required_argument, 0, 1002},
        {"verbose",          no_argument,       0, 'v'},
        {"help",             no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* input_file = NULL;
    char* output_str = NULL;
    char* smiles_col = NULL;
    char* ratio_str = NULL;
    char* method_str = "random";
    int ncpu = 0;
    unsigned int seed = 0;
    bool verbose = false;
    bool stratified = false;
    (void)argc;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "f:o:s:n:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'f':
                input_file = optarg;
                break;
            case 'o':
                output_str = optarg;
                break;
            case 's':
                smiles_col = optarg;
                break;
            case 1000:  /* --split-ratios */
                ratio_str = optarg;
                break;
            case 1001:  /* --splitting-method */
                method_str = optarg;
                break;
            case 'n':
                ncpu = atoi(optarg);
                if (ncpu < 0) ncpu = 0;
                break;
            case 1002:  /* --seed */
                seed = (unsigned int)atoi(optarg);
                break;
            case 1003:  /* --stratified */
                stratified = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_split_usage(argv[0]);
                return 0;
            default:
                print_split_usage(argv[0]);
                return 1;
        }
    }

    /* Validate required arguments */
    if (!input_file) {
        fprintf(stderr, "Error: Input file required (-f/--file)\n");
        return 1;
    }
    if (!output_str) {
        fprintf(stderr, "Error: Output files required (-o/--output)\n");
        return 1;
    }
    if (!smiles_col) {
        fprintf(stderr, "Error: SMILES column name required (-s/--col)\n");
        return 1;
    }

    /* Parse output files */
    char* output_files_buf[SPLITTER_MAX_SPLITS];
    int num_outputs = 0;

    char* output_copy = strdup(output_str);
    if (!output_copy) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }

    char* token = strtok(output_copy, ",");
    while (token && num_outputs < SPLITTER_MAX_SPLITS) {
        output_files_buf[num_outputs++] = strdup(token);
        token = strtok(NULL, ",");
    }
    free(output_copy);

    if (num_outputs < 2) {
        fprintf(stderr, "Error: At least 2 output files required\n");
        for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
        return 1;
    }

    /* Parse splitting method */
    split_method_t method = SPLIT_METHOD_RANDOM;
    if (strcmp(method_str, "scaffold") == 0) {
        method = SPLIT_METHOD_SCAFFOLD;
    } else if (strcmp(method_str, "random") != 0) {
        fprintf(stderr, "Error: Unknown splitting method '%s'. Use 'random' or 'scaffold'.\n", method_str);
        for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
        return 1;
    }

    /* Set up split options */
    split_options_t options = SPLIT_OPTIONS_DEFAULT;
    options.method = method;
    options.num_splits = num_outputs;
    options.num_threads = ncpu;
    options.seed = seed;
    options.verbose = verbose;
    options.stratified = stratified;

    /* Parse ratios */
    if (ratio_str) {
        int num_ratios;
        if (splitter_parse_ratios(ratio_str, options.ratios, &num_ratios) != CCHEM_OK) {
            fprintf(stderr, "Error: Invalid split ratios '%s'\n", ratio_str);
            for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
            return 1;
        }
        if (num_ratios != num_outputs) {
            fprintf(stderr, "Error: Number of ratios (%d) must match number of output files (%d)\n",
                    num_ratios, num_outputs);
            for (int i = 0; i < num_outputs; i++) free(output_files_buf[i]);
            return 1;
        }
    } else {
        /* Default ratios based on number of outputs */
        if (num_outputs == 2) {
            options.ratios[0] = 0.8;
            options.ratios[1] = 0.2;
        } else if (num_outputs == 3) {
            options.ratios[0] = 0.8;
            options.ratios[1] = 0.1;
            options.ratios[2] = 0.1;
        } else {
            /* Equal splits for more than 3 */
            for (int i = 0; i < num_outputs; i++) {
                options.ratios[i] = 1.0 / num_outputs;
            }
        }
    }

    /* Run the split */
    char error_buf[512];
    const char* output_files[SPLITTER_MAX_SPLITS];
    for (int i = 0; i < num_outputs; i++) {
        output_files[i] = output_files_buf[i];
    }

    cchem_status_t status = splitter_split_csv(
        input_file, output_files, num_outputs, smiles_col,
        &options, error_buf, sizeof(error_buf));

    /* Cleanup */
    for (int i = 0; i < num_outputs; i++) {
        free(output_files_buf[i]);
    }

    if (status != CCHEM_OK) {
        fprintf(stderr, "Error: %s\n", error_buf);
        return 1;
    }

    return 0;
}
