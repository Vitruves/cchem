/**
 * @file cmd_depict.c
 * @brief Depict command implementation
 */

#include "cchem/utils/commands.h"

void print_depict_usage(const char* prog_name) {
    printf("Usage: %s depict [options]\n\n", prog_name);
    printf("Generate publication-quality molecular structure images from SMILES\n\n");
    printf("Required:\n");
    printf("  -S, --smiles <SMILES>   SMILES string to depict\n");
    printf("  -o, --output <file>     Output file (.png, .jpg, .svg)\n");
    printf("\n");
    printf("Render Style:\n");
    printf("  -s, --style <style>     wireframe (default), sticks, balls-sticks, spacefill, surface\n");
    printf("  -m, --mode <mode>       2d (default) or 3d (MMFF94 optimized geometry)\n");
    printf("\n");
    printf("Dimensions:\n");
    printf("  -W, --width <pixels>    Image width (default: 800)\n");
    printf("  -H, --height <pixels>   Image height (default: 800)\n");
    printf("  --margin <px>           Margin (default: 50)\n");
    printf("  --scale <factor>        Resolution multiplier (e.g., 2.0 for 2x)\n");
    printf("\n");
    printf("Bond Rendering:\n");
    printf("  --bond-length <px>      Bond length (default: 35)\n");
    printf("  --bond-width <px>       Bond line width (default: 1.8)\n");
    printf("  --heteroatom-gap <0-1>  Gap at heteroatoms (default: 1.0, 0 to disable)\n");
    printf("  --line-cap <style>      round, butt (default), square\n");
    printf("\n");
    printf("Atom Labels:\n");
    printf("  --show-carbons          Show 'C' labels on carbons\n");
    printf("  --show-hydrogens        Show hydrogen labels\n");
    printf("  --terminal-carbons      Show CH3 on terminal carbons\n");
    printf("  --font-size <scale>     Font size scale (default: 3.5)\n");
    printf("\n");
    printf("3D/Surface Options:\n");
    printf("  --atom-filling          Draw atoms as filled CPK spheres\n");
    printf("  --proportional-atoms    Scale atoms by VDW radius (default: on)\n");
    printf("  --no-proportional       Disable proportional sizing\n");
    printf("  --surface-color <mode>  uniform, atom, polarity (for -s surface)\n");
    printf("  --max-iter <N>          3D optimization iterations (default: 500)\n");
    printf("\n");
    printf("Other:\n");
    printf("  --toggle-aromaticity    Draw aromatic circles\n");
    printf("  --quality <1-100>       JPEG quality (default: 95)\n");
    printf("  -v, --verbose           Verbose output\n");
    printf("  --debug                 Print debug information\n");
    printf("  -h, --help              Print this help\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s depict -S \"c1ccccc1\" -o benzene.png\n", prog_name);
    printf("  %s depict -S \"CCO\" -o ethanol.svg\n", prog_name);
    printf("  %s depict -S \"CCO\" -o ethanol.png -m 3d -s balls-sticks\n", prog_name);
    printf("  %s depict -S \"CCO\" -o ethanol.png --scale 2\n", prog_name);
}

int cmd_depict(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"smiles",         required_argument, 0, 'S'},
        {"output",         required_argument, 0, 'o'},
        {"mode",           required_argument, 0, 'm'},
        {"style",          required_argument, 0, 's'},
        {"width",          required_argument, 0, 'W'},
        {"height",         required_argument, 0, 'H'},
        {"bond-length",    required_argument, 0, 1000},
        {"bond-width",     required_argument, 0, 1001},
        {"margin",         required_argument, 0, 1002},
        {"show-carbons",   no_argument,       0, 1003},
        {"show-hydrogens", no_argument,       0, 1004},
        {"toggle-aromaticity", no_argument,   0, 1005},
        {"atom-filling",   no_argument,       0, 1006},
        {"font-size",      required_argument, 0, 1007},
        {"quality",        required_argument, 0, 1008},
        {"max-iter",       required_argument, 0, 1009},
        {"proportional-atoms", no_argument,   0, 1010},
        {"no-proportional", no_argument,      0, 1011},
        {"surface-color",  required_argument, 0, 1012},
        {"modern",         no_argument,       0, 1013},
        {"heteroatom-gap", required_argument, 0, 1014},
        {"scale",          required_argument, 0, 1015},
        {"line-cap",       required_argument, 0, 1016},
        {"terminal-carbons", no_argument,     0, 1017},
        {"verbose",        no_argument,       0, 'v'},
        {"debug",          no_argument,       0, 1018},
        {"help",           no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    char* smiles = NULL;
    char* output_file = NULL;
    bool verbose = false;

    depictor_options_t options = DEPICTOR_OPTIONS_DEFAULT;

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "S:o:m:s:W:H:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'S':
                smiles = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'm':
                if (strcmp(optarg, "2d") == 0 || strcmp(optarg, "2D") == 0) {
                    options.mode = DEPICT_MODE_2D;
                } else if (strcmp(optarg, "3d") == 0 || strcmp(optarg, "3D") == 0) {
                    options.mode = DEPICT_MODE_3D;
                } else {
                    fprintf(stderr, "Error: Invalid mode '%s'. Use '2d' or '3d'.\n", optarg);
                    return 1;
                }
                break;
            case 's':  /* --style */
                if (strcmp(optarg, "wireframe") == 0) {
                    options.render_style = RENDER_STYLE_WIREFRAME;
                } else if (strcmp(optarg, "sticks") == 0) {
                    options.render_style = RENDER_STYLE_STICKS;
                } else if (strcmp(optarg, "balls-sticks") == 0 || strcmp(optarg, "balls") == 0) {
                    options.render_style = RENDER_STYLE_BALLS_AND_STICKS;
                } else if (strcmp(optarg, "spacefill") == 0 || strcmp(optarg, "cpk") == 0) {
                    options.render_style = RENDER_STYLE_SPACEFILL;
                } else if (strcmp(optarg, "surface") == 0) {
                    options.render_style = RENDER_STYLE_SURFACE;
                } else {
                    fprintf(stderr, "Error: Invalid style '%s'. Use wireframe, sticks, balls-sticks, spacefill, or surface.\n", optarg);
                    return 1;
                }
                break;
            case 'W':
                options.width = atoi(optarg);
                if (options.width < 50 || options.width > 4096) {
                    fprintf(stderr, "Error: Width must be between 50 and 4096.\n");
                    return 1;
                }
                break;
            case 'H':
                options.height = atoi(optarg);
                if (options.height < 50 || options.height > 4096) {
                    fprintf(stderr, "Error: Height must be between 50 and 4096.\n");
                    return 1;
                }
                break;
            case 1000:  /* --bond-length */
                options.bond_length = atof(optarg);
                break;
            case 1001:  /* --bond-width */
                options.bond_width = atof(optarg);
                break;
            case 1002:  /* --margin */
                options.margin = atoi(optarg);
                break;
            case 1003:  /* --show-carbons */
                options.show_carbons = true;
                break;
            case 1004:  /* --show-hydrogens */
                options.show_hydrogens = true;
                break;
            case 1005:  /* --toggle-aromaticity */
                options.draw_aromatic_circles = true;
                break;
            case 1006:  /* --atom-filling */
                options.atom_filling = true;
                break;
            case 1007:  /* --font-size */
                options.font_size = atof(optarg);
                if (options.font_size < 0.5) options.font_size = 0.5;
                if (options.font_size > 10.0) options.font_size = 10.0;
                break;
            case 1008:  /* --quality */
                options.jpeg_quality = atoi(optarg);
                if (options.jpeg_quality < 1) options.jpeg_quality = 1;
                if (options.jpeg_quality > 100) options.jpeg_quality = 100;
                break;
            case 1009:  /* --max-iter */
                options.max_iterations = atoi(optarg);
                break;
            case 1010:  /* --proportional-atoms */
                options.proportional_atoms = true;
                break;
            case 1011:  /* --no-proportional */
                options.proportional_atoms = false;
                break;
            case 1012:  /* --surface-color */
                if (strcmp(optarg, "uniform") == 0 || strcmp(optarg, "blue") == 0) {
                    options.surface_color = SURFACE_COLOR_UNIFORM;
                } else if (strcmp(optarg, "atom") == 0 || strcmp(optarg, "element") == 0) {
                    options.surface_color = SURFACE_COLOR_ATOM;
                } else if (strcmp(optarg, "polarity") == 0 || strcmp(optarg, "charge") == 0) {
                    options.surface_color = SURFACE_COLOR_POLARITY;
                } else {
                    fprintf(stderr, "Error: Invalid surface color '%s'. Use uniform, atom, or polarity.\n", optarg);
                    return 1;
                }
                break;
            case 1013:  /* --modern - now default, kept for compatibility */
                /* Already default, no changes needed */
                break;
            case 1014:  /* --heteroatom-gap */
                options.heteroatom_gap = atof(optarg);
                if (options.heteroatom_gap < 0.0) options.heteroatom_gap = 0.0;
                if (options.heteroatom_gap > 1.0) options.heteroatom_gap = 1.0;
                break;
            case 1015:  /* --scale */
                options.scale_factor = atof(optarg);
                if (options.scale_factor < 0.1) options.scale_factor = 0.1;
                if (options.scale_factor > 10.0) options.scale_factor = 10.0;
                break;
            case 1016:  /* --line-cap */
                if (strcmp(optarg, "round") == 0) {
                    options.line_cap = LINE_CAP_ROUND;
                } else if (strcmp(optarg, "butt") == 0) {
                    options.line_cap = LINE_CAP_BUTT;
                } else if (strcmp(optarg, "square") == 0) {
                    options.line_cap = LINE_CAP_SQUARE;
                } else {
                    fprintf(stderr, "Error: Invalid line cap '%s'. Use round, butt, or square.\n", optarg);
                    return 1;
                }
                break;
            case 1017:  /* --terminal-carbons */
                options.terminal_carbon_labels = true;
                break;
            case 1018:  /* --debug */
                options.debug = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_depict_usage(argv[0]);
                return 0;
            default:
                print_depict_usage(argv[0]);
                return 1;
        }
    }

    if (!smiles) {
        fprintf(stderr, "Error: SMILES string required (-S/--smiles)\n");
        print_depict_usage(argv[0]);
        return 1;
    }

    if (!output_file) {
        fprintf(stderr, "Error: Output file required (-o/--output)\n");
        print_depict_usage(argv[0]);
        return 1;
    }

    /* Auto-detect output format from file extension */
    const char* ext = strrchr(output_file, '.');
    if (ext) {
        if (strcasecmp(ext, ".svg") == 0) {
            options.format = IMG_FORMAT_SVG;
        } else if (strcasecmp(ext, ".png") == 0) {
            options.format = IMG_FORMAT_PNG;
        } else if (strcasecmp(ext, ".jpg") == 0 || strcasecmp(ext, ".jpeg") == 0) {
            options.format = IMG_FORMAT_JPEG;
        }
    }

    /* Get render style name */
    const char* style_name = "sticks";
    switch (options.render_style) {
        case RENDER_STYLE_WIREFRAME: style_name = "wireframe"; break;
        case RENDER_STYLE_STICKS: style_name = "sticks"; break;
        case RENDER_STYLE_BALLS_AND_STICKS: style_name = "balls-and-sticks"; break;
        case RENDER_STYLE_SPACEFILL: style_name = "spacefill"; break;
        case RENDER_STYLE_SURFACE: style_name = "surface"; break;
    }

    char error_buf[256];
    depict_info_t info = {0};
    cchem_status_t status;

    if (verbose) {
        status = depict_smiles_verbose(smiles, output_file, &options,
                                       &info, error_buf, sizeof(error_buf));
    } else {
        status = depict_smiles(smiles, output_file, &options,
                               error_buf, sizeof(error_buf));
    }

    if (status == CCHEM_OK) {
        if (verbose) {
            printf("-- Depiction Summary --\n");
            printf("Input SMILES:     %s\n", smiles);
            printf("Canonical SMILES: %s\n", info.canonical_smiles);
            printf("Output file:      %s\n", output_file);
            printf("Mode:             %s\n", options.mode == DEPICT_MODE_2D ? "2D" : "3D (MMFF94)");
            printf("Render style:     %s\n", style_name);
            printf("Dimensions:       %dx%d\n", options.width, options.height);
            printf("Molecule:         %d atoms, %d bonds, %d rings\n",
                   info.num_atoms, info.num_bonds, info.num_rings);
            if (options.mode == DEPICT_MODE_3D) {
                printf("MMFF94 Energy:    %.2f -> %.2f kcal/mol (delta: %.2f)\n",
                       info.energy_initial, info.energy_final,
                       info.energy_initial - info.energy_final);
            }
        }
        return 0;
    } else {
        fprintf(stderr, "Error: %s\n", error_buf);
        return 1;
    }
}
