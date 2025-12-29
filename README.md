# cchem

<img width="197" height="168" alt="cchem_logo" src="https://github.com/user-attachments/assets/461c8341-390e-4568-a803-a4137342f2b8" />

High-performance cheminformatics library written in pure C.

**Key Features:**
- 1600+ molecular descriptors
- SMILES canonicalization with stereochemistry support
- Molecular sanitization (salt removal, aromatization, neutralization, tautomers)
- 2D/3D molecular visualization with MMFF94 force field
- Multi-threaded batch processing with SIMD optimization

## Features

### SMILES Canonicalization
- Morgan algorithm for canonical ordering
- Stereochemistry preservation (@ chirality, E/Z double bonds)
- Smallest Set of Smallest Rings (SSSR) detection
- Aromaticity perception
- Isotope and hydrogen handling

### Molecular Sanitization
- **Salt removal**: Identify and remove counter-ions (Na+, K+, Cl-, etc.)
- **Aromatization**: Perceive aromaticity using Hückel's rule (4n+2 π electrons)
- **Kekulization**: Convert aromatic bonds to alternating single/double
- **Neutralization**: Remove protonation states (R-NH3+ → R-NH2, R-O- → R-OH)
- **Normalization**: Standardize functional groups (nitro, sulfoxide, phosphate)
- **Tautomer enumeration**: Generate keto-enol, amide-imidic acid tautomers
- **Cleanup options**: Remove stereochemistry, isotopes, explicit hydrogens

### Molecular Descriptors (1600+)

| Category | Count | Description |
|----------|-------|-------------|
| Counts | 83 | Element counts (C, H, N, O, S, halogens), bond types, ring counts |
| Ratios | 30+ | Elemental ratios, hybridization ratios, electronegativity ratios |
| Topology | 20+ | Kier-Hall Chi indices, Zagreb indices, Wiener index, Balaban J |
| Electronic | 30+ | Gasteiger-Marsili charges (PEOE), electrotopological states |
| Steric | 20+ | Van der Waals volume, McGowan volume, TPSA |
| Energetic | 30+ | Born solvation proxies, FMO hardness, Hansen parameters |
| Fractional | 60+ | Molecular weight fractions, bond fractions |
| Hash | 30+ | SMILES hashing, n-gram hashes, MinHash signatures |
| Graph | 30+ | Density, centrality, clustering, spectral properties |
| Autocorrelation | 54 | Broto-Moreau 2D autocorrelations (ATS) lags 0-8 |
| Solubility | 1 | CLogS aqueous solubility |
| LogP/LogD | 2 | Wildman-Crippen LogP, Neural network LogD 7.4 (R²=0.94) |
| MQN | 42 | Molecular Quantum Numbers |
| VSA | 47 | SlogP_VSA, SMR_VSA, PEOE_VSA, EState_VSA |
| BCUT | 48 | Burden-CAS-University of Texas eigenvalues |
| Zagreb | 24 | Zagreb indices and variants |
| Information | 24 | Information content descriptors |
| Walk Counts | 36 | Molecular walk counts |
| E-State Sums | 32 | Electrotopological state sums |
| ETA | 24 | Extended topochemical atom indices |
| Ring Complexity | 18 | Ring system complexity measures |
| CPSA | 70 | Charged partial surface area |
| Moments | 42 | Molecular geometry moments |
| Aromatic | 64 | Aromatic system descriptors |
| Atom Pairs | 56 | Distance-based atom pair fingerprints |
| Framework | 40 | Molecular framework descriptors |
| Constitutional | 34 | Constitutional descriptors |
| Functional | 50+ | Carbonyl, nitrogen, sulfur, oxygen, heterocyclic scaffolds |
| Pharmacophore | 30+ | Pharmacophore points, density, drug-likeness metrics |

### 2D/3D Visualization
- 2D coordinate generation with automatic layout
- 3D coordinate generation with MMFF94 force field optimization
- Render styles: wireframe, sticks, ball-and-stick, spacefill, surface
- JPEG output with configurable quality
- Atom coloring by element (CPK scheme)

### Dataset Utilities
- CSV batch processing with parallel execution
- Dataset splitting: random, scaffold-based (Murcko), stratified
- Configurable train/validation/test ratios

### Performance
- SIMD optimization (`-march=native`)
- Multi-threaded processing (pthreads)
- Link Time Optimization (LTO)
- Arena allocators and thread-local caches
- Batch compute functions per descriptor category
- Pipeline streaming for constant memory on large datasets

## Installation

### Dependencies

| Library | Purpose |
|---------|---------|
| cairo | 2D vector graphics rendering |
| libjpeg | JPEG image output |
| pthreads | Parallel processing |
| libm | Mathematics functions |

### macOS (Homebrew)

```bash
brew install cairo jpeg
```

### Linux (apt)

```bash
sudo apt-get install libcairo2-dev libjpeg-dev
```

### Build

```bash
git clone https://github.com/your-repo/cchem.git
cd cchem
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Quick Start

### Canonicalize a SMILES

```bash
# Single molecule
./cchem canonicalize -S "c1ccccc1"
# Output: c1ccccc1

# Batch processing
./cchem canonicalize -f molecules.csv -s smiles -c canonical -o output.csv -n 4
```

### Sanitize Molecules

```bash
# Complete sanitization (unsalt + aromatize + neutralize + normalize)
./cchem canonicalize -S "[Na+].CC(=O)[O-]" --sanitize complete
# Output: C(=O)(O)C  (acetic acid, Na+ removed, carboxylate neutralized)

# Remove salts only
./cchem canonicalize -S "CCO.[Cl-].[Na+]" --sanitize unsalt
# Output: CCO

# Aromatize Kekule form
./cchem canonicalize -S "C1=CC=CC=C1" --sanitize aromatize
# Output: c1ccccc1

# Multiple operations
./cchem canonicalize -S "[NH3+]Cc1ccccc1.[Cl-]" --sanitize unsalt,neutralize,aromatize
# Output: c1ccc(cc1)C[NH2]  (benzylamine)

# List tautomers
./cchem canonicalize -S "CC(=O)CC" --list-tautomers
# Output: CCC(=O)C  (keto-enol tautomers)
```

### Compute Descriptors

```bash
# Single molecule with specific descriptors
./cchem compute -S "CCO" -d CarbonCount,HydrogenCount,MolecularWeight

# All descriptors for a dataset
./cchem compute -f data.csv -s smiles -d all -o descriptors.csv -n 8

# List available descriptors
./cchem compute --list
```

### Generate Molecular Images

```bash
# 2D depiction
./cchem depict -S "c1ccccc1" -o benzene.jpg -m 2d

# 3D with MMFF94 optimization
./cchem depict -S "CCO" -o ethanol.jpg -m 3d -s balls-sticks --max-iter 500
```

### Split Dataset

```bash
# Train/test split (80/20)
./cchem split -f data.csv -s smiles -o train.csv,test.csv

# Train/validation/test with scaffold splitting
./cchem split -f data.csv -s smiles -o train.csv,val.csv,test.csv \
    --split-ratios 80,10,10 --splitting-method scaffold
```

## CLI Reference

### canonicalize

Convert SMILES to canonical form (aromatized by default) with optional sanitization.

```
Usage: cchem canonicalize [options]

Options:
  -S, --smiles <string>     Single SMILES string to canonicalize
  -f, --file <path>         Input CSV file
  -s, --smiles-col <name>   Column name containing SMILES (default: smiles)
  -c, --canon-col <name>    Output column name for canonical SMILES
  -o, --output <path>       Output file path
  -n, --threads <num>       Number of threads (default: auto)
  --sanitize <opts>         Apply sanitization before canonicalization
                            Values: "complete" or comma-separated list of:
                              unsalt          - Remove salts, keep largest fragment
                              aromatize       - Perceive and apply aromaticity (default)
                              kekulize        - Convert aromatic to Kekule form
                              neutralize      - Neutralize charges
                              normalize       - Normalize functional groups
                              remove-stereo   - Remove stereochemistry
                              remove-isotopes - Remove isotope labels
                              remove-h        - Remove explicit hydrogens
                              validate        - Validate structure
  --list-tautomers          List all tautomeric forms
  -v, --verbose             Enable verbose output
  -h, --help                Print this help message
```

### compute

Calculate molecular descriptors.

```
Usage: cchem compute [options]

Options:
  -S, --smiles <string>     Single SMILES string
  -f, --file <path>         Input CSV file
  -s, --smiles-col <name>   Column name containing SMILES (default: smiles)
  -d, --descriptors <list>  Comma-separated descriptor names or "all"
  -o, --output <path>       Output file path
  -n, --threads <num>       Number of threads (default: 1)
  -l, --list                List all available descriptors
  --no-canonicalization     Skip SMILES canonicalization
  -v, --verbose             Enable verbose output
```

### depict

Generate 2D/3D molecular structure images.

```
Usage: cchem depict [options]

Options:
  -S, --smiles <string>     SMILES string to depict
  -o, --output <path>       Output JPEG file path
  -m, --mode <2d|3d>        Rendering mode (default: 2d)
  -s, --style <style>       Render style:
                              wireframe    - Simple lines
                              sticks       - Colored bond sticks
                              balls-sticks - Ball and stick model
                              spacefill    - CPK space-filling
                              surface      - Molecular surface
  -W, --width <pixels>      Image width (default: 800)
  -H, --height <pixels>     Image height (default: 800)
  --bond-length <float>     Bond length in pixels (default: 35)
  --bond-width <float>      Bond width in pixels (default: 2)
  --margin <pixels>         Image margin (default: 20)
  --show-carbons            Show carbon atom labels
  --show-hydrogens          Show hydrogen atoms
  --toggle-aromaticity      Toggle aromatic ring display
  --atom-filling <float>    Atom sphere filling (0.0-1.0)
  --quality <int>           JPEG quality (1-100, default: 90)
  --max-iter <int>          Max MMFF94 iterations (default: 200)
  --surface-color <hex>     Surface color (e.g., 0x808080)
```

### split

Split datasets for machine learning.

```
Usage: cchem split [options]

Options:
  -f, --file <path>              Input CSV file
  -o, --output <paths>           Comma-separated output file paths
  -s, --smiles-col <name>        Column containing SMILES (default: smiles)
  --split-ratios <ratios>        Comma-separated percentages (default: 80,20)
  --splitting-method <method>    Splitting method:
                                   random   - Random split
                                   scaffold - Murcko scaffold-based
  --stratified                   Use stratified splitting
  --seed <int>                   Random seed for reproducibility
```

### validate

Check SMILES syntax validity.

```
Usage: cchem validate [options]

Options:
  -S, --smiles <string>     SMILES string to validate
  -v, --verbose             Show detailed validation info
```

### version / help

```bash
cchem version    # Show version information
cchem help       # Show help message
```

## Building from Source

### Requirements

- CMake 3.16+
- C11 compatible compiler (GCC, Clang)
- Dependencies: igraph, cairo, libjpeg

### Build Commands

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Build Options

```bash
# Debug build with sanitizers
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release build with LTO
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Running Tests

```bash
cd build
ctest                    # Run all tests
ctest -V                 # Verbose output
ctest -R test_canon      # Run specific test
```

## Project Structure

```
cchem/
├── include/cchem/
│   ├── cchem.h              # Main library header
│   ├── descriptors.h        # Descriptor system
│   ├── canonicalizer/       # SMILES parsing & canonicalization
│   │   ├── sanitize.h       # Sanitization API
│   │   └── ...              # Parser, molecule, stereo, etc.
│   └── depictor/            # 2D/3D visualization
├── src/
│   ├── main.c               # CLI entry point
│   ├── canonicalizer/       # SMILES implementation
│   │   ├── sanitize.c       # Sanitization implementation
│   │   └── ...
│   ├── descriptors/         # 500+ descriptor implementations
│   └── depictor/            # Visualization & MMFF94
├── tests/                   # Test suite
├── data/                    # Training data & test files
└── doc/                     # Documentation
```

## License

Licensed under the Apache License, Version 2.0.

```
Copyright 2024

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```
