# cchem Benchmarks

Performance benchmarks comparing cchem against other cheminformatics libraries.

## Canonicalization Benchmark

Two benchmark implementations are available:
1. **Python benchmark** - Uses Python bindings (RDKit, OpenBabel)
2. **C++ benchmark** - Direct C/C++ API calls for fairer comparison

Compares SMILES canonicalization speed across:
- **cchem** (this library)
- **RDKit** (industry standard)
- **OpenBabel** (open-source toolkit)
- **CDK** (Chemistry Development Kit, Java - Python benchmark only)

## Results Summary (27k molecules)

### C/C++ Native Performance

| Library | 1-thread | 8-thread | Scaling |
|---------|----------|----------|---------|
| **cchem** | 33,210/s | 139,265/s | 4.19x |
| RDKit | 10,617/s | 48,740/s | 4.59x |
| OpenBabel | 9,611/s | n/a | not thread-safe |

**cchem is 3.1x faster single-threaded and 2.9x faster multi-threaded vs RDKit**

### Python Bindings Performance

| Library | 1-thread | 8-thread | Scaling |
|---------|----------|----------|---------|
| **cchem** | 32,006/s | 124,833/s | 3.90x |
| RDKit | 11,320/s | 21,721/s | 1.92x |
| OpenBabel | 9,137/s | 20,008/s | 2.19x |

Python multiprocessing has higher overhead than native threading.

## C++ Benchmark

### Build

```bash
cd build

# Build with RDKit only (recommended)
cmake .. -DBUILD_BENCHMARKS=ON -DWITH_RDKIT=ON
make canonicalization_benchmark

# Build with RDKit and OpenBabel
cmake .. -DBUILD_BENCHMARKS=ON -DWITH_RDKIT=ON -DWITH_OPENBABEL=ON
make canonicalization_benchmark
```

### Run

```bash
./benchmarks/canonicalization_benchmark -f data.csv -s smiles -n 8 -i 3
```

### Requirements

```bash
# macOS
brew install rdkit open-babel

# Linux
# See RDKit/OpenBabel documentation for installation
```

## Python Benchmark

### Quick Start

```bash
# Basic usage
python benchmarks/canonicalization_benchmark.py -f molecules.csv -s smiles

# With multiple threads
python benchmarks/canonicalization_benchmark.py -f molecules.csv -s smiles -n 8

# Using shell wrapper
./benchmarks/run_benchmark.sh molecules.csv smiles 8
```

### Installation

Install comparison libraries (optional):

```bash
# RDKit (recommended)
pip install rdkit

# OpenBabel
pip install openbabel

# CDK (requires Java)
pip install jpype1
export CDK_JAR=/path/to/cdk-2.9.jar
```

### CLI Options

```
usage: canonicalization_benchmark.py [-h] -f FILE [-s SMILES_COL] [-n THREADS]
                                     [--iterations N] [--warmup N] [--limit N]
                                     [--baseline LIB] [--cchem-path PATH]
                                     [--single-mode] [--skip-rdkit]
                                     [--skip-openbabel] [--skip-cdk] [-v]

Arguments:
  -f, --file FILE         Input CSV file with SMILES (required)
  -s, --smiles-col COL    Column name containing SMILES (default: smiles)
  -n, --threads N         Number of threads for cchem batch mode (default: 1)
  --iterations N          Number of benchmark iterations (default: 3)
  --warmup N              Number of warmup iterations (default: 1)
  --limit N               Limit number of molecules to process
  --baseline LIB          Library to use as baseline for speedup (default: rdkit)
  --cchem-path PATH       Path to cchem binary
  --single-mode           Also benchmark cchem single-molecule mode
  --skip-rdkit            Skip RDKit benchmark
  --skip-openbabel        Skip OpenBabel benchmark
  --skip-cdk              Skip CDK benchmark
  -v, --verbose           Verbose output
```

### Example Output

```
================================================================================
CANONICALIZATION BENCHMARK RESULTS
================================================================================

Library              Time (s)     Throughput      Avg (us)     Failed   Speedup
--------------------------------------------------------------------------------
cchem (n=8)          0.234        427,350         2.34         0        4.27x
cchem (n=1)          1.021        97,942          10.21        0        0.98x
RDKit                1.012        100,000         10.00        12       1.00x
OpenBabel            3.456        28,935          34.56        5        0.29x
--------------------------------------------------------------------------------
Baseline: RDKit

Dataset: molecules.csv
Molecules: 100,000
Iterations: 3
```

### Metrics

- **Time (s)**: Total wall-clock time to canonicalize all molecules
- **Throughput**: Molecules processed per second
- **Avg (us)**: Average time per molecule in microseconds
- **Failed**: Number of molecules that failed to canonicalize
- **Speedup**: Relative to baseline library (default: RDKit)

### Sample Datasets

You can use publicly available molecular datasets:

```bash
# ChEMBL subset
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_33/chembl_33_chemreps.txt.gz
gunzip chembl_33_chemreps.txt.gz

# ZINC subset
# Download from https://zinc.docking.org/

# PubChem subset
# Download from https://pubchem.ncbi.nlm.nih.gov/
```

### Notes

- cchem uses native pthreads for parallel execution
- Python RDKit/OpenBabel use multiprocessing (process pool) to work around GIL
- C++ RDKit uses native pthreads and scales well
- C++ OpenBabel is not thread-safe (single-threaded only)
- Results may vary based on molecule complexity and system configuration
