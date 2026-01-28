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

**Test System:** MacBook Air M3 (2024), 8-core CPU @ 4.06 GHz, 16 GB RAM, macOS Tahoe 26.2

### C/C++ Native Performance

| Library | 1-thread | vs RDKit | 8-thread | vs RDKit |
|---------|----------|----------|----------|----------|
| **cchem** | 30,052/s | **2.82x** | 132,335/s | **2.66x** |
| RDKit | 10,654/s | 1.00x | 49,778/s | 1.00x |
| OpenBabel | 9,733/s | 0.91x | n/a | not thread-safe |

### Python Bindings Performance

| Library | 1-thread | vs RDKit | 8-thread | vs RDKit |
|---------|----------|----------|----------|----------|
| CDK | 32,609/s | **2.90x** | n/a | single-threaded only |
| **cchem** | 28,311/s | **2.52x** | 112,815/s | **9.87x** |
| RDKit | 11,229/s | 1.00x | 11,434/s | 1.00x |
| OpenBabel | 9,092/s | 0.81x | 10,561/s | 0.92x |

*Note: CDK (Java) has excellent single-threaded performance due to JVM JIT optimization. RDKit/OpenBabel multi-threading limited by Python GIL.*

### Parallel Scaling Efficiency

| Library | 1→8 threads | Efficiency |
|---------|-------------|------------|
| **cchem** | 3.98x | 50% |
| OpenBabel | 1.16x | 15% |
| RDKit | 1.02x | 13% |

cchem achieves better parallel scaling due to native pthreads vs Python multiprocessing overhead.

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
