#!/bin/bash
# Convenience wrapper for running canonicalization benchmarks
#
# Usage:
#   ./run_benchmark.sh molecules.csv smiles
#   ./run_benchmark.sh molecules.csv smiles 8  # with 8 threads

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_SCRIPT="$SCRIPT_DIR/canonicalization_benchmark.py"

if [ $# -lt 1 ]; then
    echo "Usage: $0 <csv_file> [smiles_column] [threads] [extra_args...]"
    echo ""
    echo "Arguments:"
    echo "  csv_file       Input CSV file with SMILES"
    echo "  smiles_column  Column name containing SMILES (default: smiles)"
    echo "  threads        Number of threads for cchem (default: 1)"
    echo ""
    echo "Examples:"
    echo "  $0 molecules.csv"
    echo "  $0 molecules.csv smiles 8"
    echo "  $0 data.csv SMILES 4 --iterations 5"
    exit 1
fi

CSV_FILE="$1"
SMILES_COL="${2:-smiles}"
THREADS="${3:-1}"

# Shift away the first 3 positional args to pass remaining to Python
shift 3 2>/dev/null || shift $#

if [ ! -f "$CSV_FILE" ]; then
    echo "Error: File not found: $CSV_FILE"
    exit 1
fi

# Check for Python
PYTHON=""
for cmd in python3 python; do
    if command -v "$cmd" &>/dev/null; then
        PYTHON="$cmd"
        break
    fi
done

if [ -z "$PYTHON" ]; then
    echo "Error: Python not found"
    exit 1
fi

# Find cchem binary
CCHEM_PATH=""
for path in "$SCRIPT_DIR/../build/cchem" "$SCRIPT_DIR/../cchem" "$(which cchem 2>/dev/null)"; do
    if [ -x "$path" ]; then
        CCHEM_PATH="$path"
        break
    fi
done

if [ -z "$CCHEM_PATH" ]; then
    echo "Error: cchem binary not found. Build with: cd build && cmake .. && make"
    exit 1
fi

echo "Running benchmark..."
echo "  CSV file: $CSV_FILE"
echo "  SMILES column: $SMILES_COL"
echo "  Threads: $THREADS"
echo ""

exec "$PYTHON" "$BENCHMARK_SCRIPT" \
    -f "$CSV_FILE" \
    -s "$SMILES_COL" \
    -n "$THREADS" \
    --cchem-path "$CCHEM_PATH" \
    "$@"
