#!/usr/bin/env python3
"""
Benchmark script for comparing SMILES canonicalization speed across cheminformatics libraries.

Compares both single-threaded and multi-threaded performance:
- cchem (this library)
- RDKit
- OpenBabel (pybel)
- CDK (via jpype, optional)

Usage:
    python canonicalization_benchmark.py -f molecules.csv -s smiles_column
    python canonicalization_benchmark.py -f molecules.csv -s smiles -n 8 --iterations 3
"""

import argparse
import csv
import multiprocessing as mp
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Optional

# Try to import optional dependencies
RDKIT_AVAILABLE = False
OPENBABEL_AVAILABLE = False
CDK_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')  # Suppress RDKit warnings
    RDKIT_AVAILABLE = True
except ImportError:
    pass

try:
    from openbabel import pybel
    pybel.ob.obErrorLog.SetOutputLevel(0)  # Suppress warnings
    OPENBABEL_AVAILABLE = True
except ImportError:
    pass

try:
    import jpype
    import jpype.imports
    if not jpype.isJVMStarted():
        cdk_jar = os.environ.get('CDK_JAR', None)
        if cdk_jar and os.path.exists(cdk_jar):
            jpype.startJVM(classpath=[cdk_jar])
            from org.openscience.cdk.smiles import SmilesParser, SmilesGenerator
            from org.openscience.cdk.silent import SilentChemObjectBuilder
            CDK_AVAILABLE = True
except (ImportError, Exception):
    pass


@dataclass
class BenchmarkResult:
    """Result of a benchmark run."""
    library: str
    total_time: float
    molecules_processed: int
    molecules_failed: int
    throughput: float  # molecules per second
    avg_time_per_mol: float  # microseconds
    threads: int = 1


def load_smiles_from_csv(filepath: str, smiles_col: str, limit: Optional[int] = None) -> list[str]:
    """Load SMILES strings from a CSV file."""
    smiles_list = []
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if smiles_col not in reader.fieldnames:
            available = ', '.join(reader.fieldnames)
            raise ValueError(f"Column '{smiles_col}' not found. Available columns: {available}")

        for i, row in enumerate(reader):
            if limit and i >= limit:
                break
            smi = row[smiles_col].strip()
            if smi:
                smiles_list.append(smi)

    return smiles_list


# =============================================================================
# cchem benchmarks
# =============================================================================

def benchmark_cchem(smiles_list: list[str], cchem_path: str, threads: int = 1) -> BenchmarkResult:
    """Benchmark cchem canonicalization."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_in:
        tmp_in.write('smiles\n')
        for smi in smiles_list:
            tmp_in.write(f'{smi}\n')
        tmp_in_path = tmp_in.name

    tmp_out_path = tmp_in_path.replace('.csv', '_out.csv')

    try:
        cmd = [
            cchem_path, 'canonicalize',
            '-f', tmp_in_path,
            '-s', 'smiles',
            '-c', 'canonical',
            '-o', tmp_out_path,
            '-n', str(threads)
        ]

        start = time.perf_counter()
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.perf_counter() - start

        if result.returncode != 0:
            print(f"cchem error: {result.stderr}", file=sys.stderr)
            return BenchmarkResult(
                library='cchem',
                total_time=elapsed,
                molecules_processed=0,
                molecules_failed=len(smiles_list),
                throughput=0,
                avg_time_per_mol=0,
                threads=threads
            )

        success_count = 0
        if os.path.exists(tmp_out_path):
            with open(tmp_out_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row.get('canonical', '').strip():
                        success_count += 1

        failed = len(smiles_list) - success_count
        throughput = success_count / elapsed if elapsed > 0 else 0
        avg_time = (elapsed / success_count * 1e6) if success_count > 0 else 0

        return BenchmarkResult(
            library=f'cchem',
            total_time=elapsed,
            molecules_processed=success_count,
            molecules_failed=failed,
            throughput=throughput,
            avg_time_per_mol=avg_time,
            threads=threads
        )
    finally:
        os.unlink(tmp_in_path)
        if os.path.exists(tmp_out_path):
            os.unlink(tmp_out_path)


# =============================================================================
# RDKit benchmarks
# =============================================================================

def _rdkit_canonicalize_single(smi: str) -> int:
    """Canonicalize a single SMILES with RDKit. Returns 1 for success, 0 for failure."""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            canonical = Chem.MolToSmiles(mol, canonical=True)
            return 1 if canonical else 0
        return 0
    except Exception:
        return 0


def _rdkit_worker_init():
    """Initialize RDKit in worker process."""
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')


def benchmark_rdkit(smiles_list: list[str], threads: int = 1) -> Optional[BenchmarkResult]:
    """Benchmark RDKit canonicalization (single or multi-threaded)."""
    if not RDKIT_AVAILABLE:
        return None

    start = time.perf_counter()

    if threads == 1:
        # Single-threaded
        results = [_rdkit_canonicalize_single(smi) for smi in smiles_list]
    else:
        # Multi-threaded using multiprocessing
        with mp.Pool(threads, initializer=_rdkit_worker_init) as pool:
            results = pool.map(_rdkit_canonicalize_single, smiles_list, chunksize=1000)

    elapsed = time.perf_counter() - start

    success_count = sum(results)
    failed = len(smiles_list) - success_count
    throughput = success_count / elapsed if elapsed > 0 else 0
    avg_time = (elapsed / success_count * 1e6) if success_count > 0 else 0

    return BenchmarkResult(
        library='RDKit',
        total_time=elapsed,
        molecules_processed=success_count,
        molecules_failed=failed,
        throughput=throughput,
        avg_time_per_mol=avg_time,
        threads=threads
    )


# =============================================================================
# OpenBabel benchmarks
# =============================================================================

def _openbabel_canonicalize_single(smi: str) -> int:
    """Canonicalize a single SMILES with OpenBabel. Returns 1 for success, 0 for failure."""
    try:
        from openbabel import pybel
        mol = pybel.readstring('smi', smi)
        if mol:
            canonical = mol.write('can').strip()
            return 1 if canonical else 0
        return 0
    except Exception:
        return 0


def _openbabel_worker_init():
    """Initialize OpenBabel in worker process."""
    from openbabel import pybel
    pybel.ob.obErrorLog.SetOutputLevel(0)


def benchmark_openbabel(smiles_list: list[str], threads: int = 1) -> Optional[BenchmarkResult]:
    """Benchmark OpenBabel canonicalization (single or multi-threaded)."""
    if not OPENBABEL_AVAILABLE:
        return None

    start = time.perf_counter()

    if threads == 1:
        # Single-threaded
        results = [_openbabel_canonicalize_single(smi) for smi in smiles_list]
    else:
        # Multi-threaded using multiprocessing
        with mp.Pool(threads, initializer=_openbabel_worker_init) as pool:
            results = pool.map(_openbabel_canonicalize_single, smiles_list, chunksize=1000)

    elapsed = time.perf_counter() - start

    success_count = sum(results)
    failed = len(smiles_list) - success_count
    throughput = success_count / elapsed if elapsed > 0 else 0
    avg_time = (elapsed / success_count * 1e6) if success_count > 0 else 0

    return BenchmarkResult(
        library='OpenBabel',
        total_time=elapsed,
        molecules_processed=success_count,
        molecules_failed=failed,
        throughput=throughput,
        avg_time_per_mol=avg_time,
        threads=threads
    )


# =============================================================================
# CDK benchmarks
# =============================================================================

def benchmark_cdk(smiles_list: list[str], threads: int = 1) -> Optional[BenchmarkResult]:
    """Benchmark CDK canonicalization (single-threaded only due to JVM)."""
    if not CDK_AVAILABLE:
        return None

    # CDK via JPype doesn't support multiprocessing well
    if threads > 1:
        print("  Warning: CDK benchmark runs single-threaded (JVM limitation)")

    from org.openscience.cdk.smiles import SmilesParser, SmilesGenerator
    from org.openscience.cdk.silent import SilentChemObjectBuilder

    builder = SilentChemObjectBuilder.getInstance()
    parser = SmilesParser(builder)
    generator = SmilesGenerator.unique()

    success_count = 0
    failed = 0

    start = time.perf_counter()
    for smi in smiles_list:
        try:
            mol = parser.parseSmiles(smi)
            if mol:
                canonical = generator.create(mol)
                if canonical:
                    success_count += 1
                else:
                    failed += 1
            else:
                failed += 1
        except Exception:
            failed += 1
    elapsed = time.perf_counter() - start

    throughput = success_count / elapsed if elapsed > 0 else 0
    avg_time = (elapsed / success_count * 1e6) if success_count > 0 else 0

    return BenchmarkResult(
        library='CDK',
        total_time=elapsed,
        molecules_processed=success_count,
        molecules_failed=failed,
        throughput=throughput,
        avg_time_per_mol=avg_time,
        threads=1  # Always single-threaded
    )


# =============================================================================
# Utility functions
# =============================================================================

def find_cchem_binary() -> Optional[str]:
    """Find the cchem binary."""
    script_dir = Path(__file__).parent.parent
    candidates = [
        script_dir / 'build' / 'cchem',
        script_dir / 'cchem',
        Path('/usr/local/bin/cchem'),
        Path('/usr/bin/cchem'),
    ]

    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)

    import shutil
    path_cchem = shutil.which('cchem')
    if path_cchem:
        return path_cchem

    return None


def print_results(results: list[BenchmarkResult], baseline_name: str = None, title: str = "BENCHMARK RESULTS"):
    """Print benchmark results in a formatted table."""
    if not results:
        return

    print("\n" + "=" * 85)
    print(title)
    print("=" * 85)

    # Find baseline for speedup calculation
    baseline_result = None
    if baseline_name:
        for r in results:
            if baseline_name.lower() in r.library.lower() and '(n=1)' not in r.library:
                baseline_result = r
                break
    if not baseline_result:
        # Use single-threaded RDKit as default baseline
        for r in results:
            if 'rdkit' in r.library.lower():
                baseline_result = r
                break
        if not baseline_result and results:
            baseline_result = results[0]

    # Table header
    print(f"\n{'Library':<20} {'Threads':<8} {'Time (s)':<10} {'Throughput':<14} {'Avg (us)':<10} {'Failed':<8} {'Speedup':<10}")
    print("-" * 85)

    # Sort by throughput (fastest first)
    results_sorted = sorted(results, key=lambda x: x.throughput, reverse=True)

    for r in results_sorted:
        speedup = ""
        if baseline_result and baseline_result.throughput > 0:
            speedup_val = r.throughput / baseline_result.throughput
            speedup = f"{speedup_val:.2f}x"

        print(f"{r.library:<20} {r.threads:<8} {r.total_time:<10.3f} {r.throughput:<14,.0f} {r.avg_time_per_mol:<10.2f} {r.molecules_failed:<8} {speedup:<10}")

    print("-" * 85)
    if baseline_result:
        print(f"Baseline: {baseline_result.library} (n={baseline_result.threads})")
    print()


def run_warmup(smiles_sample: list[str], cchem_path: str, threads: int):
    """Run warmup iterations to stabilize timings."""
    print("Running warmup...")

    # Warmup RDKit
    if RDKIT_AVAILABLE:
        for smi in smiles_sample[:100]:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    Chem.MolToSmiles(mol, canonical=True)
            except:
                pass

    # Warmup cchem
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
        tmp.write('smiles\n')
        for smi in smiles_sample[:100]:
            tmp.write(f'{smi}\n')
        tmp_path = tmp.name

    try:
        subprocess.run([cchem_path, 'canonicalize', '-f', tmp_path, '-s', 'smiles',
                       '-c', 'canonical', '-o', '/dev/null', '-n', str(threads)],
                       capture_output=True)
    except:
        pass
    finally:
        os.unlink(tmp_path)


def run_benchmark_iteration(smiles_list: list[str], cchem_path: str, threads: int,
                           skip_rdkit: bool, skip_openbabel: bool, skip_cdk: bool,
                           verbose: bool) -> list[BenchmarkResult]:
    """Run one iteration of all benchmarks."""
    results = []

    # Single-threaded benchmarks
    print("  Single-threaded:")

    print("    cchem...", end=" ", flush=True)
    result = benchmark_cchem(smiles_list, cchem_path, threads=1)
    results.append(result)
    print(f"{result.throughput:,.0f} mol/s")

    if RDKIT_AVAILABLE and not skip_rdkit:
        print("    RDKit...", end=" ", flush=True)
        result = benchmark_rdkit(smiles_list, threads=1)
        if result:
            results.append(result)
            print(f"{result.throughput:,.0f} mol/s")

    if OPENBABEL_AVAILABLE and not skip_openbabel:
        print("    OpenBabel...", end=" ", flush=True)
        result = benchmark_openbabel(smiles_list, threads=1)
        if result:
            results.append(result)
            print(f"{result.throughput:,.0f} mol/s")

    if CDK_AVAILABLE and not skip_cdk:
        print("    CDK...", end=" ", flush=True)
        result = benchmark_cdk(smiles_list, threads=1)
        if result:
            results.append(result)
            print(f"{result.throughput:,.0f} mol/s")

    # Multi-threaded benchmarks (if threads > 1)
    if threads > 1:
        print(f"  Multi-threaded (n={threads}):")

        print("    cchem...", end=" ", flush=True)
        result = benchmark_cchem(smiles_list, cchem_path, threads=threads)
        results.append(result)
        print(f"{result.throughput:,.0f} mol/s")

        if RDKIT_AVAILABLE and not skip_rdkit:
            print("    RDKit...", end=" ", flush=True)
            result = benchmark_rdkit(smiles_list, threads=threads)
            if result:
                results.append(result)
                print(f"{result.throughput:,.0f} mol/s")

        if OPENBABEL_AVAILABLE and not skip_openbabel:
            print("    OpenBabel...", end=" ", flush=True)
            result = benchmark_openbabel(smiles_list, threads=threads)
            if result:
                results.append(result)
                print(f"{result.throughput:,.0f} mol/s")

    return results


def average_results(all_results: list[list[BenchmarkResult]]) -> list[BenchmarkResult]:
    """Average results across multiple iterations."""
    # Group by (library, threads)
    grouped = {}
    for iteration_results in all_results:
        for r in iteration_results:
            key = (r.library, r.threads)
            if key not in grouped:
                grouped[key] = []
            grouped[key].append(r)

    averaged = []
    for (library, threads), results in grouped.items():
        avg_result = BenchmarkResult(
            library=library,
            total_time=sum(r.total_time for r in results) / len(results),
            molecules_processed=results[0].molecules_processed,
            molecules_failed=results[0].molecules_failed,
            throughput=sum(r.throughput for r in results) / len(results),
            avg_time_per_mol=sum(r.avg_time_per_mol for r in results) / len(results),
            threads=threads
        )
        averaged.append(avg_result)

    return averaged


def main():
    parser = argparse.ArgumentParser(
        description='Benchmark SMILES canonicalization across cheminformatics libraries',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f molecules.csv -s smiles
  %(prog)s -f data.csv -s SMILES -n 8 --iterations 3
  %(prog)s -f large_dataset.csv -s smiles -n 4 --baseline rdkit

Available libraries:
  - cchem (this library, always available)
  - RDKit (pip install rdkit)
  - OpenBabel (pip install openbabel)
  - CDK (requires jpype and CDK_JAR environment variable)

Notes:
  - Single-threaded benchmarks always run for fair comparison
  - Multi-threaded benchmarks run when -n/--threads > 1
  - RDKit/OpenBabel use Python multiprocessing for parallelization
  - CDK is always single-threaded due to JVM limitations
"""
    )

    parser.add_argument('-f', '--file', required=True, help='Input CSV file with SMILES')
    parser.add_argument('-s', '--smiles-col', default='smiles', help='Column name containing SMILES (default: smiles)')
    parser.add_argument('-n', '--threads', type=int, default=mp.cpu_count(),
                       help=f'Number of threads for multi-threaded benchmarks (default: {mp.cpu_count()})')
    parser.add_argument('--iterations', type=int, default=3, help='Number of benchmark iterations (default: 3)')
    parser.add_argument('--warmup', type=int, default=1, help='Number of warmup iterations (default: 1)')
    parser.add_argument('--limit', type=int, default=None, help='Limit number of molecules to process')
    parser.add_argument('--baseline', default='rdkit', help='Library to use as baseline for speedup (default: rdkit)')
    parser.add_argument('--cchem-path', default=None, help='Path to cchem binary')
    parser.add_argument('--skip-rdkit', action='store_true', help='Skip RDKit benchmark')
    parser.add_argument('--skip-openbabel', action='store_true', help='Skip OpenBabel benchmark')
    parser.add_argument('--skip-cdk', action='store_true', help='Skip CDK benchmark')
    parser.add_argument('--single-only', action='store_true', help='Only run single-threaded benchmarks')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Find cchem binary
    cchem_path = args.cchem_path or find_cchem_binary()
    if not cchem_path:
        print("Error: Could not find cchem binary. Use --cchem-path to specify.", file=sys.stderr)
        sys.exit(1)

    print(f"Using cchem binary: {cchem_path}")

    # Print available libraries
    print("\nAvailable libraries:")
    print(f"  - cchem: Yes")
    print(f"  - RDKit: {'Yes' if RDKIT_AVAILABLE else 'No (pip install rdkit)'}")
    print(f"  - OpenBabel: {'Yes' if OPENBABEL_AVAILABLE else 'No (pip install openbabel)'}")
    print(f"  - CDK: {'Yes' if CDK_AVAILABLE else 'No (requires jpype + CDK_JAR)'}")

    # Load SMILES
    print(f"\nLoading SMILES from {args.file}...")
    try:
        smiles_list = load_smiles_from_csv(args.file, args.smiles_col, args.limit)
    except Exception as e:
        print(f"Error loading CSV: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(smiles_list):,} SMILES")

    if len(smiles_list) == 0:
        print("Error: No SMILES loaded", file=sys.stderr)
        sys.exit(1)

    threads = 1 if args.single_only else args.threads

    # Warmup
    for _ in range(args.warmup):
        run_warmup(smiles_list, cchem_path, threads)

    # Run benchmarks
    all_results = []

    for iteration in range(args.iterations):
        print(f"\n--- Iteration {iteration + 1}/{args.iterations} ---")
        iteration_results = run_benchmark_iteration(
            smiles_list, cchem_path, threads,
            args.skip_rdkit, args.skip_openbabel, args.skip_cdk,
            args.verbose
        )
        all_results.append(iteration_results)

    # Average and print results
    if args.iterations > 1:
        print("\nAveraging results across iterations...")
        final_results = average_results(all_results)
    else:
        final_results = all_results[0]

    # Separate single and multi-threaded for clearer output
    single_results = [r for r in final_results if r.threads == 1]
    multi_results = [r for r in final_results if r.threads > 1]

    if single_results:
        print_results(single_results, args.baseline, "SINGLE-THREADED RESULTS")

    if multi_results:
        print_results(multi_results, args.baseline, f"MULTI-THREADED RESULTS (n={threads})")

    # Print scaling summary if we have both
    if single_results and multi_results:
        print("=" * 85)
        print("SCALING EFFICIENCY")
        print("=" * 85)
        print(f"\n{'Library':<20} {'1-thread':<14} {f'{threads}-thread':<14} {'Scaling':<10}")
        print("-" * 60)

        for sr in sorted(single_results, key=lambda x: x.library):
            mr = next((r for r in multi_results if r.library == sr.library), None)
            if mr:
                scaling = mr.throughput / sr.throughput if sr.throughput > 0 else 0
                efficiency = (scaling / threads) * 100
                print(f"{sr.library:<20} {sr.throughput:<14,.0f} {mr.throughput:<14,.0f} {scaling:.2f}x ({efficiency:.0f}%)")

        print("-" * 60)
        print(f"Perfect scaling would be {threads:.2f}x (100%)")
        print()

    # Print summary
    print(f"Dataset: {args.file}")
    print(f"Molecules: {len(smiles_list):,}")
    print(f"Iterations: {args.iterations}")
    print(f"Threads: {threads}")


if __name__ == '__main__':
    main()
