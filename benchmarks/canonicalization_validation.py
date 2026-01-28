#!/usr/bin/env python3
"""
Validation script for cchem SMILES canonicalization.

Compares cchem canonicalization against RDKit (gold standard) to verify:
1. Canonical equivalence - Do both produce equivalent molecules?
2. Round-trip consistency - Is mol → smiles → mol → smiles stable?
3. Stereochemistry preservation - Is chirality and E/Z geometry preserved?
4. Failure rate - How many molecules fail to parse/canonicalize?

Usage:
    python canonicalization_validation.py -f molecules.csv -s smiles
    python canonicalization_validation.py -f molecules.csv -s smiles --verbose --limit 1000
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Error: RDKit is required for validation. Install with: pip install rdkit", file=sys.stderr)
    sys.exit(1)


@dataclass
class ValidationResult:
    """Results of validation."""
    total_molecules: int = 0

    # Parsing failures
    cchem_parse_failures: int = 0
    rdkit_parse_failures: int = 0
    both_parse_failures: int = 0

    # Canonicalization results
    canonical_matches: int = 0  # Exact string match
    equivalent_molecules: int = 0  # Different strings but equivalent molecules
    non_equivalent: int = 0  # Different molecules (ERROR)

    # Round-trip consistency
    cchem_roundtrip_stable: int = 0
    cchem_roundtrip_unstable: int = 0
    rdkit_roundtrip_stable: int = 0
    rdkit_roundtrip_unstable: int = 0

    # Stereochemistry
    stereo_preserved: int = 0
    stereo_lost: int = 0
    stereo_added: int = 0  # cchem added stereo not in original
    stereo_molecules: int = 0  # molecules with stereochemistry

    # Detailed discrepancies
    discrepancies: list = field(default_factory=list)
    stereo_issues: list = field(default_factory=list)


def find_cchem_binary() -> Optional[str]:
    """Find the cchem binary."""
    script_dir = Path(__file__).parent.parent
    candidates = [
        script_dir / 'build' / 'cchem',
        script_dir / 'cchem',
        Path('/usr/local/bin/cchem'),
    ]
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)
    import shutil
    return shutil.which('cchem')


def canonicalize_with_cchem(smiles: str, cchem_path: str) -> Optional[str]:
    """Canonicalize a single SMILES with cchem."""
    try:
        result = subprocess.run(
            [cchem_path, 'canonicalize', '-S', smiles],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return None
    except Exception:
        return None


def canonicalize_with_cchem_batch(smiles_list: list[str], cchem_path: str) -> dict[str, Optional[str]]:
    """Canonicalize multiple SMILES with cchem using batch mode."""
    results = {}

    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_in:
        tmp_in.write('smiles\n')
        for smi in smiles_list:
            tmp_in.write(f'{smi}\n')
        tmp_in_path = tmp_in.name

    tmp_out_path = tmp_in_path.replace('.csv', '_out.csv')

    try:
        subprocess.run(
            [cchem_path, 'canonicalize', '-f', tmp_in_path, '-s', 'smiles',
             '-c', 'canonical', '-o', tmp_out_path],
            capture_output=True, timeout=300
        )

        if os.path.exists(tmp_out_path):
            with open(tmp_out_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    orig = row.get('smiles', '').strip()
                    canon = row.get('canonical', '').strip()
                    results[orig] = canon if canon else None
    except Exception as e:
        print(f"Batch canonicalization error: {e}", file=sys.stderr)
    finally:
        os.unlink(tmp_in_path)
        if os.path.exists(tmp_out_path):
            os.unlink(tmp_out_path)

    return results


def canonicalize_with_rdkit(smiles: str) -> Optional[str]:
    """Canonicalize a single SMILES with RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True)
        return None
    except Exception:
        return None


def molecules_equivalent(smiles1: str, smiles2: str) -> bool:
    """Check if two SMILES represent equivalent molecules using InChI."""
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 is None or mol2 is None:
            return False

        # Use canonical SMILES comparison (RDKit as arbiter)
        canon1 = Chem.MolToSmiles(mol1, canonical=True)
        canon2 = Chem.MolToSmiles(mol2, canonical=True)
        return canon1 == canon2
    except Exception:
        return False


def has_stereochemistry(smiles: str) -> bool:
    """Check if SMILES contains stereochemistry markers."""
    return '@' in smiles or '/' in smiles or '\\' in smiles


def count_stereocenters(smiles: str) -> tuple[int, int]:
    """Count chiral centers and double bond stereo in a molecule."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (0, 0)

        # Count chiral centers
        chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))

        # Count E/Z double bonds
        ez_bonds = 0
        for bond in mol.GetBonds():
            if bond.GetStereo() in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                ez_bonds += 1

        return (chiral_centers, ez_bonds)
    except Exception:
        return (0, 0)


def load_smiles_from_csv(filepath: str, smiles_col: str, limit: Optional[int] = None) -> list[str]:
    """Load SMILES strings from a CSV file."""
    smiles_list = []
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if smiles_col not in reader.fieldnames:
            available = ', '.join(reader.fieldnames)
            raise ValueError(f"Column '{smiles_col}' not found. Available: {available}")

        for i, row in enumerate(reader):
            if limit and i >= limit:
                break
            smi = row[smiles_col].strip()
            if smi:
                smiles_list.append(smi)

    return smiles_list


def validate_canonicalization(smiles_list: list[str], cchem_path: str,
                             verbose: bool = False) -> ValidationResult:
    """Run full validation comparing cchem vs RDKit."""
    result = ValidationResult()
    result.total_molecules = len(smiles_list)

    print(f"Validating {len(smiles_list):,} molecules...")
    print("Step 1/4: Canonicalizing with cchem (batch)...")

    # Batch canonicalize with cchem
    cchem_results = canonicalize_with_cchem_batch(smiles_list, cchem_path)

    print("Step 2/4: Canonicalizing with RDKit...")
    rdkit_results = {}
    for i, smi in enumerate(smiles_list):
        rdkit_results[smi] = canonicalize_with_rdkit(smi)
        if (i + 1) % 5000 == 0:
            print(f"  {i+1:,}/{len(smiles_list):,}")

    print("Step 3/4: Comparing results...")

    for smi in smiles_list:
        cchem_canon = cchem_results.get(smi)
        rdkit_canon = rdkit_results.get(smi)

        # Check parsing failures
        if cchem_canon is None and rdkit_canon is None:
            result.both_parse_failures += 1
            continue
        elif cchem_canon is None:
            result.cchem_parse_failures += 1
            if verbose:
                result.discrepancies.append({
                    'type': 'cchem_parse_failure',
                    'original': smi,
                    'rdkit': rdkit_canon
                })
            continue
        elif rdkit_canon is None:
            result.rdkit_parse_failures += 1
            continue

        # Compare canonical forms
        if cchem_canon == rdkit_canon:
            result.canonical_matches += 1
        elif molecules_equivalent(cchem_canon, rdkit_canon):
            result.equivalent_molecules += 1
            if verbose:
                result.discrepancies.append({
                    'type': 'different_but_equivalent',
                    'original': smi,
                    'cchem': cchem_canon,
                    'rdkit': rdkit_canon
                })
        else:
            result.non_equivalent += 1
            result.discrepancies.append({
                'type': 'NON_EQUIVALENT',
                'original': smi,
                'cchem': cchem_canon,
                'rdkit': rdkit_canon
            })

        # Check stereochemistry
        if has_stereochemistry(smi):
            result.stereo_molecules += 1
            orig_stereo = count_stereocenters(smi)
            cchem_stereo = count_stereocenters(cchem_canon)

            if orig_stereo == cchem_stereo:
                result.stereo_preserved += 1
            elif sum(cchem_stereo) < sum(orig_stereo):
                result.stereo_lost += 1
                result.stereo_issues.append({
                    'type': 'stereo_lost',
                    'original': smi,
                    'cchem': cchem_canon,
                    'orig_stereo': orig_stereo,
                    'cchem_stereo': cchem_stereo
                })
            else:
                result.stereo_added += 1

    print("Step 4/4: Testing round-trip consistency...")

    # Test round-trip on a sample
    sample_size = min(1000, len(smiles_list))
    sample = smiles_list[:sample_size]

    for smi in sample:
        # cchem round-trip
        cchem1 = cchem_results.get(smi)
        if cchem1:
            cchem2 = canonicalize_with_cchem(cchem1, cchem_path)
            if cchem2 == cchem1:
                result.cchem_roundtrip_stable += 1
            else:
                result.cchem_roundtrip_unstable += 1
                if verbose and cchem2:
                    result.discrepancies.append({
                        'type': 'cchem_roundtrip_unstable',
                        'original': smi,
                        'first': cchem1,
                        'second': cchem2
                    })

        # RDKit round-trip
        rdkit1 = rdkit_results.get(smi)
        if rdkit1:
            rdkit2 = canonicalize_with_rdkit(rdkit1)
            if rdkit2 == rdkit1:
                result.rdkit_roundtrip_stable += 1
            else:
                result.rdkit_roundtrip_unstable += 1

    return result


def print_results(result: ValidationResult, verbose: bool = False):
    """Print validation results."""
    print("\n" + "=" * 80)
    print("CANONICALIZATION VALIDATION RESULTS")
    print("=" * 80)

    print(f"\nTotal molecules: {result.total_molecules:,}")

    # Parsing
    print("\n--- Parsing ---")
    print(f"  cchem failures:  {result.cchem_parse_failures:,} ({100*result.cchem_parse_failures/result.total_molecules:.2f}%)")
    print(f"  RDKit failures:  {result.rdkit_parse_failures:,} ({100*result.rdkit_parse_failures/result.total_molecules:.2f}%)")
    print(f"  Both failed:     {result.both_parse_failures:,}")

    # Canonicalization comparison
    comparable = result.total_molecules - result.cchem_parse_failures - result.rdkit_parse_failures - result.both_parse_failures
    print("\n--- Canonical Comparison ---")
    print(f"  Comparable molecules: {comparable:,}")
    print(f"  Exact match:          {result.canonical_matches:,} ({100*result.canonical_matches/comparable:.2f}%)")
    print(f"  Equivalent (diff str):{result.equivalent_molecules:,} ({100*result.equivalent_molecules/comparable:.2f}%)")
    print(f"  NON-EQUIVALENT:       {result.non_equivalent:,} ({100*result.non_equivalent/comparable:.2f}%)")

    total_correct = result.canonical_matches + result.equivalent_molecules
    print(f"\n  ** Correctness rate: {100*total_correct/comparable:.4f}% **")

    # Stereochemistry
    if result.stereo_molecules > 0:
        print("\n--- Stereochemistry ---")
        print(f"  Molecules with stereo: {result.stereo_molecules:,}")
        print(f"  Preserved:             {result.stereo_preserved:,} ({100*result.stereo_preserved/result.stereo_molecules:.2f}%)")
        print(f"  Lost:                  {result.stereo_lost:,} ({100*result.stereo_lost/result.stereo_molecules:.2f}%)")
        print(f"  Added:                 {result.stereo_added:,}")

    # Round-trip
    roundtrip_total = result.cchem_roundtrip_stable + result.cchem_roundtrip_unstable
    if roundtrip_total > 0:
        print("\n--- Round-trip Consistency (sample) ---")
        print(f"  cchem stable:   {result.cchem_roundtrip_stable:,}/{roundtrip_total} ({100*result.cchem_roundtrip_stable/roundtrip_total:.2f}%)")
        print(f"  RDKit stable:   {result.rdkit_roundtrip_stable:,}/{roundtrip_total} ({100*result.rdkit_roundtrip_stable/roundtrip_total:.2f}%)")

    # Show discrepancies
    if result.non_equivalent > 0:
        print("\n" + "=" * 80)
        print("NON-EQUIVALENT MOLECULES (ERRORS)")
        print("=" * 80)
        non_equiv = [d for d in result.discrepancies if d['type'] == 'NON_EQUIVALENT']
        for i, d in enumerate(non_equiv[:20]):  # Show first 20
            print(f"\n{i+1}. Original: {d['original']}")
            print(f"   cchem:    {d['cchem']}")
            print(f"   RDKit:    {d['rdkit']}")
        if len(non_equiv) > 20:
            print(f"\n... and {len(non_equiv) - 20} more")

    if verbose and result.stereo_issues:
        print("\n" + "=" * 80)
        print("STEREOCHEMISTRY ISSUES")
        print("=" * 80)
        for i, d in enumerate(result.stereo_issues[:10]):
            print(f"\n{i+1}. Original: {d['original']} (stereo: {d['orig_stereo']})")
            print(f"   cchem:    {d['cchem']} (stereo: {d['cchem_stereo']})")

    print("\n" + "=" * 80)

    # Summary verdict
    if result.non_equivalent == 0 and result.cchem_parse_failures < result.total_molecules * 0.01:
        print("VERDICT: PASS - cchem canonicalization is correct")
    elif result.non_equivalent < comparable * 0.001:
        print(f"VERDICT: MOSTLY PASS - {result.non_equivalent} errors ({100*result.non_equivalent/comparable:.4f}%)")
    else:
        print(f"VERDICT: FAIL - {result.non_equivalent} non-equivalent results")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Validate cchem SMILES canonicalization against RDKit'
    )
    parser.add_argument('-f', '--file', required=True, help='Input CSV file')
    parser.add_argument('-s', '--smiles-col', default='smiles', help='SMILES column name')
    parser.add_argument('--limit', type=int, help='Limit number of molecules')
    parser.add_argument('--cchem-path', help='Path to cchem binary')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show detailed discrepancies')
    parser.add_argument('-o', '--output', help='Output discrepancies to CSV file')

    args = parser.parse_args()

    # Find cchem
    cchem_path = args.cchem_path or find_cchem_binary()
    if not cchem_path:
        print("Error: cchem binary not found", file=sys.stderr)
        sys.exit(1)
    print(f"Using cchem: {cchem_path}")

    # Load SMILES
    print(f"Loading SMILES from {args.file}...")
    smiles_list = load_smiles_from_csv(args.file, args.smiles_col, args.limit)
    print(f"Loaded {len(smiles_list):,} SMILES")

    # Run validation
    result = validate_canonicalization(smiles_list, cchem_path, args.verbose)

    # Print results
    print_results(result, args.verbose)

    # Save discrepancies if requested
    if args.output and result.discrepancies:
        with open(args.output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['type', 'original', 'cchem', 'rdkit', 'first', 'second'])
            writer.writeheader()
            for d in result.discrepancies:
                writer.writerow({k: d.get(k, '') for k in writer.fieldnames})
        print(f"\nDiscrepancies saved to {args.output}")

    # Exit code
    sys.exit(0 if result.non_equivalent == 0 else 1)


if __name__ == '__main__':
    main()
