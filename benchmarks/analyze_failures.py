#!/usr/bin/env python3
"""
Analyze canonicalization failures to identify patterns.
"""

import csv
import subprocess
import tempfile
import os
from collections import defaultdict, Counter
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def find_cchem():
    script_dir = Path(__file__).parent.parent
    for p in [script_dir / 'build' / 'cchem', script_dir / 'cchem']:
        if p.exists():
            return str(p)
    return None


def canonicalize_cchem_batch(smiles_list, cchem_path):
    results = {}
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write('smiles\n')
        for smi in smiles_list:
            f.write(f'{smi}\n')
        tmp_in = f.name
    tmp_out = tmp_in.replace('.csv', '_out.csv')

    try:
        subprocess.run([cchem_path, 'canonicalize', '-f', tmp_in, '-s', 'smiles',
                       '-c', 'canonical', '-o', tmp_out], capture_output=True, timeout=300)
        if os.path.exists(tmp_out):
            with open(tmp_out) as f:
                for row in csv.DictReader(f):
                    results[row['smiles']] = row.get('canonical', '').strip() or None
    finally:
        os.unlink(tmp_in)
        if os.path.exists(tmp_out):
            os.unlink(tmp_out)
    return results


def analyze_molecule_difference(orig, cchem_smi, rdkit_smi):
    """Analyze why two SMILES are different."""
    issues = []

    # Parse molecules
    mol_orig = Chem.MolFromSmiles(orig)
    mol_cchem = Chem.MolFromSmiles(cchem_smi)
    mol_rdkit = Chem.MolFromSmiles(rdkit_smi)

    if mol_cchem is None:
        return ['cchem_invalid_smiles']

    # Compare atom counts
    if mol_cchem.GetNumAtoms() != mol_rdkit.GetNumAtoms():
        issues.append(f'atom_count_diff:{mol_cchem.GetNumAtoms()}vs{mol_rdkit.GetNumAtoms()}')

    # Compare bond counts
    if mol_cchem.GetNumBonds() != mol_rdkit.GetNumBonds():
        issues.append(f'bond_count_diff:{mol_cchem.GetNumBonds()}vs{mol_rdkit.GetNumBonds()}')

    # Compare molecular formula
    from rdkit.Chem import rdMolDescriptors
    formula_cchem = rdMolDescriptors.CalcMolFormula(mol_cchem)
    formula_rdkit = rdMolDescriptors.CalcMolFormula(mol_rdkit)
    if formula_cchem != formula_rdkit:
        issues.append(f'formula_diff:{formula_cchem}vs{formula_rdkit}')

    # Check for aromatic nitrogen differences
    def count_aromatic_nh(mol):
        count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7 and atom.GetIsAromatic():
                if atom.GetTotalNumHs() > 0:
                    count += 1
        return count

    nh_cchem = count_aromatic_nh(mol_cchem)
    nh_rdkit = count_aromatic_nh(mol_rdkit)
    if nh_cchem != nh_rdkit:
        issues.append(f'aromatic_nh_diff:{nh_cchem}vs{nh_rdkit}')

    # Check ring count
    from rdkit.Chem import rdMolDescriptors
    ring_cchem = rdMolDescriptors.CalcNumRings(mol_cchem)
    ring_rdkit = rdMolDescriptors.CalcNumRings(mol_rdkit)
    if ring_cchem != ring_rdkit:
        issues.append(f'ring_count_diff:{ring_cchem}vs{ring_rdkit}')

    # Check aromatic ring count
    aromatic_cchem = rdMolDescriptors.CalcNumAromaticRings(mol_cchem)
    aromatic_rdkit = rdMolDescriptors.CalcNumAromaticRings(mol_rdkit)
    if aromatic_cchem != aromatic_rdkit:
        issues.append(f'aromatic_ring_diff:{aromatic_cchem}vs{aromatic_rdkit}')

    # Check for specific heteroaromatic patterns
    def has_pattern(mol, smarts):
        pat = Chem.MolFromSmarts(smarts)
        return mol.HasSubstructMatch(pat) if pat else False

    # Imidazole pattern [nH]1ccnc1
    if has_pattern(mol_rdkit, '[nH]1ccnc1') != has_pattern(mol_cchem, '[nH]1ccnc1'):
        issues.append('imidazole_diff')

    # Pyrazole pattern [nH]1nccc1
    if has_pattern(mol_rdkit, '[nH]1nccc1') != has_pattern(mol_cchem, '[nH]1nccc1'):
        issues.append('pyrazole_diff')

    # Triazole patterns
    if has_pattern(mol_rdkit, '[nH]1nncn1') != has_pattern(mol_cchem, '[nH]1nncn1'):
        issues.append('triazole_diff')

    # Benzimidazole
    if has_pattern(mol_rdkit, '[nH]1cnc2ccccc12') != has_pattern(mol_cchem, '[nH]1cnc2ccccc12'):
        issues.append('benzimidazole_diff')

    # Check InChI
    from rdkit.Chem.inchi import MolToInchi
    try:
        inchi_cchem = MolToInchi(mol_cchem)
        inchi_rdkit = MolToInchi(mol_rdkit)
        if inchi_cchem != inchi_rdkit:
            issues.append('inchi_diff')
            # Check if only hydrogen layer differs
            if inchi_cchem and inchi_rdkit:
                # InChI format: InChI=1S/formula/c.../h...
                parts_cchem = inchi_cchem.split('/')
                parts_rdkit = inchi_rdkit.split('/')
                if len(parts_cchem) >= 2 and len(parts_rdkit) >= 2:
                    if parts_cchem[1] == parts_rdkit[1]:  # Same formula
                        # Find h layer
                        h_cchem = [p for p in parts_cchem if p.startswith('h')]
                        h_rdkit = [p for p in parts_rdkit if p.startswith('h')]
                        if h_cchem != h_rdkit:
                            issues.append('inchi_h_layer_diff')
    except:
        pass

    if not issues:
        issues.append('unknown_diff')

    return issues


def categorize_failures(smiles_list, cchem_path, limit=1000):
    """Categorize all failures."""
    print(f"Analyzing {min(len(smiles_list), limit)} molecules...")

    # Batch canonicalize
    sample = smiles_list[:limit]
    cchem_results = canonicalize_cchem_batch(sample, cchem_path)

    categories = defaultdict(list)
    issue_counts = Counter()

    for i, smi in enumerate(sample):
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(sample)}")

        cchem_canon = cchem_results.get(smi)
        rdkit_canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True) if Chem.MolFromSmiles(smi) else None

        if not cchem_canon or not rdkit_canon:
            continue

        if cchem_canon == rdkit_canon:
            continue

        # Check equivalence
        mol_cchem = Chem.MolFromSmiles(cchem_canon)
        mol_rdkit = Chem.MolFromSmiles(rdkit_canon)

        if mol_cchem is None:
            categories['invalid_cchem_smiles'].append((smi, cchem_canon, rdkit_canon))
            issue_counts['invalid_cchem_smiles'] += 1
            continue

        # Re-canonicalize with RDKit to check equivalence
        cchem_recanon = Chem.MolToSmiles(mol_cchem, canonical=True)

        if cchem_recanon == rdkit_canon:
            categories['equivalent_different_string'].append((smi, cchem_canon, rdkit_canon))
            continue

        # Non-equivalent - analyze why
        issues = analyze_molecule_difference(smi, cchem_canon, rdkit_canon)
        for issue in issues:
            issue_counts[issue] += 1
            categories[issue].append((smi, cchem_canon, rdkit_canon))

    return categories, issue_counts


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True)
    parser.add_argument('-s', '--smiles-col', default='smiles')
    parser.add_argument('--limit', type=int, default=5000)
    args = parser.parse_args()

    cchem_path = find_cchem()
    if not cchem_path:
        print("cchem not found")
        return

    # Load SMILES
    smiles_list = []
    with open(args.file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            smi = row.get(args.smiles_col, '').strip()
            if smi:
                smiles_list.append(smi)

    print(f"Loaded {len(smiles_list)} SMILES")

    categories, issue_counts = categorize_failures(smiles_list, cchem_path, args.limit)

    print("\n" + "=" * 80)
    print("FAILURE ANALYSIS")
    print("=" * 80)

    print("\nIssue counts (sorted by frequency):")
    for issue, count in issue_counts.most_common():
        print(f"  {issue}: {count}")

    print("\n" + "=" * 80)
    print("EXAMPLES BY CATEGORY")
    print("=" * 80)

    for category in ['inchi_h_layer_diff', 'aromatic_nh_diff', 'imidazole_diff',
                     'pyrazole_diff', 'triazole_diff', 'benzimidazole_diff',
                     'formula_diff', 'invalid_cchem_smiles', 'unknown_diff']:
        examples = categories.get(category, [])
        if examples:
            print(f"\n--- {category} ({len(examples)} cases) ---")
            for orig, cchem, rdkit in examples[:5]:
                print(f"  Original: {orig}")
                print(f"  cchem:    {cchem}")
                print(f"  RDKit:    {rdkit}")
                print()


if __name__ == '__main__':
    main()
