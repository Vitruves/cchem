"""
Command-line interface for pycchem.

Provides similar functionality to the cchem C CLI.
"""

import argparse
import csv
import sys
from typing import List, Optional


def cmd_canonicalize(args: argparse.Namespace) -> int:
    """Canonicalize SMILES strings."""
    import pycchem

    smiles_list = []

    if args.smiles:
        smiles_list.append(args.smiles)

    if args.file:
        try:
            with open(args.file, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        smiles_list.append(line)
        except IOError as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            return 1

    if not smiles_list:
        print("Error: No SMILES provided. Use -S or -f option.", file=sys.stderr)
        return 1

    errors = 0
    for smi in smiles_list:
        try:
            result = pycchem.canonicalize(
                smi,
                isomeric=not args.no_stereo,
                kekulize=args.kekulize,
            )
            print(result)
        except pycchem.CChemError as e:
            if args.verbose:
                print(f"Error: {smi} -> {e}", file=sys.stderr)
            errors += 1

    return 1 if errors == len(smiles_list) else 0


def cmd_validate(args: argparse.Namespace) -> int:
    """Validate SMILES strings."""
    import pycchem

    smiles_list = []

    if args.smiles:
        smiles_list.append(args.smiles)

    if args.file:
        try:
            with open(args.file, "r") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        smiles_list.append(line)
        except IOError as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            return 1

    if not smiles_list:
        print("Error: No SMILES provided. Use -S or -f option.", file=sys.stderr)
        return 1

    all_valid = True
    for smi in smiles_list:
        valid = pycchem.validate(smi)
        if args.verbose or not valid:
            status = "valid" if valid else "invalid"
            print(f"{smi}: {status}")
        if not valid:
            all_valid = False

    return 0 if all_valid else 1


def cmd_compute(args: argparse.Namespace) -> int:
    """Compute molecular descriptors."""
    import pycchem

    if args.list:
        # List available descriptors (no SMILES needed)
        descriptors = pycchem.list_descriptors()
        for name in descriptors:
            print(name)
        print(f"\nTotal: {len(descriptors)} descriptors", file=sys.stderr)
        return 0

    if not args.smiles:
        print("Error: SMILES required. Use -S option.", file=sys.stderr)
        return 1

    try:
        mol = pycchem.Molecule(args.smiles)
    except pycchem.CChemError as e:
        print(f"Error parsing SMILES: {e}", file=sys.stderr)
        return 1

    if args.descriptor:
        # Compute specific descriptor(s)
        descriptors = [d.strip() for d in args.descriptor.split(",")]
        try:
            results = mol.descriptors(descriptors)
            for name in descriptors:
                value = results.get(name)
                if value is not None:
                    if isinstance(value, float):
                        print(f"{name}: {value:.6f}")
                    else:
                        print(f"{name}: {value}")
                else:
                    print(f"{name}: N/A", file=sys.stderr)
        except pycchem.CChemError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
    elif args.all:
        # Compute all descriptors
        try:
            results = mol.descriptors()
            if args.output:
                # Write to CSV
                with open(args.output, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["descriptor", "value"])
                    for name in sorted(results.keys()):
                        writer.writerow([name, results[name]])
                print(f"Wrote {len(results)} descriptors to {args.output}")
            else:
                # Print to stdout
                for name in sorted(results.keys()):
                    value = results[name]
                    if isinstance(value, float):
                        print(f"{name}: {value:.6f}")
                    else:
                        print(f"{name}: {value}")
        except pycchem.CChemError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
    else:
        print("Error: Specify -d DESCRIPTOR, --all, or --list", file=sys.stderr)
        return 1

    return 0


def cmd_sanitize(args: argparse.Namespace) -> int:
    """Sanitize SMILES strings."""
    import pycchem

    if not args.smiles:
        print("Error: SMILES required. Use -S option.", file=sys.stderr)
        return 1

    operations = []
    if args.unsalt:
        operations.append("unsalt")
    if args.neutralize:
        operations.append("neutralize")
    if args.aromatize:
        operations.append("aromatize")
    if args.kekulize:
        operations.append("kekulize")
    if args.normalize:
        operations.append("normalize")

    if not operations:
        operations = "complete"

    try:
        result = pycchem.sanitize(args.smiles, operations=operations)
        print(result)
        return 0
    except pycchem.CChemError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def cmd_compare(args: argparse.Namespace) -> int:
    """Compare two SMILES strings."""
    import pycchem

    if not args.smiles1 or not args.smiles2:
        print("Error: Two SMILES required.", file=sys.stderr)
        return 1

    try:
        equivalent = pycchem.are_equivalent(args.smiles1, args.smiles2)
        if equivalent:
            print("equivalent")
            return 0
        else:
            print("different")
            return 1
    except pycchem.CChemError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def cmd_depict(args: argparse.Namespace) -> int:
    """Depict molecules (requires native cchem binary)."""
    print("Error: The 'depict' command requires the native cchem binary.", file=sys.stderr)
    print("", file=sys.stderr)
    print("To use molecular depiction, install cchem from source:", file=sys.stderr)
    print("  git clone https://github.com/Vitruves/cchem.git", file=sys.stderr)
    print("  cd cchem && mkdir build && cd build", file=sys.stderr)
    print("  cmake .. && make && sudo cmake --install .", file=sys.stderr)
    print("", file=sys.stderr)
    print("Then use: cchem depict -S 'c1ccccc1' -o molecule.png", file=sys.stderr)
    return 1


def cmd_version(args: argparse.Namespace) -> int:
    """Show version information."""
    import pycchem

    print(f"pycchem {pycchem.__version__}")
    print(f"cchem library {pycchem.version()}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        prog="cchem",
        description="High-performance cheminformatics toolkit",
    )
    parser.add_argument(
        "--version", "-V",
        action="store_true",
        help="Show version information",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # canonicalize command
    p_canon = subparsers.add_parser(
        "canonicalize",
        aliases=["canon", "c"],
        help="Canonicalize SMILES strings",
    )
    p_canon.add_argument("-S", "--smiles", help="SMILES string to canonicalize")
    p_canon.add_argument("-f", "--file", help="File containing SMILES (one per line)")
    p_canon.add_argument("--no-stereo", action="store_true", help="Ignore stereochemistry")
    p_canon.add_argument("--kekulize", "-k", action="store_true", help="Output Kekule form")
    p_canon.add_argument("-v", "--verbose", action="store_true", help="Show errors")
    p_canon.set_defaults(func=cmd_canonicalize)

    # validate command
    p_validate = subparsers.add_parser(
        "validate",
        aliases=["v"],
        help="Validate SMILES syntax",
    )
    p_validate.add_argument("-S", "--smiles", help="SMILES string to validate")
    p_validate.add_argument("-f", "--file", help="File containing SMILES (one per line)")
    p_validate.add_argument("-v", "--verbose", action="store_true", help="Show all results")
    p_validate.set_defaults(func=cmd_validate)

    # compute command
    p_compute = subparsers.add_parser(
        "compute",
        aliases=["desc", "d"],
        help="Compute molecular descriptors",
    )
    p_compute.add_argument("-S", "--smiles", help="SMILES string")
    p_compute.add_argument("-d", "--descriptor", help="Descriptor name(s), comma-separated")
    p_compute.add_argument("--all", "-a", action="store_true", help="Compute all descriptors")
    p_compute.add_argument("--list", "-l", action="store_true", help="List available descriptors")
    p_compute.add_argument("-o", "--output", help="Output CSV file (with --all)")
    p_compute.set_defaults(func=cmd_compute)

    # sanitize command
    p_sanitize = subparsers.add_parser(
        "sanitize",
        aliases=["s"],
        help="Sanitize SMILES strings",
    )
    p_sanitize.add_argument("-S", "--smiles", required=True, help="SMILES string")
    p_sanitize.add_argument("--unsalt", action="store_true", help="Remove salts")
    p_sanitize.add_argument("--neutralize", action="store_true", help="Neutralize charges")
    p_sanitize.add_argument("--aromatize", action="store_true", help="Aromatize rings")
    p_sanitize.add_argument("--kekulize", action="store_true", help="Kekulize rings")
    p_sanitize.add_argument("--normalize", action="store_true", help="Normalize functional groups")
    p_sanitize.set_defaults(func=cmd_sanitize)

    # compare command
    p_compare = subparsers.add_parser(
        "compare",
        aliases=["cmp"],
        help="Compare two SMILES strings",
    )
    p_compare.add_argument("smiles1", help="First SMILES string")
    p_compare.add_argument("smiles2", help="Second SMILES string")
    p_compare.set_defaults(func=cmd_compare)

    # depict command (requires native binary)
    p_depict = subparsers.add_parser(
        "depict",
        help="Depict molecules (requires native cchem binary)",
    )
    p_depict.add_argument("-S", "--smiles", help="SMILES string to depict")
    p_depict.add_argument("-o", "--output", help="Output file (PNG, JPG, SVG)")
    p_depict.add_argument("-m", "--mode", choices=["2d", "3d"], default="2d", help="Depiction mode")
    p_depict.set_defaults(func=cmd_depict)

    args = parser.parse_args(argv)

    if args.version:
        return cmd_version(args)

    if not args.command:
        parser.print_help()
        return 0

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
