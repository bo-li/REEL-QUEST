#!/usr/bin/env python3
"""
Query Ligands CLI Tool
Command-line interface for querying the REEL-QUEST ligand database
"""

import os
import sys
import argparse
import json

# Add parent directory to path to import rel_quest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rel_quest import LigandDatabase, LigandQuery
# Vocabularies come from rel_quest.schema, not literals: this file's --trend
# choices were left at five values when the classifier gained Peaked and
# U_Shaped, so 18 ligands were unreachable from the command line.
from rel_quest.schema import TREND_CLASSIFICATIONS, canonical_mechanisms


def format_json(data: dict, indent: int = 2) -> str:
    """Format dictionary as pretty JSON"""
    return json.dumps(data, indent=indent)


def format_table(ligand: LigandQuery, format_type: str) -> str:
    """
    Format ligand data as ASCII table

    Args:
        ligand: LigandQuery object
        format_type: One of 'raw', 'trend', 'stats', 'full'

    Returns:
        Formatted string
    """
    lines = []
    lines.append("=" * 70)
    lines.append(f"Ligand: {ligand.name}")
    lines.append("=" * 70)

    if format_type in ['full', 'trend']:
        lines.append(f"Mechanism: {ligand.mechanism}")
        if ligand.structural_class:   # legacy field; the public release omits it
            lines.append(f"Structural Class: {ligand.structural_class}")
        lines.append(f"")
        lines.append(f"Selectivity Trend: {ligand.selectivity_trend} ({ligand.trend_confidence} confidence)")

        trend_summary = ligand.get_trend_summary()
        lines.append(f"Direction: {trend_summary['direction']}")
        lines.append(f"Elements tested: {', '.join(ligand.element_coverage)} ({len(ligand.element_coverage)} total)")

    if format_type in ['full', 'stats']:
        stats = ligand.get_statistics()
        lines.append(f"")
        lines.append(f"Statistics:")
        lines.append(f"  D_min: {stats.get('D_min')}")
        lines.append(f"  D_max: {stats.get('D_max')}")
        lines.append(f"  D_mean: {stats.get('D_mean')}")
        lines.append(f"  D_std: {stats.get('D_std')}")
        lines.append(f"  log_D_mean: {stats.get('log_D_mean')}")
        lines.append(f"  log_D_std: {stats.get('log_D_std')}")
        lines.append(f"  Element count: {stats.get('element_count')}")

    if format_type in ['full', 'raw']:
        d_values = ligand.get_d_values()
        lines.append(f"")
        lines.append(f"Distribution Ratios (D):")
        for elem, d_val in sorted(d_values.items()):
            lines.append(f"  {elem}: {d_val:.3f}")

    lines.append("=" * 70)
    return "\n".join(lines)


def format_csv(ligand: LigandQuery, format_type: str) -> str:
    """Format ligand data as CSV"""
    if format_type == 'raw':
        # CSV format for D-values
        d_values = ligand.get_d_values()
        header = "Ligand,Element,D_value"
        rows = [header]
        for elem, d_val in sorted(d_values.items()):
            rows.append(f"{ligand.name},{elem},{d_val}")
        return "\n".join(rows)

    elif format_type == 'trend':
        header = "Ligand,Mechanism,Trend,Confidence,Elements_Tested"
        row = f"{ligand.name},{ligand.mechanism},{ligand.selectivity_trend},{ligand.trend_confidence},{len(ligand.element_coverage)}"
        return f"{header}\n{row}"

    elif format_type == 'stats':
        stats = ligand.get_statistics()
        header = "Ligand,D_min,D_max,D_mean,D_std,Element_Count"
        row = f"{ligand.name},{stats.get('D_min')},{stats.get('D_max')},{stats.get('D_mean')},{stats.get('D_std')},{stats.get('element_count')}"
        return f"{header}\n{row}"

    else:  # full
        # Combine all info
        stats = ligand.get_statistics()
        header = "Ligand,Mechanism,Trend,Confidence,D_mean,D_std,Elements_Tested"
        row = f"{ligand.name},{ligand.mechanism},{ligand.selectivity_trend},{ligand.trend_confidence},{stats.get('D_mean')},{stats.get('D_std')},{len(ligand.element_coverage)}"
        return f"{header}\n{row}"


def main():
    parser = argparse.ArgumentParser(
        description='Query rare earth ligand database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Query by name with full output
  python scripts/query_ligands.py --name TODGA --format full

  # Get raw D-values only
  python scripts/query_ligands.py --name TODGA --format raw --output table

  # Find all ascending trend ligands
  python scripts/query_ligands.py --trend Ascending --format trend

  # Search for DGA ligands
  python scripts/query_ligands.py --search DGA

  # Get statistics for specific ligand
  python scripts/query_ligands.py --name DMDODGA --format stats
        """
    )

    # Query type (mutually exclusive)
    query_group = parser.add_mutually_exclusive_group(required=True)
    query_group.add_argument('--name', type=str, help='Ligand name (e.g., TODGA)')
    query_group.add_argument('--smiles', type=str, help='Ligand SMILES string')
    query_group.add_argument('--mechanism', type=str,
                           choices=list(canonical_mechanisms()),
                           help='Filter by mechanism')
    query_group.add_argument('--trend', type=str,
                           choices=list(TREND_CLASSIFICATIONS),
                           help='Filter by selectivity trend')
    query_group.add_argument('--community-name', dest='community_name', type=str,
                           help='Filter by community-form Type label (e.g., DAPhen, TODGA, BTBP, CMPO)')
    query_group.add_argument('--search', type=str, help='Search term (searches name, SMILES, subtype, structural class)')
    query_group.add_argument('--list-all', action='store_true', help='List all ligands')
    query_group.add_argument('--stats', action='store_true', help='Show database statistics')

    # Output format
    parser.add_argument('--format', choices=['raw', 'trend', 'stats', 'full'], default='full',
                       help='Output format (default: full)')
    parser.add_argument('--output', choices=['json', 'table', 'csv'], default='json',
                       help='Output style (default: json)')
    parser.add_argument('--condition', type=str, help='Specific condition ID (optional, for --format raw/stats)')
    parser.add_argument('--fuzzy', action='store_true', help='Enable fuzzy name matching')

    args = parser.parse_args()

    # Initialize database
    try:
        db = LigandDatabase()
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        print("\nPlease ensure the database has been built first:")
        print("  python scripts/build_database.py")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to load database: {e}")
        sys.exit(1)

    # Handle database statistics
    if args.stats:
        stats = db.get_statistics()
        print("Database Statistics:")
        print("=" * 50)
        print(f"Total ligands: {stats['total_ligands']}")
        print(f"\nMechanisms:")
        for mech, count in sorted(stats['mechanisms'].items()):
            print(f"  {mech}: {count}")
        print(f"\nSelectivity Trends:")
        for trend, count in sorted(stats['trends'].items()):
            print(f"  {trend}: {count}")
        if stats['structural_classes']:   # legacy field; the public release omits it
            print(f"\nStructural Classes:")
            for sc, count in sorted(stats['structural_classes'].items()):
                print(f"  {sc}: {count}")
        return

    # Execute query
    results = []

    if args.name:
        results = db.query_by_name(args.name, fuzzy=args.fuzzy)
    elif args.smiles:
        result = db.query_by_smiles(args.smiles)
        results = [result] if result else []
    elif args.mechanism:
        results = db.filter_by_mechanism(args.mechanism)
    elif args.trend:
        results = db.filter_by_trend(args.trend)
    elif args.community_name:
        results = db.filter_by_community_name(args.community_name)
    elif args.search:
        results = db.search(args.search)
    elif args.list_all:
        results = db.get_all_ligands()

    # Handle no results
    if not results:
        print(f"No ligands found matching query")
        sys.exit(0)

    # Print header for multiple results
    if len(results) > 1:
        print(f"Found {len(results)} ligand(s):")
        print()

    # Format and print output
    for i, ligand in enumerate(results):
        if args.output == 'json':
            # Build output dict based on format
            if args.format == 'raw':
                d_values = ligand.get_d_values(args.condition)
                output = {'ligand': ligand.name, 'D_values': d_values}
            elif args.format == 'trend':
                trend_summary = ligand.get_trend_summary()
                output = {
                    'ligand': ligand.name,
                    'mechanism': ligand.mechanism,
                    'trend': trend_summary
                }
            elif args.format == 'stats':
                stats = ligand.get_statistics(args.condition)
                output = {'ligand': ligand.name, 'statistics': stats}
            else:  # full
                output = {
                    'ligand': ligand.name,
                    'mechanism': ligand.mechanism,
                    'structural_class': ligand.structural_class,
                    'trend': ligand.get_trend_summary(),
                    'statistics': ligand.get_statistics(args.condition),
                    'D_values': ligand.get_d_values(args.condition),
                    'conditions_available': ligand.get_condition_list()
                }

            print(format_json(output))

        elif args.output == 'table':
            print(format_table(ligand, args.format))

        elif args.output == 'csv':
            print(format_csv(ligand, args.format))

        # Add separator between multiple results
        if len(results) > 1 and i < len(results) - 1:
            print()


if __name__ == '__main__':
    main()
