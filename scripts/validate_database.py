#!/usr/bin/env python3
"""
Validate Database
Checks every record with rel_quest.validators.validate_ligand_record (cross-field
invariants) and, for current-version records, against the JSON Schema in
schema/ligand_record.schema.json.
"""

import os
import sys
import argparse

# Add parent directory to path to import rel_quest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rel_quest.validators import validate_database_file


def main():
    parser = argparse.ArgumentParser(
        description='Validate REEL-QUEST database file',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Default to the canonical database (same file LigandDatabase resolves to)
    db_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    default_db = os.path.join(db_dir, 'data', 'sx_ligands_merged_v14_flat.jsonl')

    parser.add_argument('--database', '-d', type=str, default=default_db,
                       help=f'Path to JSONL database (default: {default_db})')
    parser.add_argument('--max-records', '-m', type=int, default=None,
                       help='Maximum number of records to validate (default: all)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Show all error details')
    parser.add_argument('--no-json-schema', action='store_true',
                       help='Skip the JSON Schema check (schema/ligand_record.schema.json)')

    args = parser.parse_args()

    # Check database exists
    if not os.path.exists(args.database):
        print(f"ERROR: Database file not found: {args.database}")
        print("\nPlease build the database first:")
        print("  python scripts/build_database.py")
        sys.exit(1)

    print("="*60)
    print("REEL-QUEST Database Validation")
    print("="*60)
    print(f"Database: {args.database}")
    if args.max_records:
        print(f"Max records: {args.max_records}")
    print("="*60)
    print()

    # Run validation
    valid_count, invalid_count, error_summary = validate_database_file(
        args.database,
        max_records=args.max_records,
        use_json_schema=not args.no_json_schema,
    )

    total = valid_count + invalid_count

    print(f"Validation Results:")
    print(f"  Total records: {total}")
    print(f"  Valid: {valid_count}")
    print(f"  Invalid: {invalid_count}")

    if invalid_count > 0:
        print(f"\nValidation FAILED - {invalid_count} error(s) found")
        if args.verbose or invalid_count <= 10:
            print("\nError Details:")
            for error in error_summary:
                print(f"  - {error}")
        else:
            print(f"\nShowing first 10 errors (use --verbose for all):")
            for error in error_summary[:10]:
                print(f"  - {error}")
        sys.exit(1)
    else:
        print(f"\nValidation PASSED - All records are valid")
        sys.exit(0)


if __name__ == '__main__':
    main()
