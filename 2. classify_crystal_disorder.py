"""
    classify_crystal_disorder.py - Check and classify crystals for disordered parts.
"""

import csv
from ccdc import io
import traceback

# Configuration parameters
INPUT_CSV = "1. disordered_crystals_no_solvent.csv"
OUTPUT_HAS_DISORDER = "2. crystal structures with disorder.csv"
OUTPUT_NO_DISORDER = "2. crystal structures without disorder.csv"
IDENTIFIER_COLUMN = "Identifier"


def has_disorder_structure(identifier: str) -> bool:
    """Check if a crystal structure contains valid disordered groups (with atoms)"""
    try:
        with io.CrystalReader("CSD") as cry_reader:
            crystal = cry_reader.crystal(identifier)

        # Check for existence of disorder assemblies
        if not (crystal.disorder and crystal.disorder.assemblies):
            return False

        # Verify if any disorder group contains atoms
        for assembly in crystal.disorder.assemblies:
            for group in assembly.groups:
                if group.atoms:
                    return True

        return False

    except Exception:
        return False


def main():
    total = 0
    has_disorder_count = 0
    no_disorder_count = 0

    try:
        with open(INPUT_CSV, "r", encoding="utf-8", newline="") as in_file:
            reader = csv.DictReader(in_file)
            fieldnames = reader.fieldnames

            if IDENTIFIER_COLUMN not in fieldnames:
                print(f"Available columns: {', '.join(fieldnames)}")
                raise ValueError(f"Identifier column [{IDENTIFIER_COLUMN}] not found")

            with open(OUTPUT_HAS_DISORDER, "w", encoding="utf-8", newline="") as has_file, \
                    open(OUTPUT_NO_DISORDER, "w", encoding="utf-8", newline="") as no_file:

                has_writer = csv.DictWriter(has_file, fieldnames=fieldnames)
                no_writer = csv.DictWriter(no_file, fieldnames=fieldnames)
                has_writer.writeheader()
                no_writer.writeheader()

                for row in reader:
                    total += 1
                    identifier = row[IDENTIFIER_COLUMN].strip()

                    if not identifier:
                        continue

                    if has_disorder_structure(identifier):
                        has_writer.writerow(row)
                        has_disorder_count += 1
                    else:
                        no_writer.writerow(row)
                        no_disorder_count += 1

                    # Print progress every 1000 rows
                    if total % 1000 == 0:
                        print(
                            f"Processed {total} rows | With disorder: {has_disorder_count} | Without disorder: {no_disorder_count}")

    except FileNotFoundError:
        print(f"Input file [{INPUT_CSV}] not found")
    except Exception as e:
        print(f"Program error: {str(e)}")
        traceback.print_exc()

    # Final statistics
    print("\n" + "=" * 60)
    print(f"Processing complete | Total rows: {total}")
    print(f"→ With disorder: {has_disorder_count} rows (saved to {OUTPUT_HAS_DISORDER})")
    print(f"→ Without disorder: {no_disorder_count} rows (saved to {OUTPUT_NO_DISORDER})")
    print("=" * 60)


if __name__ == "__main__":
    main()
