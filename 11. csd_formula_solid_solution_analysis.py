"""
csd_formula_different_same_analysis.py - Classifies CSD refcodes into different/same formulas
by comparing molecular formulas of regular and disordered molecule groups.
"""

import csv
from ccdc.io import EntryReader


def classify_formula_comparison(input_csv):
    """Classify refcodes into same/different formulas via CSD formula comparison"""
    same_path = "11. same_molecular_formula_in_same_atom_disorder.csv"
    different_path = "11. different_molecular_formula_in_same_atom_disorder.csv"
    output_headers = ['refcode', 'molecule_formula', 'disordered_molecule_formula']

    # Initialize CSD reader
    try:
        csd_reader = EntryReader('CSD')
    except Exception as e:
        print(f"CSD EntryReader init failed: {e}")
        return

    # Step 1: Read and filter valid refcodes (single file read)
    try:
        with open(input_csv, 'r', encoding='utf-8-sig') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            refcodes = [row[0].strip() for row in reader if row and row[0].strip()]

        if not refcodes:
            print("No valid refcodes found (check input CSV)")
            return
        total = len(refcodes)
    except FileNotFoundError:
        print(f"Input file not found: {input_csv}")
        return
    except Exception as e:
        print(f"Failed to read input: {e}")
        return

    # Step 2: Process refcodes
    error_count = 0
    with open(same_path, 'w', newline='', encoding='utf-8-sig') as same_f, \
            open(different_path, 'w', newline='', encoding='utf-8-sig') as different_f:

        # Init writers
        same_writer = csv.DictWriter(same_f, fieldnames=output_headers)
        different_writer = csv.DictWriter(different_f, fieldnames=output_headers)
        same_writer.writeheader()
        different_writer.writeheader()

        # Process each refcode
        for refcode in refcodes:
            try:
                # Get formulas from CSD
                entry = csd_reader.entry(refcode)
                reg_formula = entry.molecule.formula.strip()
                dis_formula = entry.disordered_molecule.formula.strip()

                # Classify and write (same: same molecular formula; different: different molecular formula)
                data = {
                    'refcode': refcode,
                    'molecule_formula': reg_formula,
                    'disordered_molecule_formula': dis_formula
                }
                if reg_formula == dis_formula:
                    same_writer.writerow(data)
                    print(f"Processed: {refcode} → Same formula")
                else:
                    different_writer.writerow(data)
                    print(f"Processed: {refcode} → Different formula")

            except Exception as e:
                error_count += 1
                err_msg = f"{str(e)[:30]}..." if len(str(e)) > 30 else str(e)
                print(f"Error processing {refcode}: {err_msg}")

    # Print summary
    success = total - error_count
    print(f"\n===== Done =====")
    print(f"Total Refcodes: {total} | Success: {success} | Errors: {error_count}")
    print(f"\nOutputs:")
    print(f"- Same formula: {same_path}")
    print(f"- Different formula: {different_path}")


if __name__ == "__main__":
    # Run classification (configure input path here)
    INPUT_CSV = "10. different smiles refcodes.csv"
    classify_formula_comparison(INPUT_CSV)
