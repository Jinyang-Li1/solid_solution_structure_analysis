"""
filter_solvent_entries.py - Filters disordered, solvent-free molecules data based on solvent terms and keywords.
It separates entries containing solvent-related information from those without, saving results to two CSV files.
"""

import csv
import re
from tqdm import tqdm

def load_solvent_terms(solvent_file):
    """Read solvent file and return a set of lowercase terms"""
    solvent_terms = set()
    try:
        with open(solvent_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header row
            for row in reader:
                for term in row:
                    term = term.strip()
                    if term:  # Skip empty values
                        solvent_terms.add(term.lower())
    except FileNotFoundError:
        print(f"Error: Solvent file {solvent_file} not found")
        exit()
    return solvent_terms

def process_data():
    # Configuration parameters
    input_file = "4. disorder details.csv"
    output_match = "6. disordered entry details with solvent.csv"
    output_unmatch = "6. disordered entry details without solvent.csv"
    solvent_file = "5. solvent.txt"

    # Load solvent data
    solvent_terms = load_solvent_terms(solvent_file)
    print(f"Loaded {len(solvent_terms)} solvent terms")

    # Precompile regex patterns (keyword matching)
    keyword_pattern = re.compile(r'\b(solvent|solvents|solvate|mask|olex2|squeeze|platon)\b', flags=re.IGNORECASE)

    # Build regex patterns for solvent names, including plural forms
    solvent_patterns = []
    for term in solvent_terms:
        # Create regex pattern for each solvent name, allowing optional plural 's'
        pattern = r'\b' + re.escape(term) + r's?\b'
        solvent_patterns.append(re.compile(pattern, flags=re.IGNORECASE))

    # Process CSV files
    try:
        with open(input_file, 'r', encoding='utf-8') as f_in, \
             open(output_match, 'w', encoding='utf-8', newline='') as f_match, \
             open(output_unmatch, 'w', encoding='utf-8', newline='') as f_unmatch:

            reader = csv.reader(f_in)
            writer_match = csv.writer(f_match)
            writer_unmatch = csv.writer(f_unmatch)

            # Process header row
            headers = next(reader)
            writer_match.writerow(headers)
            writer_unmatch.writerow(headers)

            # Get total rows for progress bar
            total = sum(1 for _ in reader)
            f_in.seek(0)
            next(reader)

            pbar = tqdm(total=total, desc="Processing progress", unit="row")
            match_count = 0

            for row in reader:
                if len(row) < 4:
                    continue  # Skip incomplete rows

                # Condition 1: Check if columns 4 onwards contain solvent terms
                has_solvent = any(
                    any(pattern.search(cell) for pattern in solvent_patterns)
                    for cell in row[3:]
                    if cell.strip()
                )

                # Condition 2: Check if any column contains keywords
                has_keyword = any(
                    keyword_pattern.search(cell)
                    for cell in row
                    if cell.strip()
                )

                # Write to corresponding file
                if has_solvent or has_keyword:
                    writer_match.writerow(row)
                    match_count += 1
                else:
                    writer_unmatch.writerow(row)

                pbar.update(1)

            pbar.close()
            print(f"\nMatching results: {match_count} rows matched in total")

    except Exception as e:
        print(f"Error occurred while processing files: {str(e)}")

if __name__ == "__main__":
    process_data()