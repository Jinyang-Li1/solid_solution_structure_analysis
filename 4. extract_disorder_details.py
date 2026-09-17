"""
    extract_disorder_details.py - Extract disorder details from CSD for identifiers.
"""

import csv
from ccdc import io
from tqdm import tqdm

# Read input CSV with CSD identifiers
with open('3. disorder structure deduplicate.csv', 'r', encoding='utf-8') as infile:
    csv_reader = csv.reader(infile)
    headers = next(csv_reader)
    # Define output headers (first 3 columns + disorder detail)
    new_headers = headers[:3] + ['disorder detail']

    csd_reader = io.CrystalReader('CSD')
    data_rows = list(csv_reader)  # Load all data rows

    # Write output CSV with extracted disorder details
    with open('4. disorder details.csv', 'w', newline='', encoding='utf-8') as outfile:
        csv_writer = csv.writer(outfile)
        csv_writer.writerow(new_headers)  # Write header

        # Process each row with progress tracking
        for row in tqdm(data_rows, desc="Processing rows"):
            identifier = row[0]
            try:
                # Retrieve CSD entry and extract full disorder details
                entry = csd_reader.entry(identifier)
                disorder_details = str(entry.disorder_details).strip()
                # Construct output row with preserved columns + disorder details
                new_row = row[:3] + [disorder_details]
                csv_writer.writerow(new_row)
            except Exception as e:
                print(f"Error processing {identifier}: {e}")
                # Write row with empty disorder detail for invalid entries
                csv_writer.writerow(row[:3] + [''])