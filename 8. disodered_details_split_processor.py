'''分割
disodered_details_split_processor.py - Processes CSV data of disordered, solvent-free molecules
that contain disorder occupancy information.
'''

import csv
import re
from ccdc import io
from tqdm import tqdm


def split_sentence(sentence):
    """Split sentences using semicolons and periods followed by space+uppercase letter (to avoid splitting abbreviations)"""
    if not sentence:
        return []

    # Regex pattern: split at semicolons or periods followed by space and uppercase letter
    pattern = re.compile(r';|\.(?=\s+[A-Z])')

    # Split text and filter empty strings
    parts = [part.strip() for part in pattern.split(sentence) if part.strip()]

    # Remove residual trailing periods from segments
    processed_parts = []
    for part in parts:
        if part.endswith('.'):
            processed_parts.append(part[:-1].strip())
        else:
            processed_parts.append(part)

    return processed_parts


def process_data():
    input_file = '7. disordered entries with occupancy.csv'
    output_file = '8. disordered entries with occupancy split.csv'

    # Read header row from input file
    with open(input_file, 'r', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        headers = next(reader)

    csd_reader = io.CrystalReader('CSD')
    max_new_columns = 0

    # First pass: determine maximum number of split disorder segments
    with open(input_file, 'r', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        next(reader)
        rows = list(reader)

        for row in tqdm(rows, desc="Analyzing max segments"):
            identifier = row[0]
            try:
                entry = csd_reader.entry(identifier)
                details = str(entry.disorder_details) if entry.disorder_details else ""
                parts = split_sentence(details)
                max_new_columns = max(max_new_columns, len(parts))
            except Exception as e:
                print(f"Error：{identifier} - {str(e)}")

    # Generate dynamic headers for split disorder details
    new_headers = headers[:3] + [f'disorder{i + 1}' for i in range(max_new_columns)]

    # Second pass: process and write structured data
    with open(input_file, 'r', encoding='utf-8') as infile, \
            open(output_file, 'w', newline='', encoding='utf-8') as outfile:

        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        next(reader)
        writer.writerow(new_headers)

        for row in tqdm(reader, desc="Processing entries", unit="entry"):
            identifier = row[0]
            try:
                entry = csd_reader.entry(identifier)
                details = str(entry.disorder_details) if entry.disorder_details else ""
                parts = split_sentence(details)

                # Build row with aligned columns (pad with empty strings if needed)
                new_row = row[:3] + parts + [''] * (max_new_columns - len(parts))
                writer.writerow(new_row)

            except Exception as e:
                print(f"Error：{identifier} - {str(e)}")
                # Mark entries with retrieval errors
                error_row = row[:3] + ['ERROR'] + [''] * (max_new_columns - 1)
                writer.writerow(error_row)


if __name__ == "__main__":
    process_data()