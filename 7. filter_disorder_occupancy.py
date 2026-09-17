"""
filter_disorder_occupancy.py - Filters disordered, solvent-free molecules (without solvent clues)
based on the presence of disorder occupancy information in the disorder details column.
Splits data into two CSV files: one containing entries with occupancy information,
and another with entries lacking such information.
"""

import csv
import re
from tqdm import tqdm


def is_valid_number(num_str):
    try:
        num = float(num_str)
        return 0 <= num <= 1
    except ValueError:
        return False


def filter_rows():
    # File paths (input and outputs)
    input_file = '6. disordered entry details without solvent.csv'
    output_file_with_ratio = '7. disordered entries with occupancy.csv'
    output_file_without_ratio = '7. disordered entries without occupancy.csv'

    # Read input data
    with open(input_file, 'r', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        headers = next(reader)
        rows = list(reader)
    # Write filtered results to output files
    with open(output_file_with_ratio, 'w', newline='', encoding='utf-8') as outfile_with_ratio, \
            open(output_file_without_ratio, 'w', newline='', encoding='utf-8') as outfile_without_ratio:
        writer_with_ratio = csv.writer(outfile_with_ratio)
        writer_without_ratio = csv.writer(outfile_without_ratio)

        # Write header rows to both output files
        writer_with_ratio.writerow(headers)
        writer_without_ratio.writerow(headers)

        # Process each row with progress bar
        for row in tqdm(rows, desc="Filtering data", unit="row"):
            if len(row) >= 4:  # Ensure row has at least 4 columns (disorder details in column 4)
                fourth_column = row[3]  # 4th column (index 3) contains disorder details

                # Check for various patterns indicating occupancy information:
                # 1. Valid numbers between 0 and 1 (common for occupancy ratios)
                has_valid_num = any(is_valid_number(num) for num in re.findall(r'\d+\.?\d*', fourth_column))

                # 2. Percentage symbols (e.g., "50% occupancy")
                has_percentage = '%' in fourth_column

                # 3. Percentage-related terms (e.g., "per cent", "percent")
                has_percent_terms = any(term in fourth_column.lower() for term in ['per cent', 'percent'])

                # 4. Terms indicating equal/partial occupancy (e.g., "half-occupancy", "equal", or other custom descriptions)
                has_equal_terms = any(term in fourth_column.lower() for term in ['equal', 'equatorial',
                                                                                 'half-occupancy', 'half occupancy',
                                                                                 'half site occupancy', 'half-occupied',
                                                                                 'half occupancies', 'major occupancy (78.7)','ratio'])

                # 5. Ratio formats (e.g., "1:1", "2:3")
                has_number_colon_number = re.search(r'\d+:\d+', fourth_column)

                # 6. Fraction formats (e.g., "1/2", "3/4")
                has_number_slash_number = re.search(r'\d+/\d+', fourth_column)

                # 7. Decimals without leading digits (e.g., ".5" for 0.5)
                has_dot_without_leading_digit = re.search(r'(?<!\d)\.\d+', fourth_column)

                # If any occupancy pattern is found, write to "with occupancy" file
                if (has_valid_num or has_percentage or has_percent_terms or has_equal_terms or has_number_colon_number
                    or has_number_slash_number or has_dot_without_leading_digit):
                    writer_with_ratio.writerow(row)
                else:
                    writer_without_ratio.writerow(row)


if __name__ == "__main__":
    filter_rows()