"""
disordered_molecule_analysis.py - Analyzes disordered molecule data to identify disorder patterns and
classify molecular similarity.
"""

import csv
import re
from tqdm import tqdm

# Define keyword lists for disorder detection
DISORDER_TERMS = [
    'partially disordered', 'two-fold disordered', 'equally disordered',
    'each disordered', 'disordered', 'disorder', 'refined'
]
TWO_SITES_TERMS = [
    'two sites', 'two configurations', 'two positions disorder', 'two positions', 'two equal sites', 'two sets',
    'two possible sites', '2 positions', '2 configurations', '2 positions disorder', '2 equal sites',
    'between the two', 'disordered between two', 'two orientations', 'two enantiomers', 'twofold disordered',
    'mixed occupancy', 'alternate positions', 'twofold site', 'equal occupancy', 'over sites', '0.5:0.5', '1:1 ratio',
    'two orientations', '2 position', '2 site', '2 sites', '2 configuration', '2 equal site', 'two position',
    'one sites', 'one configurations', 'one positions', 'one equal sites', 'one site', 'one configuration',
    'one position', 'one equal site', '1 equal sites', '1 equal site',
    '1 sites', '1 configurations', '1 positions', '1 equal sites', '1 site', '1 configuration', '1 position',
]
THREE_SITES_TERMS = [
    'three sites', 'three configurations', 'three positions', 'three sets',
    'triple sites', 'three-fold disorder', 'three alternative positions', '3 ring sites',
    '3 positions', '3 sites', '3 configurations', '3 positions disorder', '3 equal sites'
]
FOUR_SITES_TERMS = [
    'four sites', 'four configurations', 'four positions', 'four sets',
    'quadruple sites', 'four-fold disorder', 'four alternative positions',
    '4 positions', '4 sites', '4 configurations', '4 positions disorder', '4 equal sites'
]
FIVE_SITES_TERMS = [
    'five sites', 'five configurations', 'five positions', 'five sets',
    'quintuple sites', 'five-fold disorder', 'five alternative positions',
    '5 positions', '5 sites', '5 configurations', '5 positions disorder', '5 equal sites'
]
OCCUPANCY_KEYS = ['occupancies', 'occupancy', 's.o.f.s']

# Compile regex patterns
# Pattern to detect disorder-related keywords
DISORDER_PATTERN = re.compile(
    r'(?i)\b(' + '|'.join(map(re.escape, DISORDER_TERMS)) + r')\b'
)
OCCUPANCY_PATTERN = re.compile(
    r'(?:occupancy of|occupation factors of|occupancies of|occupancies|occupancy|s\.o\.f\.s|ratio|mixed|disordered|equal|0\.\d+|1:1)\s*[:,\s-]?\s*((?:\d+\.?\d*[:,\s-]?)+)',
    flags=re.IGNORECASE
)
# Patterns to clean text (remove redundant terms/punctuation)
REMOVE_PATTERNS = [
    (r',\s*respectively', '', re.IGNORECASE),
    (r'(?i)\s*from\s+[^,]+,', ',', re.IGNORECASE),
    (r'\b(group|groups|on|both|atom|atoms|all|a|an|the|some|any|is|was|were|are|one)\b', '', re.IGNORECASE)
]

# Set of chemical element symbols (for molecular similarity analysis)
ELEMENT_SYMBOLS = {
    'H', 'Li', 'Be', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'P', 'S', 'Cl', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
}


def clean_text(text):
    """
        Clean input text by applying predefined regex patterns and normalizing punctuation/spaces.
        """
    if not isinstance(text, str):
        return ""  # Handle non-string inputs (e.g., None, numbers)

    # Apply all text cleaning patterns
    for pattern, repl, flags in REMOVE_PATTERNS:
        text = re.sub(pattern, repl, text, flags=flags)

    # Additional cleaning: merge consecutive commas and spaces
    text = re.sub(r'[,]+', ',', text) # Merge consecutive commas
    text = re.sub(r'\s+', ' ', text).strip() # Merge consecutive spaces and trim edges

    return text


def analyze_molecule(parts):
    """Analyze molecular similarity across split disorder segments."""
    if not parts:
        return None

    # Split each segment by comma (handle ", " or "," delimiters)
    split_parts = [re.split(r',\s?', part) for part in parts]
    # Get maximum length of split segments to align for comparison
    max_length = max(len(sub_parts) for sub_parts in split_parts)

    # Pad shorter split lists with empty strings to match max length
    for sub_parts in split_parts:
        sub_parts.extend([''] * (max_length - len(sub_parts)))

    # Compare each position across all segments
    for i in range(max_length):
        signatures = []
        for sub_parts in split_parts:
            part = sub_parts[i]
            if part:
                # Extract signature: prioritize element symbols, then leading letters
                first_two_chars = part[:2]
                if first_two_chars in ELEMENT_SYMBOLS:
                    signature = first_two_chars
                elif first_two_chars[0] in ELEMENT_SYMBOLS:
                    signature = first_two_chars[0]
                else:
                    # Extract only letters for non-element segments
                    only_letters = ''.join(filter(str.isalpha, first_two_chars))
                    # Early exit if current signature differs from existing ones
                    if signatures and ''.join(filter(str.isalpha, signatures[0])) != only_letters:
                        return 'different'
                    signature = only_letters
                signatures.append(signature)

        # If any position has mismatched signatures, molecules are different
        if signatures and not all(sig == signatures[0] for sig in signatures):
            return 'different'

    # All positions match: molecules are the same
    return 'same'


def process_cell(cell):
    """
    Process a single CSV cell to extract disorder segments, occupancy, and molecular similarity.
    """
    cell = clean_text(cell)
    parts = [cell]  # Default: single segment if no disorder detected
    occupancy = None
    same_molecule = None
    different_molecule = None
    is_matched = False

    # Check if cell contains any disorder-related keywords
    disorder_positions = [cell.lower().find(t) for t in DISORDER_TERMS if cell.lower().find(t) != -1]
    if disorder_positions:
        is_matched = True
        first_disorder_pos = min(disorder_positions)
        sub_str = cell[:first_disorder_pos]

        # Detect two-site patterns (occupancy format or keywords)
        occupancy_match = OCCUPANCY_PATTERN.search(cell)
        has_two_occupancy = occupancy_match and len(
            occupancy_match.group(1).split(':')) == 2 if occupancy_match else False
        has_two_sites_keyword = any(term.lower() in cell.lower() for term in TWO_SITES_TERMS)
        has_implicit_two = re.search(r'\btwo\b.*?(disorder|occupancy|ratio|mixed|alternate)', cell.lower())
        has_two_sites = has_two_sites_keyword or has_two_occupancy or bool(has_implicit_two)

        # Detect multi-site patterns (3/4/5 sites)
        has_three_sites = any(term.lower() in cell.lower() for term in THREE_SITES_TERMS)
        has_four_sites = any(term.lower() in cell.lower() for term in FOUR_SITES_TERMS)
        has_five_sites = any(term.lower() in cell.lower() for term in FIVE_SITES_TERMS)

        # Split segments and extract/infer occupancy based on site count
        if ' and ' in sub_str or has_two_sites:
            if ' and ' in sub_str:
                parts = [p.strip() for p in sub_str.split(' and ')]
                # Extract or infer occupancy
                if occupancy_match:
                    occupancy = occupancy_match.group(1)
                elif 'equal' in cell.lower():
                    occupancy = 'equal (0.5:0.5)'
            else:
                # Duplicate pre-disorder text for two-site cases without "and"
                before_part = sub_str.strip()
                parts = [before_part] * 2
                if occupancy_match:
                    occupancy = occupancy_match.group(1)
                elif 'equal' in cell.lower():
                    occupancy = 'equal (0.5:0.5)'

        elif has_three_sites:
            before_part = sub_str.strip()
            parts = [before_part] * 3
            num_match = OCCUPANCY_PATTERN.search(cell)
            if num_match:
                occupancy = num_match.group(1)
            if not occupancy:
                if 'equal' in cell.lower() or 'equatorial' in cell.lower():
                    occupancy = 'equal (0.33:0.33:0.33)'

        elif has_four_sites:
            before_part = sub_str.strip()
            parts = [before_part] * 4
            num_match = OCCUPANCY_PATTERN.search(cell)
            if num_match:
                occupancy = num_match.group(1)
            if not occupancy:
                if 'equal' in cell.lower() or 'equatorial' in cell.lower():
                    occupancy = 'equal (0.25:0.25:0.25:0.25)'

        elif has_five_sites:
            before_part = sub_str.strip()
            parts = [before_part] * 5
            num_match = OCCUPANCY_PATTERN.search(cell)
            if num_match:
                occupancy = num_match.group(1)
            if not occupancy:
                if 'equal' in cell.lower() or 'equatorial' in cell.lower():
                    occupancy = 'equal (0.2:0.2:0.2:0.2:0.2)'

        # Analyze molecular similarity (filter out empty segments)
        valid_parts = [p for p in parts if p]
        if valid_parts:
            similarity = analyze_molecule(valid_parts)
            same_molecule = 'same' if similarity == 'same' else None
            different_molecule = 'different' if similarity != 'same' else None

    return parts, occupancy, same_molecule, different_molecule, is_matched


def process_file():
    """Main function: process input CSV, generate analysis results, and export outputs."""
    input_file = "8. disordered entries with occupancy split.csv"
    unmatched_file = "9. unmatched rows atom disorder analysis.csv"
    same_refcode_file = "9. refcodes all same.csv"
    different_refcode_file = "9. refcodes any different.csv"
    unmatched_rows = []

    # Cache to track final status of each refcode (priority: different > same > None)
    refcode_status = {}

    try:
        with open(input_file, 'r', newline='') as infile:
            reader = csv.reader(infile)
            header = next(reader)
            data = [row for row in reader]

        # Track which rows are matched across any column
        row_matched = [False] * len(data)

        # Process columns 4-15 (0-indexed: 3 to 14)
        for col_index in tqdm(range(3, 15), desc="Processing columns"):
            # Step 1: Calculate max number of segments for dynamic column headers
            max_parts = 0
            for row in data:
                cell = row[col_index]
                if DISORDER_PATTERN.search(cell):
                    parts, _, _, _, _ = process_cell(cell)
                    max_parts = max(max_parts, len(parts))

            # Dynamic header: original column + split segments + analysis results
            dynamic_header = header[:3] + [header[col_index]] + [f"Part {i + 1}" for i in range(max_parts)] + [
                "occupancy", "same molecule", "different molecule"
            ]

            processed_data = []

            # Step 2: Process each row in the current column
            for row_idx, row in enumerate(tqdm(data, desc=f"Processing column {col_index + 1} ", leave=False)):
                refcode = row[0]
                cell = row[col_index]
                parts, occupancy, same, different, is_matched = process_cell(cell)

                # Mark row as matched if valid analysis is available
                current_column_matched = is_matched and (same is not None or different is not None)
                if current_column_matched:
                    row_matched[row_idx] = True

                # Prepare processed row (pad segments to match max length)
                if current_column_matched:
                    parts += [''] * (max_parts - len(parts))
                    processed_row = row[:3] + [row[col_index]] + parts + [
                        occupancy or '',
                        same or '',
                        different or ''
                    ]
                    processed_data.append(processed_row)

                # Update refcode status (prioritize "different" over "same")
                if current_column_matched:
                    if different == 'different':
                        refcode_status[refcode] = 'different'
                    elif same == 'same' and refcode not in refcode_status:
                        refcode_status[refcode] = 'same'


            def sort_key(row):
                """Sort key: 'same' (0) > 'different' (1) > other (2)"""
                return 0 if row[-2] == 'same' else 1 if row[-1] == 'different' else 2

            processed_data.sort(key=sort_key)

            # Export column-specific results
            output_file = f"9. column {col_index + 1} disorder analysis.csv"
            with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
                writer = csv.writer(outfile)
                writer.writerow(dynamic_header)
                writer.writerows(processed_data)

        # Step 3: Collect and export unmatched rows (no matches in any column)
        for row_idx, row in enumerate(data):
            if not row_matched[row_idx]:
                unmatched_rows.append(row)
        with open(unmatched_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(header)
            writer.writerows(unmatched_rows)

        #  Step 4: Split refcodes by final status (ensure mutual exclusivity)
        same_refcodes = set()
        different_refcodes = set()
        for refcode, status in refcode_status.items():
            if status == 'same':
                same_refcodes.add(refcode)
            elif status == 'different':
                different_refcodes.add(refcode)

        # Export "all same" refcodes
        same_file_path = same_refcode_file
        with open(same_file_path, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(["refcode"])
            for refcode in sorted(same_refcodes):
                writer.writerow([refcode])

        # Export "any different" refcodes
        different_file_path = different_refcode_file
        with open(different_file_path, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(["refcode"])
            for refcode in sorted(different_refcodes):
                writer.writerow([refcode])

        # 输出统计信息
        print(f"\nProcessing Complete!")
        print(f"1. Column-specific results: '9. column X disorder analysis.csv' (columns 4-15)")
        print(f"2. Unmatched rows: {unmatched_file} ({len(unmatched_rows)} rows)")
        print(f"3. Refcodes with all 'same' matches: {same_refcode_file} ({len(same_refcodes)} refcodes)")
        print(
            f"4. Refcodes with any 'different' matches: {different_refcode_file} ({len(different_refcodes)} refcodes)")
        print(f"Note: Refcodes with both 'same' and 'different' matches are classified as 'different'.")


    except FileNotFoundError:

        print(f"Error: Input file '{input_file}' not found. Check the file path.")

    except Exception as e:

        print(f"Unexpected error during processing: {str(e)}.")


if __name__ == "__main__":
    process_file()