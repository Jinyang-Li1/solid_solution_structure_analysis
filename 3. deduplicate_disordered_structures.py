"""
    deduplicate_disordered_structures.py - Based on the refcode, remove duplicate entries from the unordered structure.
"""


import csv

unique_codes = {}

with open('2. crystal structures with disorder.csv', 'r', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)

    for row in reader:
        identifier = row[0]
        if len(identifier) > 6 and identifier[6:].isdigit():
            code = identifier[:6]
            number = int(identifier[6:])
            if code in unique_codes:
                existing_identifier = unique_codes[code][0]
                if len(existing_identifier) > 6 and existing_identifier[6:].isdigit():
                    existing_number = int(existing_identifier[6:])
                    if number < existing_number:
                        unique_codes[code] = row
            else:
                unique_codes[code] = row
        else:
            unique_codes[identifier] = row

with open('3. disorder structure deduplicate.csv', 'w', newline='',
          encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    for row in unique_codes.values():
        writer.writerow(row)