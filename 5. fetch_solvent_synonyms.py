"""
fetch_solvent_synonyms.py - A script to fetch the top 5 English-containing synonyms for a list of solvents from PubChem.

"""

import pubchempy as pcp
from tqdm import tqdm
import time
import csv
import re

def get_top5_synonyms():
    """
        Fetch top 5 English-containing synonyms for each solvent from PubChem.

        Reads solvent names from 'solvent.txt', queries PubChem for each,
        filters synonyms to include only those with English letters,
        deduplicates, takes top 5, and saves to 'solvent_top5.csv'.

        Returns:
            set: All unique English-containing synonyms collected
        """
    # Read solvent names from input file
    with open('solvent.txt', 'r') as f:
        solvents = [line.strip() for line in f if line.strip()]

    synonym_db = {}
    print("\n[Fetching top 5 English-containing synonyms for solvents]")
    for solvent in tqdm(solvents, desc="Querying PubChem"):
        max_retries = 3 # Max retry attempts for API calls
        retries = 0
        while retries < max_retries:
            try:
                compounds = pcp.get_compounds(solvent, 'name')
                all_synonyms = []
                for comp in compounds:
                    if comp.synonyms:
                        for syn in comp.synonyms:
                            syn = syn.lower()
                            #  Filter synonyms with at least one English letter
                            if re.search('[a-zA-Z]', syn):
                                all_synonyms.append(syn)
                # 去重并取前5个
                top5 = list(dict.fromkeys(all_synonyms))[:5]
                synonym_db[solvent] = top5

                time.sleep(0.5)  # Rate limit to avoid overwhelming PubChem API
                break  # Exit retry loop on success
            except Exception as e:
                retries += 1
                if retries < max_retries:
                    tqdm.write(f"Retry {retries} for {solvent} failed: {str(e)}")
                    time.sleep(1)  # Longer wait before next retry
                else:
                    tqdm.write(f"Max retries reached for {solvent}: {str(e)}")
                    synonym_db[solvent] = []

    # Save results to CSV
    with open('solvent_top5.csv', 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['solvent', '1', '2', '3', '4', '5'])
        for name, syns in synonym_db.items():
            row = [name] + syns + [''] * (5 - len(syns))
            writer.writerow(row)

    # Collect all unique synonyms across all solvents
    all_synonyms = set()
    for syns in synonym_db.values():
        all_synonyms.update(syns)
    return all_synonyms


if __name__ == "__main__":
    # Execute main function and get results
    collected_synonyms = get_top5_synonyms()

    # Print summary
    print(f"\nCollected {len(collected_synonyms)} unique solvent synonyms")
    print("\nProcessing complete! Output file:")
    print(" - solvent_top5.csv (top 5 synonyms for each solvent)")
