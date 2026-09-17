"""
    search_in_csd.py - Filters entries that meet specific criteria from the CSD.
"""

from ccdc import io
from ccdc.search import Search
from ccdc.solid_form import SolvateAnalyser
import csv
from tqdm import tqdm
from rdkit import Chem

# -------------------------- Base filtering criteria --------------------------
base_settings = Search.Settings()
# Core structure filters
base_settings.has_3d_coordinates = True  # With 3D coordinates
base_settings.only_organic = True        # Organic compounds only
base_settings.not_polymeric = True       # Non-polymeric
base_settings.no_powder = True           # No powder diffraction structures
# Data quality filters
base_settings.max_r_factor = 15          # R-factor ≤ 15 (adjust as needed)
# Element filters
base_settings.no_metals = True           # No metals
base_settings.must_not_have_elements = [
    'As', 'Te', 'At', 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn',  # Toxic/rare gases
    'B', 'Al', 'Ga', 'In', 'Tl', 'Si', 'Ge', 'Sn', 'Pb', 'Sb', 'Po'  # Specific main group elements
]

# -------------------------- Disorder structure settings --------------------------
disorder_settings = Search.Settings()
disorder_settings.no_disorder = True  # Test pass → no disorder; Test fail → has disorder


def search_and_filter_csd():
    """Perform CSD search and filtering, output solvent-free crystal libraries"""
    # Initialize CSD reader
    csd_reader = io.EntryReader("CSD")

    # Output file paths (modify as needed)
    output_paths = {
        "all_no_solvent": "1. small_molecule_crystals_no_solvent.csv",
        "disorder_no_solvent": "1. disordered_crystals_no_solvent.csv"
    }

    # Open output files and initialize CSV writers
    with open(output_paths["all_no_solvent"], 'w', newline='', encoding='utf-8') as f_all, \
         open(output_paths["disorder_no_solvent"], 'w', newline='', encoding='utf-8') as f_disorder:

        writer_all = csv.writer(f_all)
        writer_disorder = csv.writer(f_disorder)

        # Write headers
        header = ["Identifier", "SMILES", "Canonical_SMILES"]
        writer_all.writerow(header)
        writer_disorder.writerow(header)

        # Iterate through CSD entries with progress bar
        for entry in tqdm(csd_reader, desc="Searching & Filtering CSD", unit="entry"):
            # 1. Filter out entries that don't meet base criteria
            if not base_settings.test(entry):
                continue

            # 2. Filter out high molecular weight entries (MW ≤ 1300)
            if entry.molecule.molecular_weight > 1300:
                continue

            try:
                # 3. Check for solvates (core: using SolvateAnalyser)
                crystal = csd_reader.crystal(entry.identifier)
                solvate_analyser = SolvateAnalyser(crystal)
                detected_solvents = solvate_analyser.find_solvents()
                if detected_solvents:  # Skip if solvents detected
                    continue

                # 4. Critical fix: Check for valid SMILES (prevent NoneType errors)
                raw_smiles = entry.molecule.smiles
                if not raw_smiles:  # Covers None or empty string cases
                    print(f"Entry {entry.identifier}: No valid SMILES, skipping")
                    continue

                # 5. Generate canonical SMILES
                rdkit_mol = Chem.MolFromSmiles(raw_smiles)
                canonical_smiles = Chem.MolToSmiles(rdkit_mol, canonical=True) if rdkit_mol else ""

                # 6. Write to files: First to "all solvent-free crystals"
                writer_all.writerow([entry.identifier, raw_smiles, canonical_smiles])

                # 7. Identify disordered structures: If disorder test fails → has disorder
                if not disorder_settings.test(entry):
                    writer_disorder.writerow([entry.identifier, raw_smiles, canonical_smiles])

            # Catch all exceptions to avoid program interruption (only print errors)
            except Exception as err:
                print(f"Error processing entry {entry.identifier}: {str(err)}")


if __name__ == "__main__":
    search_and_filter_csd()

