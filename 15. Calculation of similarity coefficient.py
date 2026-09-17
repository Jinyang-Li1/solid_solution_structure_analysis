import numpy as np
from ccdc import io
from collections import defaultdict, Counter
import re
import pandas as pd
import os


class AtomInfo:
    """Encapsulates atomic information (supports atoms in disordered groups)"""

    def __init__(self, ccdc_atom, occupancy=None):
        self.label = ccdc_atom.label
        self.atomic_symbol = ccdc_atom.atomic_symbol
        self.coordinates = ccdc_atom.coordinates
        self.occupancy = occupancy if occupancy is not None else ccdc_atom.occupancy
        self.vdw_radius = ccdc_atom.vdw_radius  # Directly use the van der Waals radius attribute from the CCDC library

    def __eq__(self, other):
        """Overrides equality check for atom deduplication"""
        if not isinstance(other, AtomInfo):
            return False
        coord_eq = (self.coordinates.x == other.coordinates.x) if (self.coordinates and other.coordinates) else False
        return (self.label == other.label and
                self.atomic_symbol == other.atomic_symbol and
                coord_eq)

    def __hash__(self):
        """Overrides hash method for set-based deduplication"""
        coord_hash = (self.coordinates.x, self.coordinates.y, self.coordinates.z) if self.coordinates else ()
        return hash((self.label, self.atomic_symbol, coord_hash))


def print_atom_info(atom):
    """Safely prints atomic information (for debugging)"""
    if not atom.coordinates:
        occupancy_str = "Missing" if atom.occupancy is None else f"{atom.occupancy:.2f}"
        print(
            f"Label: {atom.label} | Element: {atom.atomic_symbol} | Coordinates: Missing | Occupancy: {occupancy_str} | van der Waals Radius: {atom.vdw_radius:.3f}Å")
        return
    occupancy_str = "Missing" if atom.occupancy is None else f"{atom.occupancy:.2f}"
    print(
        f"Label: {atom.label:5} | Element: {atom.atomic_symbol:2} | Coordinates: ({atom.coordinates.x:.3f},{atom.coordinates.y:.3f},{atom.coordinates.z:.3f}) | Occupancy: {occupancy_str} | van der Waals Radius: {atom.vdw_radius:.3f}Å")


def find_affecting_atoms(target_atoms, all_atoms, cutoff=5.0):
    """Finds surrounding atoms influencing target atoms (occupancy=1 and distance < sum of van der Waals radii)"""
    affecting = set()
    for target in target_atoms:
        if not target.coordinates:
            continue
        target_coord = np.array([target.coordinates.x, target.coordinates.y, target.coordinates.z])
        for atom in all_atoms:
            if atom == target or atom.occupancy != 1.0 or not atom.coordinates:
                continue
            dist = np.linalg.norm(target_coord - np.array([atom.coordinates.x, atom.coordinates.y, atom.coordinates.z]))
            if dist < (target.vdw_radius + atom.vdw_radius):
                affecting.add(atom)
    return list(affecting)


def calculate_sphere_volume(radius_angstrom):
    """Calculates the van der Waals sphere volume of a single atom (Unit: Å³)"""
    return (4 / 3) * np.pi * (radius_angstrom ** 3)


def calculate_overlap_volume(r1, r2, distance):
    """Calculates the overlapping volume of two spheres (Unit: Å³)"""
    if distance >= r1 + r2:
        return 0.0  # Spheres are separated, no overlap
    if distance <= abs(r1 - r2):
        return calculate_sphere_volume(min(r1, r2))  # One sphere fully contains the other, take volume of smaller sphere
    # Formula for overlapping volume when spheres intersect
    h1 = r1 - (r1 ** 2 - r2 ** 2 + distance ** 2) / (2 * distance)
    h2 = r2 - (r2 ** 2 - r1 ** 2 + distance ** 2) / (2 * distance)
    return (1 / 6) * np.pi * (h1 ** 2 * (3 * r1 - h1) + h2 ** 2 * (3 * r2 - h2))


def calculate_analytic_vdw_volume(atom_group):
    """Calculates van der Waals volume of a group using analytical geometry (accounts for intramolecular atomic overlap, Unit: Å³)"""
    total_volume = 0.0
    overlap_volume = 0.0
    # Filter atoms without coordinates/valid van der Waals radius
    valid_atoms = [a for a in atom_group if a.coordinates and a.vdw_radius and a.vdw_radius > 0]
    if not valid_atoms:
        return 0.0

    for i in range(len(valid_atoms)):
        atom1 = valid_atoms[i]
        r1 = atom1.vdw_radius
        total_volume += calculate_sphere_volume(r1)

        for j in range(i + 1, len(valid_atoms)):
            atom2 = valid_atoms[j]
            r2 = atom2.vdw_radius
            # Calculate interatomic distance
            coord1 = np.array([atom1.coordinates.x, atom1.coordinates.y, atom1.coordinates.z])
            coord2 = np.array([atom2.coordinates.x, atom2.coordinates.y, atom2.coordinates.z])
            distance = np.linalg.norm(coord1 - coord2)
            overlap_volume += calculate_overlap_volume(r1, r2, distance)

    return max(total_volume - overlap_volume, 0.0)  # Ensure non-negative volume


def calculate_vdw_volume(atom_group, env_atoms, consider_surrounding=True, n_points=50000):
    """Calculates van der Waals volume using Monte Carlo method (accounts for steric hindrance from surrounding atoms, Unit: Å³)"""
    valid_target = [atom for atom in atom_group if atom.coordinates and atom.vdw_radius and atom.vdw_radius > 0]
    if not valid_target:
        return 0.0

    # Calculate bounding box for target atoms
    target_pos = np.array([[a.coordinates.x, a.coordinates.y, a.coordinates.z] for a in valid_target])
    target_rad = np.array([a.vdw_radius for a in valid_target])
    min_corner = np.min(target_pos - target_rad[:, np.newaxis], axis=0) - 1.0
    max_corner = np.max(target_pos + target_rad[:, np.newaxis], axis=0) + 1.0
    box_volume = np.prod(max_corner - min_corner)
    if box_volume <= 0:
        return 0.0

    # Generate random sampling points
    points = np.random.uniform(min_corner, max_corner, (n_points, 3))
    covered = np.zeros(n_points, dtype=bool)

    # Mark points covered by target atoms' van der Waals spheres
    for pos, rad in zip(target_pos, target_rad):
        covered |= np.linalg.norm(points - pos, axis=1) < rad

    # Exclude points covered by surrounding atoms (steric hindrance)
    if consider_surrounding:
        valid_env = [a for a in env_atoms if a.coordinates and a.vdw_radius and a.vdw_radius > 0]
        for a in valid_env:
            env_pos = np.array([a.coordinates.x, a.coordinates.y, a.coordinates.z])
            env_rad = a.vdw_radius
            covered &= np.linalg.norm(points - env_pos, axis=1) >= env_rad

    return (covered.sum() / n_points) * box_volume


def parse_group(group_str):
    """Parses elemental composition of a single group, supports ionic forms and complex groups"""
    elem_count = defaultdict(int)
    # Remove ionic charge symbols (e.g., +, -, 2+, etc.)
    group_str = re.sub(r'[+-]\d*$', '', group_str)
    # Remove possible hyphens
    group_str = group_str.replace('-', '')

    # Regex pattern to match element symbols (e.g., O, Cl) and numbers (e.g., 3, 2)
    pattern = r"([A-Z][a-z]?)(\d*)"
    matches = re.findall(pattern, group_str)
    for elem, num in matches:
        if elem:
            count = int(num) if num else 1
            elem_count[elem] += count
    return dict(elem_count)


def parse_disorder_site(disorder_site_str):
    """Parses disorder sites from Excel, supports multiple formats:
    - Space-separated: e.g., "H -OH", "H Cl Br CH3"
    - Slash-separated: e.g., "CH3/CH3", "CH3 Br/Cl I"
    - Ionic forms: e.g., "N+ P+", "Br- Cl-"
    """
    if pd.isna(disorder_site_str) or disorder_site_str.strip() == "":
        return None, None

    # First attempt: split into two groups using slash (/)
    if '/' in disorder_site_str:
        groups = [g.strip() for g in disorder_site_str.split('/') if g.strip()]
        if len(groups) == 2:
            group1_elems = defaultdict(int)
            # Process Group A (may contain multiple space-separated components)
            for part in groups[0].split():
                part_elems = parse_group(part)
                for elem, count in part_elems.items():
                    group1_elems[elem] += count

            group2_elems = defaultdict(int)
            # Process Group B (may contain multiple space-separated components)
            for part in groups[1].split():
                part_elems = parse_group(part)
                for elem, count in part_elems.items():
                    group2_elems[elem] += count

            return dict(group1_elems), dict(group2_elems)

    # Second attempt: split into two groups using spaces (backward compatibility)
    groups = [g.strip() for g in disorder_site_str.split() if g.strip()]
    # Handle mixed single-atom/multi-atom substitution (split into two groups evenly)
    if len(groups) >= 2:
        mid = len(groups) // 2
        group1_parts = groups[:mid]
        group2_parts = groups[mid:]

        group1_elems = defaultdict(int)
        for part in group1_parts:
            part_elems = parse_group(part)
            for elem, count in part_elems.items():
                group1_elems[elem] += count

        group2_elems = defaultdict(int)
        for part in group2_parts:
            part_elems = parse_group(part)
            for elem, count in part_elems.items():
                group2_elems[elem] += count

        return dict(group1_elems), dict(group2_elems)

    return None, None  # Unparseable format


def get_atoms_by_element(atom_list, target_elems):
    """Filters atoms by elemental composition (adapts to element requirements from disorder sites)"""
    if not target_elems:
        return []
    # Count atoms by element
    elem_counter = defaultdict(list)
    for atom in atom_list:
        if atom.atomic_symbol in target_elems and atom.coordinates and atom.occupancy != 1.0:
            elem_counter[atom.atomic_symbol].append(atom)

    # Filter atoms according to element requirements (prioritize atoms with higher occupancy)
    selected_atoms = []
    for elem, req_count in target_elems.items():
        available = elem_counter.get(elem, [])
        if len(available) < req_count:
            return []  # Insufficient atoms, filter failed
        # Sort by occupancy in descending order, take top 'req_count' atoms
        available_sorted = sorted(available, key=lambda x: x.occupancy, reverse=True)
        selected_atoms.extend(available_sorted[:req_count])
    return selected_atoms


def process_single_crystal(crystal_id, disorder_site_str):
    """Processes a single crystal: configures atom groups based on Excel disorder site data and calculates similarity coefficient"""
    # Predefine result fields (includes original disorder site info from Excel)
    result = {
        "crystal_id": crystal_id,
        "disorder_site": disorder_site_str,
        "vdw_volume_group1": None,
        "vdw_volume_group2": None,
        "delta": None,  # Δ: Minimum non-overlapping volume (absolute difference between two volumes)
        "tau": None,    # τ: Maximum overlapping volume (smaller of the two volumes)
        "similarity_coefficient": None,
        "status": "Pending calculation"
    }

    try:
        # 1. Read CSD crystal data
        cry_reader = io.CrystalReader("CSD")
        cry = cry_reader.crystal(crystal_id)
        mol = cry.molecule

        # 2. Extract all atoms (main structure + disordered groups)
        # Main structure atoms (occupancy = 1.0)
        main_atoms = [AtomInfo(atom) for atom in mol.atoms if atom.coordinates]
        # Disordered group atoms (extracted from disorder.assemblies, correct atom occupancy)
        disorder_atoms = []
        if cry.disorder:
            for assembly in cry.disorder.assemblies:
                for group in assembly.groups:
                    group_occupancy = group.occupancy
                    for ccdc_atom in group.atoms:
                        disorder_atoms.append(AtomInfo(ccdc_atom, occupancy=group_occupancy))
        all_atoms = main_atoms + disorder_atoms

        # 3. Parse disorder site and configure two atom groups
        group1_elems, group2_elems = parse_disorder_site(disorder_site_str)
        if not group1_elems or not group2_elems:
            result["status"] = "Failed (invalid disorder site format)"
            return result

        # Filter atoms with occupancy ≠ 1.0 (disordered atoms)
        non_one_atoms = [atom for atom in all_atoms if atom.occupancy is not None and atom.occupancy != 1.0]
        # Filter two groups of atoms according to element requirements from disorder site
        selected_group1 = get_atoms_by_element(non_one_atoms, group1_elems)
        selected_group2 = get_atoms_by_element(non_one_atoms, group2_elems)
        if not selected_group1 or not selected_group2:
            result["status"] = "Failed (insufficient disordered atoms/element mismatch)"
            return result

        # 4. Construct target atom groups (main structure atoms + filtered disordered atoms)
        all_occ1_atoms = [atom for atom in all_atoms if atom.occupancy == 1.0]
        target_group1 = all_occ1_atoms + selected_group1
        target_group2 = all_occ1_atoms + selected_group2

        # 5. Identify environmental atoms (exclude main structure atoms in target groups)
        env_atoms = [a for a in all_occ1_atoms if a not in target_group1 and a not in target_group2]

        # 6. Calculate van der Waals volume (analytical method, balances accuracy and speed)
        vol1 = calculate_analytic_vdw_volume(target_group1)
        vol2 = calculate_analytic_vdw_volume(target_group2)
        if vol1 <= 0 or vol2 <= 0:
            result["status"] = "Failed (abnormal volume calculation, result ≤ 0)"
            return result

        # 7. Calculate similarity coefficient (ε = 1 - Δ/τ)
        delta = abs(vol1 - vol2)
        tau = min(vol1, vol2)
        similarity = round(1 - (delta / tau) if tau != 0 else 0, 4)

        # 8. Populate results
        result["vdw_volume_group1"] = round(vol1, 2)
        result["vdw_volume_group2"] = round(vol2, 2)
        result["delta"] = round(delta, 2)
        result["tau"] = round(tau, 2)
        result["similarity_coefficient"] = similarity
        result["status"] = "Success"

    except Exception as e:
        # Catch all exceptions and retain error information
        error_msg = str(e)[:60] + "..." if len(str(e)) > 60 else str(e)
        result["status"] = f"Failed ({error_msg})"

    return result


def read_excel_crystal_data(excel_path):
    """Reads crystal IDs and disorder sites from true_SS_matched_table.xlsx"""
    if not os.path.exists(excel_path):
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    # Read Excel, retain only "Identifier" (crystal ID) and "disorder site" columns
    df = pd.read_excel(excel_path, sheet_name="Sheet1")
    # Check for required columns
    required_cols = ["Identifier", "disorder site"]
    if not all(col in df.columns for col in required_cols):
        raise ValueError("Excel file missing required columns: must contain 'Identifier' and 'disorder site'")
    # Remove duplicates and null values
    crystal_data = df[["Identifier", "disorder site"]].dropna()
    crystal_data = crystal_data.drop_duplicates(subset=["Identifier"])
    return crystal_data.values.tolist()  # Format: [[crystal_id1, disorder_site1], ...]


def main():
    # 1. Configure file paths (modify according to actual file locations!)
    excel_path = "true_SS_matched_table.xlsx"  # Path to Excel input file
    output_csv = "similarity_coefficient_with_disorder_site.csv"  # Path to output CSV file

    # 2. Read crystal data from Excel
    print("Reading crystal IDs and disorder sites from Excel...")
    try:
        crystal_list = read_excel_crystal_data(excel_path)
        print(f"Successfully read {len(crystal_list)} valid crystal entries")
    except Exception as e:
        print(f"Failed to read Excel: {e}")
        return

    # 3. Batch process each crystal
    print("\nStarting batch calculation of similarity coefficients...")
    results = []
    total = len(crystal_list)
    for idx, (crystal_id, disorder_site) in enumerate(crystal_list, 1):
        print(f"Progress: {idx}/{total} | Crystal ID: {crystal_id} | Disorder Site: {disorder_site}")
        result = process_single_crystal(crystal_id, disorder_site)
        results.append(result)

    # 4. Generate output CSV
    print(f"\nGenerating result file: {output_csv}")
    df_result = pd.DataFrame(results)
    # Reorder columns for readability
    column_order = [
        "crystal_id", "disorder_site", "vdw_volume_group1", "vdw_volume_group2",
        "delta", "tau", "similarity_coefficient", "status"
    ]
    df_result = df_result[column_order]
    # Save CSV (UTF-8-sig supports opening in Chinese Excel)
    df_result.to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"Results saved to: {os.path.abspath(output_csv)}")

    # 5. Output statistical information
    success_count = len([r for r in results if r["status"] == "Success"])
    fail_count = total - success_count
    print(f"\nCalculation completed! Success: {success_count}, Failed: {fail_count}")


if __name__ == "__main__":
    main()