"""
disordered_atom_SIMILES_analysis.py - Analyzes disordered atoms from CSD refcodes to generate SMILES,
classifies refcodes by SMILES similarity, and exports structured results.
"""

import pandas as pd
import numpy as np
from ccdc import io
from collections import defaultdict, Counter
from rdkit import Chem
from rdkit.Chem import AllChem

# Output file configuration
SAME_DETAIL_FILE = "10. same smiles segments.csv"
DIFF_DETAIL_FILE = "10. different smiles segments.csv"
SAME_LIST_FILE = "10. same smiles refcodes.csv"
DIFF_LIST_FILE = "10. different smiles refcodes.csv"

# Constant dictionaries (replace small functions)
ISOTOPE_MAP = {'D': 'H', 'T': 'H'}
COVALENT_RADII = {'H': 0.31, 'C': 0.73, 'N': 0.75, 'O': 0.73, 'S': 1.05,
                  'Cl': 0.99, 'Br': 1.14, 'I': 1.33, 'F': 0.71, 'P': 1.07}
MAX_VALENCE = {'H': 1, 'C': 4, 'N': 5, 'O': 2, 'S': 6,
               'Cl': 1, 'Br': 1, 'I': 1, 'F': 1, 'P': 5}


class AtomInfo:
    """Encapsulates atom information with isotope handling and validity checks."""
    def __init__(self, ccdc_atom, occupancy=None):
        self.original_label = getattr(ccdc_atom, 'label', f"Atom_{id(ccdc_atom)}")
        self.original_elem = getattr(ccdc_atom, 'atomic_symbol', 'X').strip()
        self.elem_for_smiles = ISOTOPE_MAP.get(self.original_elem.upper(), self.original_elem)
        self.coords = getattr(ccdc_atom, 'coordinates', None)
        self.occupancy_for_filter = self._validate_occupancy(occupancy or getattr(ccdc_atom, 'occupancy', None))

    def _validate_occupancy(self, occupancy):
        """Validate occupancy (0 < occ < 1)"""
        try:
            occ = float(occupancy)
            return occ if 0 < occ < 1.0 else None
        except (TypeError, ValueError):
            return None

    def is_valid(self):
        """Check if atom is valid for analysis"""
        try:
            return (Chem.GetPeriodicTable().GetAtomicNumber(self.elem_for_smiles) > 0 and
                    self.coords and all(hasattr(self.coords, attr) for attr in ['x', 'y', 'z']) and
                    self.occupancy_for_filter is not None)
        except:
            return False

    def get_coords(self):
        """Return rounded 3D coordinates"""
        return (round(self.coords.x, 3), round(self.coords.y, 3), round(self.coords.z, 3)) if self.coords else None


def load_crystal_and_entry(refcode):
    """Load CSD Crystal (for SMILES) and Entry (for disorder details)"""
    try:
        with io.CrystalReader('CSD') as csd_reader:
            return csd_reader.crystal(refcode), csd_reader.entry(refcode), None
    except Exception as e:
        return None, None, str(e)


def extract_disordered_atoms(crystal):
    """Extract valid disordered atoms from crystal"""
    disordered_atoms = []
    if not crystal:
        return disordered_atoms

    # Extract from main molecule
    if hasattr(crystal, 'molecule'):
        disordered_atoms.extend([AtomInfo(a) for a in crystal.molecule.atoms if AtomInfo(a).is_valid()])

    # Extract from disorder assemblies
    if hasattr(crystal, 'disorder') and crystal.disorder:
        try:
            for assembly in crystal.disorder.assemblies:
                for group in assembly.groups:
                    disordered_atoms.extend([AtomInfo(a, group.occupancy) for a in group.atoms if AtomInfo(a, group.occupancy).is_valid()])
        except Exception as e:
            print(f"Disorder group extract failed: {e}")
    return disordered_atoms


def extract_raw_disorder_details(entry):
    """Extract raw disorder details from CSD Entry"""
    try:
        return str(entry.disorder_details).strip() if entry else ""
    except Exception as e:
        print(f"Disorder details extract failed: {e}")
        return ""


def group_by_occupancy(atoms):
    """Group atoms by rounded occupancy (3 decimals)"""
    groups = defaultdict(list)
    for atom in atoms:
        groups[round(atom.occupancy_for_filter, 3)].append(atom)
    return dict(groups)


def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms"""
    c1, c2 = atom1.get_coords(), atom2.get_coords()
    if not c1 or not c2:
        return float('inf')
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def create_rdkit_molecule(atoms):
    """Create RDKit molecule with bonds and conformer"""
    if not atoms:
        return None

    mol = Chem.RWMol()
    atom_map = {}  # orig_idx -> rdkit_idx
    tolerance = 0.25 if len(atoms) <= 10 else 0.2

    # Add atoms
    for idx, atom in enumerate(atoms):
        try:
            atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(atom.elem_for_smiles)
            atom_map[idx] = mol.AddAtom(Chem.Atom(atomic_num))
        except Exception as e:
            print(f"Skip atom {atom.original_label}: {e}")

    if not atom_map:
        return None

    # Add bonds
    bond_count = {idx: 0 for idx in atom_map}
    for i in atom_map:
        for j in atom_map:
            if i >= j:
                continue
            dist = calculate_distance(atoms[i], atoms[j])
            threshold = COVALENT_RADII.get(atoms[i].elem_for_smiles, 0.8) + COVALENT_RADII.get(atoms[j].elem_for_smiles, 0.8) + tolerance
            if 0.01 <= dist <= threshold and bond_count[i] < MAX_VALENCE.get(atoms[i].elem_for_smiles, 4) and bond_count[j] < MAX_VALENCE.get(atoms[j].elem_for_smiles, 4):
                mol.AddBond(atom_map[i], atom_map[j], Chem.BondType.SINGLE)
                bond_count[i] += 1
                bond_count[j] += 1

    # Add conformer
    mol = mol.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for orig_idx, rd_idx in atom_map.items():
        coords = atoms[orig_idx].get_coords()
        if coords:
            conf.SetAtomPosition(rd_idx, coords)
    mol.AddConformer(conf)
    return mol


def get_smiles(mol, atoms):
    """Generate SMILES (fallback to element string if failed)"""
    if not mol:
        return '.'.join(a.original_elem for a in atoms)

    # Handle large single-element groups
    elem_count = Counter(a.original_elem for a in atoms)
    if len(elem_count) == 1 and max(elem_count.values()) > 10:
        return '.'.join([next(iter(elem_count))] * max(elem_count.values()))

    try:
        AllChem.Compute2DCoords(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except:
        return '.'.join(a.original_elem for a in atoms)


def process_refcode(refcode):
    """Process single refcode and return results"""
    result = {'refcode': refcode, 'groups': [], 'all_smiles': [], 'unique_smiles': set(), 'raw_disorder_details': ""}
    print(f"\n=== Processing {refcode} ===")

    # Load data
    crystal, entry, load_err = load_crystal_and_entry(refcode)
    if load_err:
        print(f"Load failed: {load_err}")
        return result

    # Extract disorder details and atoms
    result['raw_disorder_details'] = extract_raw_disorder_details(entry)
    print(f"Disorder details length: {len(result['raw_disorder_details'])}")
    disordered_atoms = extract_disordered_atoms(crystal)

    if not disordered_atoms:
        print("No valid disordered atoms")
        return result
    print(f"Found {len(disordered_atoms)} valid atoms")

    # Process occupancy groups
    for occ, atoms in group_by_occupancy(disordered_atoms).items():
        print(f"  Occ {occ}: {len(atoms)} atoms | {[f'{a.original_elem}({a.original_label})' for a in atoms]}")
        mol = create_rdkit_molecule(atoms)
        smiles = get_smiles(mol, atoms)
        print(f"SMILES: {smiles}")

        # Collect atom details
        result['groups'].extend([{
            'refcode': refcode, 'occupancy': round(occ, 3),
            'original_atom_label': a.original_label, 'original_element': a.original_elem,
            'x_coord': a.get_coords()[0] if a.get_coords() else '',
            'y_coord': a.get_coords()[1] if a.get_coords() else '',
            'z_coord': a.get_coords()[2] if a.get_coords() else '',
            'group_atom_count': len(atoms), 'group_smiles': smiles
        } for a in atoms])
        result['all_smiles'].append(smiles)
        result['unique_smiles'].add(smiles)

    return result


def classify_and_save(results):
    """Classify results and export to CSV"""
    same_seg, diff_seg, same_ref, diff_ref = [], [], [], []

    for res in results:
        if not res['all_smiles'] or not res['groups']:
            continue

        # Base refcode row
        base_row = {
            'refcode': res['refcode'], 'num_groups': len(res['all_smiles']),
            'num_unique_smiles': len(res['unique_smiles']),
            'all_smiles': ' | '.join(res['all_smiles']),
            'unique_smiles': ' | '.join(res['unique_smiles'])
        }

        # Classify
        if len(res['unique_smiles']) == 1:
            same_seg.extend(res['groups'])
            same_ref.append(base_row)
        else:
            diff_seg.extend(res['groups'])
            diff_ref.append({**base_row, 'raw_disorder_details': res['raw_disorder_details']})

    # Export
    if same_seg:
        pd.DataFrame(same_seg).to_csv(SAME_DETAIL_FILE, index=False, encoding='utf-8-sig')
    if same_ref:
        pd.DataFrame(same_ref).to_csv(SAME_LIST_FILE, index=False, encoding='utf-8-sig')
    if diff_seg:
        pd.DataFrame(diff_seg).to_csv(DIFF_DETAIL_FILE, index=False, encoding='utf-8-sig')
    if diff_ref:
        pd.DataFrame(diff_ref)[['refcode', 'num_groups', 'num_unique_smiles', 'all_smiles', 'unique_smiles', 'raw_disorder_details']].to_csv(DIFF_LIST_FILE, index=False, encoding='utf-8-sig')

    # Summary
    print(f"\n=== Summary ===")
    print(f"Total valid refcodes: {len(same_ref) + len(diff_ref)}")
    print(f"Same SMILES: {len(same_ref)} | Different SMILES: {len(diff_ref)}")


def main(input_csv):
    """Main workflow: load refcodes, process, classify, export"""
    try:
        # Load refcodes (robust handling)
        try:
            refcodes = pd.read_csv(input_csv, encoding='utf-8-sig')[pd.read_csv(input_csv, encoding='utf-8-sig').columns[0]].dropna().astype(str).str.strip().tolist()
        except:
            with open(input_csv, 'r', encoding='utf-8-sig') as f:
                refcodes = [l.split(',')[0].strip() for l in f if l.strip() and not (l.strip().lower() in ['refcode', 'id', 'identifier'])]

        refcodes = list(set(filter(None, refcodes)))
        if not refcodes:
            raise ValueError("No valid refcodes found")

        # Process refcodes
        print(f"=== Starting with {len(refcodes)} unique refcodes ===")
        all_results = [process_refcode(rc) for idx, rc in enumerate(refcodes, 1) if print(f"\n{idx}/{len(refcodes)} - Processing {rc}") or True]

        # Classify and export
        classify_and_save(all_results)
        print("\n=== Done ===")

    except Exception as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    main("9. refcodes all same.csv")