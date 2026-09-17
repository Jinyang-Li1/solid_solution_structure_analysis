# solid_solution_structure_analysis
Code for CSD searches and solid solution structure analysis in 'Revisiting Kitaigorodskii's Rules via Crystal Structure Database Mining: Every Criterion Has Violations, and the Integrity of Key Interactions Governs Molecular Solid-Solution Formation' by Li et al.

Code Files:

search_in_csd.py performs customised searches within the Cambridge Structural Database to retrieve solid-solution related crystal entries.
classify_crystal_disorder.py classifies types of crystallographic disorder for extracted CSD structures.
deduplicate_disordered_structures.py removes duplicate entries from disordered crystal datasets.
extract_disorder_details.py extracts detailed crystallographic disorder information from CSD entries.
fetch_solvent_synonyms.py retrieves synonym names for crystalline solvents for subsequent filtering.
filter_solvent_entries.py filters out structures dominated by solvent disorder.
filter_disorder_occupancy.py filters crystal structures according to disorder site occupancy values.
disordered_details_split_processor.py enables batch processing and splitting of disorder metadata.
disordered_molecule_analysis.py analyses molecular geometry and disorder features of disordered crystals.
disordered_atom_SMILES_analysis.py generates and analyses SMILES strings for disordered atomic fragments.
Csd_formula_solid_solution_analysis.py screens and validates candidate solid-solution structures by chemical formula matching.
Disorder_analysis_with_diagrams.py computes disorder statistics and exports visualisation plots.
Dataset structure trace.py traces the provenance and filtering history of the curated solid-solution dataset.
Intention verification.py verifies the structural identity and compositional characteristics of candidate solid-solution entries.
similarity_coefficient_calculation.py calculates similarity coefficients for solid-solution components.
Crystal packing similarity.py calculates crystal packing similarity and evaluates isostructurality, and exports pairwise comparison results into .csv files
