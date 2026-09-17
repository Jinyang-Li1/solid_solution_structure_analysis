import pandas as pd
import os
import csv
import sys
from collections import defaultdict

import ccdc.io
import ccdc.crystal


class RealTimeReporter:
    def __init__(self, detailed_file_path):
        self.detailed_file = detailed_file_path
        self.success_pairs = defaultdict(bool)
        self._init_detailed_file()

    def _init_detailed_file(self):
        with open(self.detailed_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                "group_id", "filename1", "filename2",
                "refcode1", "refcode2",
                "spacegroup1", "a1", "b1", "c1", "alpha1", "beta1", "gamma1",
                "spacegroup2", "a2", "b2", "c2", "alpha2", "beta2", "gamma2",
                "matched_molecules", "RMSD", "status", "note"
            ])

    def _get_cell_params(self, entry):
        try:
            crystal = entry.crystal
            spacegroup = crystal.spacegroup_number_and_setting[0]
            a = crystal.cell_lengths.a
            b = crystal.cell_lengths.b
            c = crystal.cell_lengths.c
            alpha = crystal.cell_angles.alpha
            beta = crystal.cell_angles.beta
            gamma = crystal.cell_angles.gamma
            return (spacegroup, round(a, 4), round(b, 4), round(c, 4),
                    round(alpha, 2), round(beta, 2), round(gamma, 2))
        except Exception:
            return ("failed", "", "", "", "", "", "")

    def add_success(self, group_id, filename1, filename2, entry1, entry2, r1, r2, nmatch, rmsd):
        cell1 = self._get_cell_params(entry1)
        cell2 = self._get_cell_params(entry2)

        status = "isomorphous (15/15)" if nmatch == 15 else f"non-isomorphous ({nmatch}/15)"
        if nmatch == 15:
            self.success_pairs[group_id] = True

        with open(self.detailed_file, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                group_id, filename1, filename2,
                r1, r2,
                cell1[0], cell1[1], cell1[2], cell1[3], cell1[4], cell1[5], cell1[6],
                cell2[0], cell2[1], cell2[2], cell2[3], cell2[4], cell2[5], cell2[6],
                nmatch, round(rmsd, 4), status, ""
            ])
            f.flush()

    def add_failure(self, group_id, filename1, filename2, r1, r2, reason, entry1=None, entry2=None):
        cell1 = ("", "", "", "", "", "", "")
        cell2 = ("", "", "", "", "", "", "")
        if entry1:
            cell1 = self._get_cell_params(entry1)
        if entry2:
            cell2 = self._get_cell_params(entry2)

        with open(self.detailed_file, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                group_id, filename1, filename2,
                r1, r2,
                cell1[0], cell1[1], cell1[2], cell1[3], cell1[4], cell1[5], cell1[6],
                cell2[0], cell2[1], cell2[2], cell2[3], cell2[4], cell2[5], cell2[6],
                "", "", "failed", reason
            ])
            f.flush()

    def get_summary(self, all_group_info):
        summary = []
        for group_id in sorted(all_group_info.keys()):
            filename1, identifier1, filename2, identifier2 = all_group_info[group_id]
            is_iso = "True" if self.success_pairs.get(group_id, False) else "False"
            summary.append((group_id, filename1, identifier1, filename2, identifier2, is_iso))
        return summary


class PairAnalysis:
    def __init__(self, pairs_with_meta, reporter):
        self.pairs_with_meta = pairs_with_meta
        self.reporter = reporter
        if not self.pairs_with_meta:
            raise ValueError("No refcode pairs for analysis")

    def similarity(self, ref_crystal, query_crystal):
        try:
            ps = ccdc.crystal.PackingSimilarity()
            ps.settings.ignore_hydrogen_positions = True
            ps.settings.ignore_hydrogen_counts = False
            ps.settings.ignore_bond_types = False
            ps.settings.ignore_bond_counts = True
            ps.settings.match_entire_packing_shell = False
            ps.settings.allow_molecular_differences = True
            ps.settings.ignore_smallest_components = True
            ps.settings.packing_shell_size = 15
            ps.settings.timeout_ms = 20000
            result = ps.compare(ref_crystal, query_crystal)
            return (result.nmatched_molecules, result.rmsd) if result else (None, "timeout")
        except Exception as e:
            return (None, f"error: {str(e)}")

    def run(self):
        reader = ccdc.io.EntryReader('CSD')
        total_pairs = len(self.pairs_with_meta)
        for idx, (group_id, filename1, filename2, r1, r2) in enumerate(self.pairs_with_meta):
            progress = (idx + 1) / total_pairs * 100
            print(f"[Progress: {progress:.1f}%] Processing {idx + 1}/{total_pairs}: group_id={group_id}, {r1} vs {r2}",
                  end='\r')

            entry1 = None
            entry2 = None
            try:
                entry1 = reader.entry(r1)
                entry2 = reader.entry(r2)
                nmatch, rmsd_or_reason = self.similarity(entry1.crystal, entry2.crystal)
                if nmatch is None:
                    self.reporter.add_failure(
                        group_id, filename1, filename2, r1, r2,
                        reason=rmsd_or_reason, entry1=entry1, entry2=entry2
                    )
                else:
                    self.reporter.add_success(
                        group_id, filename1, filename2, entry1, entry2,
                        r1, r2, nmatch, rmsd_or_reason
                    )

            except Exception as e:
                self.reporter.add_failure(
                    group_id, filename1, filename2, r1, r2,
                    reason=f"refcode processing failed: {str(e)}",
                    entry1=entry1, entry2=entry2
                )
        print()


def main():
    input_csv = "3.smiles_identifier_matched_preprocessed.csv"
    detailed_file = "detailed_pair_results_with_cell_params.csv"
    summary_file = "filename_isomorphism_summary.csv"

    print("Reading input CSV...")
    try:
        encodings = ['gbk', 'utf-8', 'latin-1', 'utf-16']
        for encoding in encodings:
            try:
                df = pd.read_csv(input_csv, encoding=encoding)
                print(f"Successfully loaded CSV with encoding {encoding}")
                break
            except UnicodeDecodeError:
                continue
        else:
            raise Exception("All encoding attempts failed")
    except Exception as e:
        print(f"CSV read error: {e}")
        return

    group_id_col = 0
    filename1_col = 1
    identifier1_col = 3
    filename2_col = 5
    identifier2_col = 7
    id1_col = 3
    id2_col = 7

    pairs_with_meta = []
    all_group_info = {}

    for _, row in df.iterrows():
        group_id = row.iloc[group_id_col]
        filename1 = row.iloc[filename1_col]
        identifier1 = row.iloc[identifier1_col]
        filename2 = row.iloc[filename2_col]
        identifier2 = row.iloc[identifier2_col]

        all_group_info[group_id] = (filename1, identifier1, filename2, identifier2)

        r1_list = [r.strip() for r in str(row.iloc[id1_col]).split(',') if r.strip()]
        r2_list = [r.strip() for r in str(row.iloc[id2_col]).split(',') if r.strip()]

        for r1 in r1_list:
            for r2 in r2_list:
                pairs_with_meta.append((group_id, filename1, filename2, r1, r2))

    print(f"Generated {len(pairs_with_meta)} refcode pairs, covering {len(all_group_info)} groups")

    print("\nStarting packing similarity analysis...")
    reporter = RealTimeReporter(detailed_file_path=detailed_file)
    analysis = PairAnalysis(pairs_with_meta=pairs_with_meta, reporter=reporter)

    try:
        analysis.run()
    except Exception as e:
        print(f"\nAnalysis interrupted: {e}")
        return

    print("\nAnalysis finished, generating summary...")
    summary = reporter.get_summary(all_group_info)
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            "group_id",
            "filename1", "identifier1",
            "filename2", "identifier2",
            "isomorphous"
        ])
        for row in summary:
            writer.writerow(row)

    print("\nOutput files:")
    print(f"   Detailed results with cell parameters: {detailed_file}")
    print(f"   Summary with identifiers: {summary_file}")


if __name__ == "__main__":
    main()