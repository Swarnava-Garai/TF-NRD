#!/usr/bin/env python3
"""
TF-NRD Interface Analysis Suite
================================
Calculates structural interface properties (BSA, FNP, FBU, LD) for:
  - Protein-DNA complexes (PNA DNA)
  - Protein-RNA complexes (PNA RNA)
  - Dual RNA/DNA binding protein complexes (PNA DRBP)
  - Protein-Protein complexes (PP)

Processes datasets matching TF_interface_13.03.ipynb logic, preserves all metadata
from TFNRDv1.0_data_final.xlsx, and exports publication-ready multi-sheet Excel workbooks
and CSV summary files.

Author: Antigravity Team / PhD Project Suite
"""

import argparse
import logging
from math import sqrt
import os
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("TFInterfaceAnalyzer")

DEFAULT_DATA_EXCEL = Path("/home/labuser/Projects/PhD_projects/swarnava_TF_work/Interface/13.03.2026/TFNRDv1.0_data_final.xlsx")
DEFAULT_BASE_DIR = Path("/home/labuser/Projects/PhD_projects/swarnava_TF_work/Interface/13.03.2026")
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent / "results" / "Interface"


def calculate_interface_metrics(filepath: str | Path) -> dict:
    """
    Calculates interface properties from a single .int file:
      - BSA: Buried Surface Area (total)
      - FNP: Fraction of Non-Polar BSA (%)
      - FBU: Fraction of Buried Residues (%)
      - LD: Local Density (average contact count within 12 Å)
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return {'bsa': 0.0, 'fnp': 0.0, 'fbu': 0.0, 'ld': 0.0}

    total_bsa = polar_bsa = 0.0
    count = buried_count = 0.0
    coordinates = []

    with open(filepath, "r", encoding="utf-8", errors="replace") as intf:
        for line in intf:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) < 10:
                    continue
                atom = parts[2]
                try:
                    x = float(parts[-6])
                    y = float(parts[-5])
                    z = float(parts[-4])
                    sasa_free = float(parts[-3])
                    sasa_complx = float(parts[-2])
                    bsa = float(parts[-1])
                except (ValueError, IndexError):
                    continue

                coordinates.append([x, y, z])
                total_bsa += bsa
                count += 1
                if atom[0] == "C":
                    polar_bsa += bsa
                if sasa_complx == 0.0:
                    buried_count += 1

    if count == 0:
        return {'bsa': 0.0, 'fnp': 0.0, 'fbu': 0.0, 'ld': 0.0}

    ld = 0.0
    coords_arr = np.array(coordinates)
    if len(coords_arr) > 1:
        dists = np.linalg.norm(coords_arr[:, np.newaxis, :] - coords_arr[np.newaxis, :, :], axis=-1)
        np.fill_diagonal(dists, np.inf)
        ld = float(np.sum(dists <= 12.0))
        ld = round(ld / count, 2)

    fnp = round(polar_bsa / total_bsa * 100, 2) if total_bsa > 0 else 0.0
    fbu = round(buried_count / count * 100, 2)

    return {
        'bsa': round(total_bsa, 2),
        'fnp': fnp,
        'fbu': fbu,
        'ld': ld
    }


def process_pna_subset(
    df_subset: pd.DataFrame,
    pna_dir: Path,
    category: str = "DNA"
) -> pd.DataFrame:
    """
    Processes a Protein-Nucleic Acid (PNA) dataset subset ('DNA', 'RNA', or 'DRBP')
    while preserving all original metadata fields from TFNRDv1.0_data_final.xlsx.
    """
    results = []
    pna_dir = Path(pna_dir)
    partner_col = 'RNA' if category == 'RNA' else 'DNA'
    partner_ext = 'R' if category == 'RNA' else 'D'

    # Filter out existing raw BSA columns to replace with calculated metrics
    drop_raw_cols = ['Complex', 'Protein.1', 'Neucleic acid']
    meta_cols = [c for c in df_subset.columns if c not in drop_raw_cols]

    for index, row in df_subset.iterrows():
        entry = {col: row[col] for col in meta_cols}
        pdb_id = str(row['PDB_ID']).strip()
        pro_chan = str(row['Protein']).strip()
        partner_chan = str(row[partner_col]).strip() if partner_col in row and pd.notna(row[partner_col]) else ""

        dir_name = f"{pdb_id}_{pro_chan}_{partner_chan}"
        target_dir = pna_dir / dir_name

        complex_file = target_dir / f"{pdb_id}.int"
        protein_file = target_dir / f"{pdb_id}P.int"
        partner_file = target_dir / f"{pdb_id}{partner_ext}.int"

        m_c = calculate_interface_metrics(complex_file)
        m_p = calculate_interface_metrics(protein_file)
        m_partner = calculate_interface_metrics(partner_file)

        entry.update({
            "Category": category,

            "BSA_complex": m_c["bsa"],
            "BSA_protein": m_p["bsa"],
            f"BSA_{category}": m_partner["bsa"],

            "FNP_complex": m_c["fnp"],
            "FNP_protein": m_p["fnp"],
            f"FNP_{category}": m_partner["fnp"],

            "FBU_complex": m_c["fbu"],
            "FBU_protein": m_p["fbu"],
            f"FBU_{category}": m_partner["fbu"],

            "LD_complex": m_c["ld"],
            "LD_protein": m_p["ld"],
            f"LD_{category}": m_partner["ld"],
        })
        results.append(entry)

    return pd.DataFrame(results)


def process_protein_protein_dataset(
    df_pp: pd.DataFrame,
    pp_dir: Path
) -> pd.DataFrame:
    """
    Processes Protein-Protein interface dataset entries while preserving all original metadata fields.
    Handles chain-specifically named files ({pdb_id}_{c1}{c2}.int) or standard ({pdb_id}.int).
    """
    results = []
    pp_dir = Path(pp_dir)

    drop_raw_cols = ['BSA complex', 'Protein 1', 'Bprotein 2']
    meta_cols = [c for c in df_pp.columns if c not in drop_raw_cols]

    for index, row in df_pp.iterrows():
        entry = {col: row[col] for col in meta_cols}
        pdb_id = str(row['PDB ID']).strip()
        c1 = str(row['Chain 1']).split(":")[0].strip() if ":" in str(row['Chain 1']) else str(row['Chain 1']).strip()
        c2 = str(row['Chain 2']).split(":")[1].strip() if ":" in str(row['Chain 2']) else str(row['Chain 2']).strip()

        dir_name = f"{pdb_id}_{c1}_{c2}"
        target_dir = pp_dir / dir_name

        cmplx_f = target_dir / f"{pdb_id}_{c1}{c2}.int"
        c1_f = target_dir / f"{pdb_id}_{c1}.int"
        c2_f = target_dir / f"{pdb_id}_{c2}.int"

        if not (cmplx_f.exists() and c1_f.exists() and c2_f.exists()):
            cmplx_f = target_dir / f"{pdb_id}.int"
            c1_f = target_dir / f"{pdb_id}P.int"
            c2_f = target_dir / f"{pdb_id}D.int"

        m_c = calculate_interface_metrics(cmplx_f)
        m_1 = calculate_interface_metrics(c1_f)
        m_2 = calculate_interface_metrics(c2_f)

        entry.update({
            "Category": "Protein-Protein",

            "BSA_complex": m_c["bsa"],
            "BSA_chain1": m_1["bsa"],
            "BSA_chain2": m_2["bsa"],

            "FNP_complex": m_c["fnp"],
            "FNP_chain1": m_1["fnp"],
            "FNP_chain2": m_2["fnp"],

            "FBU_complex": m_c["fbu"],
            "FBU_chain1": m_1["fbu"],
            "FBU_chain2": m_2["fbu"],

            "LD_complex": m_c["ld"],
            "LD_chain1": m_1["ld"],
            "LD_chain2": m_2["ld"],
        })
        results.append(entry)

    return pd.DataFrame(results)


def generate_summary_statistics(dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Generates summary statistics (mean, std, min, max) for BSA, FNP, FBU, LD across categories.
    """
    summary_rows = []
    for cat_name, df in dfs.items():
        if df.empty:
            continue
        metric_cols = [c for c in df.columns if c.startswith(('BSA_', 'FNP_', 'FBU_', 'LD_'))]
        for col in metric_cols:
            summary_rows.append({
                "Dataset": cat_name,
                "Metric": col,
                "Count": len(df[col].dropna()),
                "Mean": round(df[col].mean(), 2),
                "Std": round(df[col].std(), 2),
                "Min": round(df[col].min(), 2),
                "Max": round(df[col].max(), 2),
            })
    return pd.DataFrame(summary_rows)


class TFInterfaceAnalyzer:
    """
    Interface Analysis Manager for TF-NRD datasets.
    """

    def __init__(
        self,
        data_excel: str | Path = DEFAULT_DATA_EXCEL,
        base_dir: str | Path = DEFAULT_BASE_DIR,
        output_dir: str | Path = DEFAULT_OUTPUT_DIR
    ):
        self.data_excel = Path(data_excel)
        self.base_dir = Path(base_dir)
        self.pna_dir = self.base_dir / "PNA"
        self.pp_dir = self.base_dir / "PP"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run_all(self) -> dict[str, pd.DataFrame]:
        """
        Runs interface property analysis for DNA, RNA, DRBP, and Protein-Protein datasets,
        preserving metadata and exporting separate Excel sheets for PNA and PP.
        """
        if not self.data_excel.exists():
            logger.error(f"Excel file not found: {self.data_excel}")
            return {}

        logger.info(f"Loading TF-NRD Excel dataset: {self.data_excel}")
        xls = pd.ExcelFile(self.data_excel, engine='openpyxl')

        # 1. Unique_Interface Sheet (PNA)
        if 'Unique_Interface' in xls.sheet_names:
            tfnrd_pna = pd.read_excel(xls, sheet_name='Unique_Interface', engine='openpyxl')
            tfnrd_pna = tfnrd_pna[~tfnrd_pna['PDB_ID'].astype(str).str.startswith('#')].dropna(subset=['PDB_ID'])
        else:
            logger.error("Sheet 'Unique_Interface' not found in Excel file.")
            tfnrd_pna = pd.DataFrame()

        # 2. pro_pro_interface Sheet (PP)
        if 'pro_pro_interface' in xls.sheet_names:
            tfnrd_pp = pd.read_excel(xls, sheet_name='pro_pro_interface', engine='openpyxl')
            tfnrd_pp = tfnrd_pp[~tfnrd_pp['PDB ID'].astype(str).str.startswith('#')].dropna(subset=['PDB ID'])
        else:
            logger.error("Sheet 'pro_pro_interface' not found in Excel file.")
            tfnrd_pp = pd.DataFrame()

        output_dfs = {}

        # Process PNA DNA (Rows 0..438)
        if not tfnrd_pna.empty:
            logger.info("--- Processing PNA DNA Interface Dataset (438 entries) ---")
            dna_subset = tfnrd_pna.iloc[:438]
            output_dfs['PNA_DNA'] = process_pna_subset(dna_subset, self.pna_dir, category="DNA")

            # Process PNA RNA (Rows 438..455)
            logger.info("--- Processing PNA RNA Interface Dataset (17 entries) ---")
            rna_subset = tfnrd_pna.iloc[438:455]
            output_dfs['PNA_RNA'] = process_pna_subset(rna_subset, self.pna_dir, category="RNA")

            # Process PNA DRBP (Rows 455+)
            logger.info("--- Processing PNA DRBP Interface Dataset (6 entries) ---")
            drbp_subset = tfnrd_pna.iloc[455:]
            output_dfs['PNA_DRBP'] = process_pna_subset(drbp_subset, self.pna_dir, category="DRBP")

            # Combined PNA Dataset
            output_dfs['PNA'] = pd.concat([output_dfs['PNA_DNA'], output_dfs['PNA_RNA'], output_dfs['PNA_DRBP']], ignore_index=True)

        # Process Protein-Protein (PP) Dataset
        if not tfnrd_pp.empty:
            logger.info("--- Processing Protein-Protein (PP) Interface Dataset (361 entries) ---")
            output_dfs['PP'] = process_protein_protein_dataset(tfnrd_pp, self.pp_dir)

        # Summary statistics
        output_dfs['Summary'] = generate_summary_statistics(output_dfs)

        # Export Multi-Sheet Excel File
        excel_out_path = self.output_dir / "TF_Interface_Properties.xlsx"
        logger.info(f"Exporting multi-sheet Excel file to: {excel_out_path}")

        with pd.ExcelWriter(excel_out_path, engine='openpyxl') as writer:
            if 'PNA' in output_dfs:
                output_dfs['PNA'].to_excel(writer, sheet_name='PNA', index=False)
            if 'PP' in output_dfs:
                output_dfs['PP'].to_excel(writer, sheet_name='PP', index=False)
            if 'PNA_DNA' in output_dfs:
                output_dfs['PNA_DNA'].to_excel(writer, sheet_name='PNA_DNA', index=False)
            if 'PNA_RNA' in output_dfs:
                output_dfs['PNA_RNA'].to_excel(writer, sheet_name='PNA_RNA', index=False)
            if 'PNA_DRBP' in output_dfs:
                output_dfs['PNA_DRBP'].to_excel(writer, sheet_name='PNA_DRBP', index=False)
            if 'Summary' in output_dfs:
                output_dfs['Summary'].to_excel(writer, sheet_name='Summary', index=False)

        # Export Individual CSV files for backward compatibility
        if 'PNA' in output_dfs:
            output_dfs['PNA'].to_csv(self.output_dir / "PNA_interface_properties.csv", index=False)
        if 'PP' in output_dfs:
            output_dfs['PP'].to_csv(self.output_dir / "PP_interface_properties.csv", index=False)
        
        all_dfs = [output_dfs[k] for k in ['PNA', 'PP'] if k in output_dfs and not output_dfs[k].empty]
        if all_dfs:
            merged_all = pd.concat(all_dfs, ignore_index=True)
            merged_all.to_csv(self.output_dir / "merged_interface_properties.csv", index=False)

        logger.info("TF-NRD Interface Analysis Pipeline Completed Successfully!")
        return output_dfs


def main():
    parser = argparse.ArgumentParser(description="TF-NRD Interface Analysis Suite")
    parser.add_argument("--excel", default=DEFAULT_DATA_EXCEL, help="Path to TF-NRD data excel file")
    parser.add_argument("--base-dir", default=DEFAULT_BASE_DIR, help="Base directory containing PNA and PP data folders")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Output directory to save Excel/CSV results")

    args = parser.parse_args()
    analyzer = TFInterfaceAnalyzer(data_excel=args.excel, base_dir=args.base_dir, output_dir=args.output_dir)
    analyzer.run_all()


if __name__ == "__main__":
    main()