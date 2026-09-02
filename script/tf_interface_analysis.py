#!/home/labuser/anaconda3/bin/python
"""
tf_interface_analysis.py
-------------------------
Unified TF-NRD Structural Interface Analysis Suite & ML Dataset Generator.

Consolidates:
  1. Older runs (swarnava_TF_work/Interface/13.03.2026 for X-ray & NMR structures)
  2. Revision runs (swarnava_TF_work/Interface/Revision for Cryo-EM structures)
  3. Feature extraction (BSA, FNP, FBU, LD) across all interface types (Protein-DNA, Protein-RNA, Protein-DRBP, Protein-Protein)
  4. Publication Supplementary Excel Workbooks (Table S1.A, Table S1.B, Table S1.C, Summary_Statistics)
  5. Machine-Learning ready consolidated datasets (merged_interface_features_ml.csv & merged_interface_features_ml.json)

Usage:
  python tf_interface_analysis.py
  python tf_interface_analysis.py --output_dir /home/labuser/Projects/PhD_projects/TF-NRD/Supplementary
"""

import os
import sys
import json
import csv
import argparse
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
PROJECT_ROOT = SCRIPT_DIR.parent
LOG_FILE = SCRIPT_PATH.with_suffix(".log")

# Configure logging to stdout and script-named log file
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_FILE, mode="a", encoding="utf-8")
    ]
)
logger = logging.getLogger(SCRIPT_PATH.stem)

# Default Global Directory Paths
INTERFACE_BASE_DIR = PROJECT_ROOT.parent / "swarnava_TF_work" / "Interface"
OLD_BASE_DIR = INTERFACE_BASE_DIR / "13.03.2026"
REVISION_BASE_DIR = INTERFACE_BASE_DIR / "Revision"

DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Supplementary"

# Standard Residue Sets & Backbone Identifier Atoms
AA_LIST = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEC', 'PYL'
}
RNA_LIST = {'A', 'U', 'C', 'G', 'RA', 'RU', 'RC', 'RG', 'I'}
DNA_LIST = {'DA', 'DT', 'DC', 'DG', 'DI', 'DU', 'T'}


# ==============================================================================
# Helper Functions & Parsers
# ==============================================================================

def clean_str(val) -> str:
    """Strips whitespace and normalizes empty/placeholder strings to ''."""
    if pd.isna(val) or val is None:
        return ""
    s = str(val).strip()
    if s.lower() in ("", "-", "nan", "none", "null", "na"):
        return ""
    return s


def clean_exp_method(val, default_method="X-RAY DIFFRACTION") -> str:
    """Normalizes shorthand experimental methods ('X' -> 'X-RAY DIFFRACTION', 'N' -> 'SOLUTION NMR', etc.)."""
    s = clean_str(val).upper()
    if s == 'X':
        return 'X-RAY DIFFRACTION'
    elif s == 'N':
        return 'SOLUTION NMR'
    elif s == 'ELECTRON MICROSCOPY' or 'CRYO' in s or 'EM' in s:
        return 'ELECTRON MICROSCOPY'
    elif s:
        return s
    return default_method


def calculate_interface_metrics(filepath: str | Path) -> dict:
    """
    Calculates interface properties from a single .int file:
      - BSA: Buried Surface Area (Å²)
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
        tree = cKDTree(coords_arr)
        neighbors = tree.query_ball_tree(tree, r=12.0)
        total_neighbors = sum(len(n) - 1 for n in neighbors)
        ld = round(total_neighbors / count, 2)

    fnp = round(polar_bsa / total_bsa * 100, 2) if total_bsa > 0 else 0.0
    fbu = round(buried_count / count * 100, 2)

    return {
        'bsa': round(total_bsa, 2),
        'fnp': fnp,
        'fbu': fbu,
        'ld': ld
    }


def parse_int_atom_chains(filepath: Path) -> tuple[set, set, set, set]:
    """Parses ATOM lines with BSA > 0 from .int file to extract active chain IDs."""
    p_chains, d_chains, r_chains, other_chains = set(), set(), set(), set()
    if not Path(filepath).exists():
        return p_chains, d_chains, r_chains, other_chains

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('ATOM'):
                parts = line.split()
                if len(parts) >= 10:
                    try:
                        bsa = float(parts[-1])
                        if bsa > 0.0:
                            raw_ch = parts[4]
                            resn = parts[3].strip()
                            base_ch = raw_ch[0]

                            if resn in AA_LIST:
                                p_chains.add(base_ch)
                            elif resn in DNA_LIST:
                                d_chains.add(base_ch)
                            elif resn in RNA_LIST:
                                r_chains.add(base_ch)
                            else:
                                other_chains.add(base_ch)
                    except (ValueError, IndexError):
                        pass
    return p_chains, d_chains, r_chains, other_chains


def determine_interface_formed(target_dir: Path, pdb_id: str, meta: dict) -> str:
    """Determines active interface chain string e.g. 'IJLM:OP' or 'NONE'."""
    complex_f = target_dir / f"{pdb_id}.int"
    if not complex_f.exists():
        c1 = clean_str(meta.get('Chain_1', meta.get('Chain 1', meta.get('Protein_1', meta.get('Subunit1_Chains', '')))))
        c2 = clean_str(meta.get('Chain_2', meta.get('Chain 2', meta.get('Protein_2', meta.get('Subunit2_Chains', '')))))
        if c1 and c2:
            alt_cmplx = target_dir / f"{pdb_id}_{c1}{c2}.int"
            if alt_cmplx.exists():
                complex_f = alt_cmplx

    if not complex_f.exists():
        return "NONE"

    p_ch, d_ch, r_ch, oth_ch = parse_int_atom_chains(complex_f)
    int_type = str(meta.get('Interface_Type', '')).strip()

    if int_type == 'Protein-Protein':
        c1_str = clean_str(meta.get('Chain_1', meta.get('Chain 1', meta.get('Protein_1', ''))))
        c2_str = clean_str(meta.get('Chain_2', meta.get('Chain 2', meta.get('Protein_2', ''))))
        s1_active = [c for c in c1_str if c in p_ch or c in oth_ch]
        s2_active = [c for c in c2_str if c in p_ch or c in oth_ch]

        if not s1_active and not s2_active and len(p_ch) >= 2:
            sorted_p = sorted(p_ch)
            s1_active = [sorted_p[0]]
            s2_active = sorted_p[1:]

        s1 = "".join(sorted(set(s1_active)))
        s2 = "".join(sorted(set(s2_active)))
        if s1 and s2:
            return f"{s1}:{s2}"
        return "NONE"

    p_str = "".join(sorted(p_ch))
    na_ch = sorted(d_ch | r_ch)
    na_str = "".join(na_ch)

    if p_str and na_str:
        return f"{p_str}:{na_str}"
    elif p_str:
        return f"{p_str}:NONE"
    return "NONE"


# ==============================================================================
# Older Runs Processing Suite (swarnava_TF_work/Interface/13.03.2026)
# ==============================================================================

def process_old_runs(old_base_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parses older 13.03.2026 runs (X-ray & NMR structures) for PNA and PP datasets.
    Returns: (df_batch, df_unique_na, df_unique_pp)
    """
    old_base_dir = Path(old_base_dir)
    excel_path = old_base_dir / "TFNRDv1.0_data_final.xlsx"

    pna_dir = old_base_dir / "PNA"
    pp_dir = old_base_dir / "PP"

    records_batch = []
    records_unique_na = []
    records_unique_pp = []

    if not excel_path.exists():
        logger.warning(f"Old excel file not found: {excel_path}")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    xls = pd.ExcelFile(excel_path, engine='openpyxl')

    # 1. PNA (Unique_Interface Sheet)
    if 'Unique_Interface' in xls.sheet_names:
        df_pna_meta = pd.read_excel(xls, sheet_name='Unique_Interface', engine='openpyxl')
        df_pna_meta = df_pna_meta[~df_pna_meta['PDB_ID'].astype(str).str.startswith('#')].dropna(subset=['PDB_ID'])

        for idx, row in df_pna_meta.iterrows():
            pdb_id = clean_str(row.get('PDB_ID')).upper()
            pro_chan = clean_str(row.get('Protein'))
            dna_chan = clean_str(row.get('DNA'))
            rna_chan = clean_str(row.get('RNA'))
            partner_chan = rna_chan if rna_chan else dna_chan

            dir_name = f"{pdb_id}_{pro_chan}_{partner_chan}"
            target_dir = pna_dir / dir_name

            cmplx_f = target_dir / f"{pdb_id}.int"
            prot_f = target_dir / f"{pdb_id}P.int"
            partner_f = target_dir / f"{pdb_id}{'R' if rna_chan else 'D'}.int"

            m_c = calculate_interface_metrics(cmplx_f)
            m_p = calculate_interface_metrics(prot_f)
            m_part = calculate_interface_metrics(partner_f)

            method = clean_exp_method(row.get('Experimental Methods'), default_method="X-RAY DIFFRACTION")
            int_type = "Protein-RNA" if rna_chan else ("Protein-DNA-RNA" if (dna_chan and rna_chan) else "Protein-DNA")

            entry = {
                'PDB_ID': pdb_id,
                'Source_Dataset': '13.03.2026',
                'Exp_Method': method,
                'Protein_Name': clean_str(row.get('Protein_Name', row.get('protein_name', ''))),
                'Source_Organism': clean_str(row.get('Source Organism', '')),
                'Protein_Length': clean_str(row.get('Protein Length', row.get('protein_length', ''))),
                'Resolution': row.get('Resolution (Å)', ''),
                'Oligomeric_State': clean_str(row.get('Oligomeric State', '')),
                'Interface_Type': int_type,
                'Protein_chain': pro_chan,
                'DNA_chain': dna_chan,
                'RNA_chain': rna_chan,
                'Interface_Formed': determine_interface_formed(target_dir, pdb_id, {'Interface_Type': int_type}),

                'BSA_complex': m_c['bsa'],
                'FNP_complex': m_c['fnp'],
                'FBU_complex': m_c['fbu'],
                'LD_complex': m_c['ld'],

                'BSA_subunit1': m_p['bsa'],
                'FNP_subunit1': m_p['fnp'],
                'FBU_subunit1': m_p['fbu'],
                'LD_subunit1': m_p['ld'],

                'BSA_subunit2': m_part['bsa'],
                'FNP_subunit2': m_part['fnp'],
                'FBU_subunit2': m_part['fbu'],
                'LD_subunit2': m_part['ld'],
                'Output_Directory': str(target_dir)
            }

            records_unique_na.append(entry)
            records_batch.append(entry)

    # 2. PP (pro_pro_interface Sheet)
    if 'pro_pro_interface' in xls.sheet_names:
        df_pp_meta = pd.read_excel(xls, sheet_name='pro_pro_interface', engine='openpyxl')
        df_pp_meta = df_pp_meta[~df_pp_meta['PDB ID'].astype(str).str.startswith('#')].dropna(subset=['PDB ID'])

        for idx, row in df_pp_meta.iterrows():
            pdb_id = clean_str(row.get('PDB ID')).upper()
            c1_raw = clean_str(row.get('Chain 1'))
            c2_raw = clean_str(row.get('Chain 2'))
            c1 = c1_raw.split(":")[0].strip() if ":" in c1_raw else c1_raw
            c2 = c2_raw.split(":")[1].strip() if ":" in c2_raw else c2_raw

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

            method = clean_exp_method(row.get('Experimental Methods'), default_method="X-RAY DIFFRACTION")

            entry = {
                'PDB_ID': pdb_id,
                'Source_Dataset': '13.03.2026',
                'Exp_Method': method,
                'Protein_Name': clean_str(row.get('Protein_Name', row.get('Protein', ''))),
                'Source_Organism': clean_str(row.get('Source Organism', '')),
                'Protein_Length': clean_str(row.get('Protein length', '')),
                'Resolution': row.get('Resolution (Å)', ''),
                'Oligomeric_State': clean_str(row.get('Oligomeric State', '')),
                'Interface_Type': 'Protein-Protein',
                'Chain_1': c1,
                'Chain_2': c2,
                'Interface_Formed': determine_interface_formed(target_dir, pdb_id, {'Interface_Type': 'Protein-Protein', 'Chain_1': c1, 'Chain_2': c2}),

                'BSA_complex': m_c['bsa'],
                'FNP_complex': m_c['fnp'],
                'FBU_complex': m_c['fbu'],
                'LD_complex': m_c['ld'],

                'BSA_subunit1': m_1['bsa'],
                'FNP_subunit1': m_1['fnp'],
                'FBU_subunit1': m_1['fbu'],
                'LD_subunit1': m_1['ld'],

                'BSA_subunit2': m_2['bsa'],
                'FNP_subunit2': m_2['fnp'],
                'FBU_subunit2': m_2['fbu'],
                'LD_subunit2': m_2['ld'],
                'Output_Directory': str(target_dir)
            }

            records_unique_pp.append(entry)

    return pd.DataFrame(records_batch), pd.DataFrame(records_unique_na), pd.DataFrame(records_unique_pp)


# ==============================================================================
# Revision Runs Processing Suite (swarnava_TF_work/Interface/Revision)
# ==============================================================================

def process_revision_runs(rev_base_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parses Revision runs (Cryo-EM structures) for batch, unique NA, and unique PP datasets,
    enriching each record with full metadata from TFNRDv1.0_EM_cleaned_merged_metadata.xlsx.
    Returns: (df_batch, df_unique_na, df_unique_pp)
    """
    rev_base_dir = Path(rev_base_dir)

    batch_dir = rev_base_dir / "prince_results"
    na_dir = rev_base_dir / "prince_results_unique_na"
    pp_dir = rev_base_dir / "prince_results_unique_pp"

    meta_excel = rev_base_dir / "TFNRDv1.0_EM_cleaned_merged_metadata.xlsx"
    meta_lookup = {}
    if meta_excel.exists():
        try:
            df_m = pd.read_excel(meta_excel, sheet_name=0)
            for pid, group in df_m.groupby('Entry ID'):
                pid = clean_str(pid).upper()
                if not pid:
                    continue

                title = ''
                for v in group['Structure Title'].dropna():
                    if clean_str(v): title = clean_str(v); break
                if not title:
                    for v in group['Macromolecule Name'].dropna():
                        if clean_str(v): title = clean_str(v); break

                organism = ''
                for v in group['Source Organism'].dropna():
                    if clean_str(v): organism = clean_str(v); break

                oligo = ''
                for v in group['Oligomeric State'].dropna():
                    if clean_str(v): oligo = clean_str(v); break

                res = ''
                for v in group['EM Resolution (Å)'].dropna():
                    if clean_str(v): res = v; break

                length = ''
                for v in group['Total Number of Polymer Residues per Deposited Model'].dropna():
                    if clean_str(v): length = v; break

                meta_lookup[pid] = {
                    'Protein_Name': title,
                    'Source_Organism': organism,
                    'Oligomeric_State': oligo,
                    'Resolution': res,
                    'Protein_Length': length,
                    'Exp_Method': 'ELECTRON MICROSCOPY'
                }
        except Exception as e:
            logger.warning(f"Could not read EM metadata file {meta_excel}: {e}")

    def process_dir(results_dir, sum_csv_name, default_int_type):
        rows = []
        if not results_dir.exists():
            return pd.DataFrame(rows)

        # Load summary CSV if present for chain metadata
        sum_csv_path = results_dir / sum_csv_name
        sum_meta = {}
        if sum_csv_path.exists():
            try:
                df_s = pd.read_csv(sum_csv_path)
                for _, r in df_s.iterrows():
                    out_name = Path(str(r.get('Output_Directory', ''))).name
                    if out_name:
                        sum_meta[out_name] = dict(r)
            except Exception as e:
                logger.warning(f"Could not parse summary CSV {sum_csv_path}: {e}")

        subdirs = [d for d in results_dir.iterdir() if d.is_dir()]
        for target_dir in sorted(subdirs, key=lambda x: x.name):
            folder_name = target_dir.name
            pdb_id = folder_name.split('_')[0].upper()
            csv_meta = sum_meta.get(folder_name, {})

            meta_entry = meta_lookup.get(pdb_id, {})
            exp_method = meta_entry.get('Exp_Method', 'ELECTRON MICROSCOPY')
            res_val = meta_entry.get('Resolution', '')
            prot_len = meta_entry.get('Protein_Length', '')
            prot_name = meta_entry.get('Protein_Name', '')
            organism = meta_entry.get('Source_Organism', '')
            oligo_state = clean_str(csv_meta.get('Oligomeric_State', meta_entry.get('Oligomeric_State', '')))

            p_chain = clean_str(csv_meta.get('Protein_chain'))
            d_chain = clean_str(csv_meta.get('DNA_chain'))
            r_chain = clean_str(csv_meta.get('RNA_chain'))

            c1 = clean_str(csv_meta.get('Protein_1', csv_meta.get('Subunit1_Target', csv_meta.get('Subunit1_Chains', ''))))
            c2 = clean_str(csv_meta.get('Protein_2', csv_meta.get('Subunit2_Candidates', csv_meta.get('Subunit2_Chains', ''))))

            if not c1 and not c2 and '_p' in folder_name and '_q' in folder_name:
                try:
                    c1 = folder_name.split('_p')[1].split('_q')[0]
                    c2 = folder_name.split('_q')[1]
                except Exception:
                    pass

            int_type = clean_str(csv_meta.get('Interface_Type', default_int_type))
            if not int_type:
                int_type = default_int_type

            complex_f = target_dir / f"{pdb_id}.int"
            protein_f = target_dir / f"{pdb_id}P.int"
            sub2_f = target_dir / f"{pdb_id}R.int"
            if not sub2_f.exists():
                sub2_f = target_dir / f"{pdb_id}D.int"

            if not complex_f.exists():
                alt_cmplx = list(target_dir.glob(f"{pdb_id}_*.int"))
                if alt_cmplx:
                    complex_f = alt_cmplx[0]

            m_c = calculate_interface_metrics(complex_f)
            m_1 = calculate_interface_metrics(protein_f)
            m_2 = calculate_interface_metrics(sub2_f) if sub2_f.exists() else {'bsa': 0.0, 'fnp': 0.0, 'fbu': 0.0, 'ld': 0.0}

            int_formed = determine_interface_formed(target_dir, pdb_id, {'Interface_Type': int_type, 'Chain_1': c1, 'Chain_2': c2})

            entry = {
                'PDB_ID': pdb_id,
                'Source_Dataset': 'Revision_CryoEM',
                'Exp_Method': exp_method,
                'Protein_Name': prot_name,
                'Source_Organism': organism,
                'Protein_Length': prot_len,
                'Resolution': res_val,
                'Oligomeric_State': oligo_state,
                'Interface_Type': int_type,
                'Protein_chain': p_chain,
                'DNA_chain': d_chain,
                'RNA_chain': r_chain,
                'Chain_1': c1,
                'Chain_2': c2,
                'Folder_Name': folder_name,
                'Interface_Formed': int_formed,

                'BSA_complex': m_c['bsa'],
                'FNP_complex': m_c['fnp'],
                'FBU_complex': m_c['fbu'],
                'LD_complex': m_c['ld'],

                'BSA_subunit1': m_1['bsa'],
                'FNP_subunit1': m_1['fnp'],
                'FBU_subunit1': m_1['fbu'],
                'LD_subunit1': m_1['ld'],

                'BSA_subunit2': m_2['bsa'],
                'FNP_subunit2': m_2['fnp'],
                'FBU_subunit2': m_2['fbu'],
                'LD_subunit2': m_2['ld'],
                'Output_Directory': str(target_dir)
            }
            rows.append(entry)

        return pd.DataFrame(rows)

    df_batch = process_dir(batch_dir, 'prince_batch_summary.csv', 'Protein-NA/PP')
    df_na = process_dir(na_dir, 'prince_unique_na_batch_summary.csv', 'Protein-NA')
    df_pp = process_dir(pp_dir, 'prince_unique_pp_batch_summary.csv', 'Protein-Protein')

    return df_batch, df_na, df_pp


# ==============================================================================
# Summary Statistics & Machine Learning Dataset Generator
# ==============================================================================

def generate_summary_statistics(dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Generates summary statistics (Count, Mean, Std, Min, Max) for BSA, FNP, FBU, LD metrics."""
    summary_rows = []
    for cat_name, df in dfs.items():
        if df.empty:
            continue
        metric_cols = [c for c in df.columns if c.startswith(('BSA_', 'FNP_', 'FBU_', 'LD_'))]
        for col in metric_cols:
            series = df[col].dropna()
            if len(series) == 0:
                continue
            summary_rows.append({
                "Dataset_Sheet": cat_name,
                "Metric": col,
                "Count": len(series),
                "Mean": round(series.mean(), 2),
                "Std": round(series.std(), 2),
                "Min": round(series.min(), 2),
                "Max": round(series.max(), 2),
            })
    return pd.DataFrame(summary_rows)


def build_machine_learning_dataset(df_batch: pd.DataFrame, df_na: pd.DataFrame, df_pp: pd.DataFrame) -> pd.DataFrame:
    """
    Constructs a clean, flat, machine-learning ready DataFrame combining all structural interface features.
    """
    all_dfs = []
    for name, df in [('Batch', df_batch), ('Unique_NA', df_na), ('Unique_PP', df_pp)]:
        if not df.empty:
            df_copy = df.copy()
            df_copy['Dataset_Source_Category'] = name
            all_dfs.append(df_copy)

    if not all_dfs:
        return pd.DataFrame()

    df_ml = pd.concat(all_dfs, ignore_index=True)

    # Standardize float formats
    float_cols = [c for c in df_ml.columns if c.startswith(('BSA_', 'FNP_', 'FBU_', 'LD_'))]
    for col in float_cols:
        df_ml[col] = pd.to_numeric(df_ml[col], errors='coerce').fillna(0.0).round(2)

    return df_ml


# ==============================================================================
# Main Pipeline Orchestration
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Unified TF-NRD Structural Interface Analysis Suite & ML Dataset Generator."
    )
    parser.add_argument(
        "--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help=f"Target output directory for Supplementary files (default: {DEFAULT_OUTPUT_DIR})."
    )
    parser.add_argument(
        "--old_base_dir", type=Path, default=OLD_BASE_DIR,
        help=f"Directory for older 13.03.2026 runs (default: {OLD_BASE_DIR})."
    )
    parser.add_argument(
        "--revision_base_dir", type=Path, default=REVISION_BASE_DIR,
        help=f"Directory for Revision Cryo-EM runs (default: {REVISION_BASE_DIR})."
    )

    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("UNIFIED TF-NRD STRUCTURAL INTERFACE ANALYSIS SUITE")
    logger.info("=" * 70)
    logger.info(f"Older Runs Dir   : {args.old_base_dir}")
    logger.info(f"Revision Runs Dir: {args.revision_base_dir}")
    logger.info(f"Output Directory : {out_dir}")
    logger.info("=" * 70)

    # 1. Process Old Runs (13.03.2026 - X-ray & NMR)
    logger.info("--- Processing Older 13.03.2026 Runs (X-ray & NMR) ---")
    old_batch, old_na, old_pp = process_old_runs(args.old_base_dir)
    logger.info(f"Old Runs Processed: Batch={len(old_batch)}, Unique_NA={len(old_na)}, Unique_PP={len(old_pp)}")

    # 2. Process Revision Runs (Revision - Cryo-EM)
    logger.info("--- Processing Revision Runs (Cryo-EM) ---")
    rev_batch, rev_na, rev_pp = process_revision_runs(args.revision_base_dir)
    logger.info(f"Revision Runs Processed: Batch={len(rev_batch)}, Unique_NA={len(rev_na)}, Unique_PP={len(rev_pp)}")

    # 3. Combine Datasets
    combined_batch = pd.concat([old_batch, rev_batch], ignore_index=True) if not (old_batch.empty and rev_batch.empty) else pd.DataFrame()
    combined_na = pd.concat([old_na, rev_na], ignore_index=True) if not (old_na.empty and rev_na.empty) else pd.DataFrame()
    combined_pp = pd.concat([old_pp, rev_pp], ignore_index=True) if not (old_pp.empty and rev_pp.empty) else pd.DataFrame()

    dfs_dict = {
        'Table S1.A': combined_batch,
        'Table S1.B': combined_na,
        'Table S1.C': combined_pp
    }

    df_stats = generate_summary_statistics(dfs_dict)
    dfs_dict['Summary_Statistics'] = df_stats

    # 4. Generate Machine Learning Ready Dataset
    logger.info("--- Generating Machine-Learning Ready Datasets ---")
    df_ml = build_machine_learning_dataset(combined_batch, combined_na, combined_pp)

    ml_csv_path = out_dir / "merged_interface_features_ml.csv"
    ml_json_path = out_dir / "merged_interface_features_ml.json"

    if not df_ml.empty:
        df_ml.to_csv(ml_csv_path, index=False)
        logger.info(f"[SUCCESS] Exported ML dataset CSV ({len(df_ml)} entries) to: {ml_csv_path}")

        df_ml.to_json(ml_json_path, orient='records', indent=2)
        logger.info(f"[SUCCESS] Exported ML dataset JSON to: {ml_json_path}")

    # 5. Export Structured JSON Datasets
    logger.info("--- Exporting JSON Datasets ---")
    consolidated_json = {
        sname: df_sheet.to_dict(orient="records")
        for sname, df_sheet in dfs_dict.items()
        if not df_sheet.empty
    }
    all_json_path = out_dir / "interface_features_all_tables.json"
    with open(all_json_path, "w", encoding="utf-8") as jf:
        json.dump(consolidated_json, jf, indent=2, ensure_ascii=False)
    logger.info(f"[SUCCESS] Exported all tables JSON to: {all_json_path}")

    # Export individual table JSONs
    table_json_mapping = {
        'Table S1.A': out_dir / "table_s1_a_batch_interfaces.json",
        'Table S1.B': out_dir / "table_s1_b_unique_na_interfaces.json",
        'Table S1.C': out_dir / "table_s1_c_unique_pp_interfaces.json",
        'Summary_Statistics': out_dir / "interface_summary_statistics.json"
    }
    for sname, json_path in table_json_mapping.items():
        if sname in dfs_dict and not dfs_dict[sname].empty:
            with open(json_path, "w", encoding="utf-8") as jf:
                json.dump(dfs_dict[sname].to_dict(orient="records"), jf, indent=2, ensure_ascii=False)
            logger.info(f"[SUCCESS] Exported {sname} JSON -> {json_path}")

    # 6. Export Supplementary Tables Excel Workbook (Supplementary_Tables.xlsx only)
    supp_tables_xlsx_path = out_dir / "Supplementary_Tables.xlsx"
    logger.info("--- Exporting Supplementary Tables Excel Workbook ---")
    try:
        with pd.ExcelWriter(supp_tables_xlsx_path, engine='openpyxl') as writer:
            for sname, df_sheet in dfs_dict.items():
                if not df_sheet.empty:
                    df_sheet.to_excel(writer, sheet_name=sname, index=False)
                    logger.info(f"Wrote {len(df_sheet)} rows to sheet '{sname}' in {supp_tables_xlsx_path.name}")
        logger.info(f"[SUCCESS] Created Supplementary Tables Excel workbook at: {supp_tables_xlsx_path}")
    except Exception as e:
        logger.warning(f"Could not write Excel file {supp_tables_xlsx_path}: {e}")

    print("\n" + "=" * 70)
    print("UNIFIED INTERFACE ANALYSIS PIPELINE COMPLETE")
    print("=" * 70)
    print(f"  Table S1.A (Batch Run)    : {len(combined_batch)} entries")
    print(f"  Table S1.B (Unique NA)    : {len(combined_na)} entries")
    print(f"  Table S1.C (Unique PP)    : {len(combined_pp)} entries")
    print(f"  ML Dataset (Combined)     : {len(df_ml)} entries")
    print(f"  Supplementary Target Dir  : {out_dir.resolve()}")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()