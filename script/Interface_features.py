#!/home/labuser/anaconda3/bin/python
"""
Interface_features.py
---------------------
Unified interface execution and feature collection framework for PRince (Biomolecular Complex Interface Calculator).

Combines:
  1. Standard batch interface calculation (Sheet: 'Interface_NA' in nr_EM_TF_with_chain.xlsx)
  2. Unique Protein-Nucleic Acid interface calculation (Sheet: 'Unique_Interface_NA')
  3. Unique Protein-Protein interface calculation (Sheet: 'Unique_Interface_PP')
  4. Automated calculation and collection of interface features (BSA, FNP, FBU, LD) into a multi-sheet Excel workbook.

Usage:
  # Default: Run standard batch mode + feature collection
  python Interface_features.py

  # Run Unique Protein-NA mode
  python Interface_features.py --mode unique_na

  # Run Unique Protein-Protein mode with custom thread count
  python Interface_features.py --mode unique_pp --num_workers 12

  # Run all interface modes (batch, unique_na, unique_pp) sequentially followed by feature collection
  python Interface_features.py --mode all

  # Shortcut flags
  python Interface_features.py --unique_na
  python Interface_features.py --unique_pp
  python Interface_features.py --all
"""

import os
import sys
import csv
import argparse
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Default global paths
PRINCE_BIN_DEFAULT = Path('/home/labuser/Projects/PRince/bin/prince')
BASE_DIR_DEFAULT = Path('/home/labuser/Projects/PhD_projects/swarnava_TF_work/Interface/Revision')
CIF_DIR_DEFAULT = BASE_DIR_DEFAULT / 'cif'
EXCEL_PATH_DEFAULT = BASE_DIR_DEFAULT / 'nr_EM_TF_with_chain.xlsx'
OUTPUT_XLSX_DEFAULT = BASE_DIR_DEFAULT / 'prince_interface_features.xlsx'


# ==============================================================================
# Helper Functions & Parsers
# ==============================================================================

def extract_chain_chars(chain_str):
    """
    Extracts unique 1-character chain IDs from chain string (preserving PDB formatting).
    Filters out delimiters like colons (:), commas (,), semicolons (;), and whitespace.
    Handles numeric chain IDs (e.g. integer or float 0 -> '0').
    """
    if pd.isna(chain_str) or chain_str is None:
        return ""

    if isinstance(chain_str, float):
        if chain_str.is_integer():
            chain_str = int(chain_str)
        else:
            return ""

    s = str(chain_str).strip()
    if not s or s.lower() in ('nan', 'none', 'null', 'na'):
        return ""
    if s.endswith('.0') and s[:-2].isdigit():
        s = s[:-2]

    seen = set()
    chars = []
    for char in s:
        if char not in seen and char not in (':', ',', ';', ' ', '\t'):
            seen.add(char)
            chars.append(char)
    return "".join(chars)


def run_prince_protein_nucleic_acid(prince_bin, cif_path, prot_chains, rna_chains, dna_chains, output_dir):
    """
    Executes PRince CLI runner for Protein-Nucleic Acid interface calculation.
    """
    cmd = [str(prince_bin), "-i", str(cif_path), "-o", str(output_dir)]
    if prot_chains:
        cmd.extend(["-p", prot_chains])
    if rna_chains:
        cmd.extend(["-r", rna_chains])
    if dna_chains:
        cmd.extend(["-d", dna_chains])

    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        return True, res.stdout
    except subprocess.CalledProcessError as e:
        return False, e.stderr


def run_prince_protein_protein(prince_bin, cif_path, prot_chains_1, prot_chains_2, output_dir):
    """
    Executes PRince CLI runner for Protein-Protein interface calculation between Subunit 1 and Subunit 2.
    """
    cmd = [str(prince_bin), "-i", str(cif_path), "-p", prot_chains_1, "-q", prot_chains_2, "-o", str(output_dir)]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        return True, res.stdout
    except subprocess.CalledProcessError as e:
        return False, e.stderr


def parse_prince_res_file(res_path):
    """
    Parses summary metrics from generated PRince .res report file.
    """
    summary = {}
    if not res_path.exists():
        return summary

    with open(res_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if 'Total Interface Area' in line or 'Buried Surface Area' in line or 'Total Buried' in line or 'Interface area' in line:
                summary['BSA_line'] = line
    return summary


def get_interacting_chains_and_bsa(filepath):
    """
    Parses a PRince .int file to extract:
      1. Total Buried Surface Area (BSA in Å²)
      2. List of unique chain IDs that have buried surface area > 0.0 Å²
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return 0.0, []

    total_bsa = 0.0
    chains = set()

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('ATOM'):
                parts = line.split()
                if len(parts) >= 10:
                    try:
                        chain = parts[4]
                        bsa = float(parts[-1])
                        total_bsa += bsa
                        if bsa > 0.0:
                            chains.add(chain)
                    except (IndexError, ValueError):
                        pass

    return round(total_bsa, 2), sorted(chains)


# ==============================================================================
# Execution Modes: Batch, Unique NA, Unique PP
# ==============================================================================

def run_batch_task(task_data):
    """Worker task for standard batch PRince interface calculations."""
    prince_bin, row_num, total_rows, pdb_id, run_type, p_chars, d_chars, r_chars, p1_chars, p2_chars, cif_path, out_dir = task_data

    if run_type == 'NA':
        na_type = "Protein-DNA-RNA" if (d_chars and r_chars) else ("Protein-DNA" if d_chars else "Protein-RNA")
        print(f"[{row_num}/{total_rows}] [{pdb_id}] Running Batch {na_type} interface (p='{p_chars}', d='{d_chars}', r='{r_chars}')...")
        success, msg = run_prince_protein_nucleic_acid(prince_bin, cif_path, p_chars, r_chars, d_chars, out_dir)
        status = "SUCCESS" if success else "FAILED"
        res_file = out_dir / f"{pdb_id}.res"
        res_summary = parse_prince_res_file(res_file) if success else {}
        return {
            'PDB_ID': pdb_id,
            'Interface_Type': na_type,
            'Protein_chain': p_chars,
            'DNA_chain': d_chars,
            'RNA_chain': r_chars,
            'Subunit1_Chains': p_chars,
            'Subunit2_Chains': f"DNA:{d_chars} RNA:{r_chars}".strip(),
            'Status': status,
            'Output_Directory': str(out_dir),
            'BSA_Summary': res_summary.get('BSA_line', '')
        }
    else:  # PP
        print(f"[{row_num}/{total_rows}] [{pdb_id}] Running Batch Protein-Protein interface (Subunit1='{p1_chars}', Subunit2='{p2_chars}')...")
        success, msg = run_prince_protein_protein(prince_bin, cif_path, p1_chars, p2_chars, out_dir)
        status = "SUCCESS" if success else "FAILED"
        res_file = out_dir / f"{pdb_id}.res"
        res_summary = parse_prince_res_file(res_file) if success else {}
        return {
            'PDB_ID': pdb_id,
            'Interface_Type': 'Protein-Protein',
            'Protein_chain': p_chars,
            'DNA_chain': d_chars,
            'RNA_chain': r_chars,
            'Subunit1_Chains': p1_chars,
            'Subunit2_Chains': p2_chars,
            'Status': status,
            'Output_Directory': str(out_dir),
            'BSA_Summary': res_summary.get('BSA_line', '')
        }


def execute_batch_mode(prince_bin, excel_path, cif_dir, output_base_dir, num_workers):
    """Executes standard batch PRince calculations based on sheet 'Interface_NA'."""
    sheet_name = 'Interface_NA'
    out_dir_base = output_base_dir / 'prince_results'
    out_dir_base.mkdir(parents=True, exist_ok=True)

    print(f"\n" + "=" * 70)
    print(f"MODE: Standard Batch Interface Calculation (Sheet: '{sheet_name}')")
    print("=" * 70)

    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except Exception as e:
        print(f"[ERROR] Failed to read sheet '{sheet_name}' from {excel_path}: {e}", file=sys.stderr)
        return

    total_rows = len(df)
    tasks = []

    for idx, row in df.iterrows():
        row_num = idx + 1
        pdb_id = str(row.get('PDB_ID', row.get('PDB ID', ''))).strip().upper()
        if not pdb_id or pdb_id.lower() in ('nan', 'none', 'null'):
            continue

        cif_path = cif_dir / f"{pdb_id}.cif"
        if not cif_path.exists():
            print(f"[{row_num}/{total_rows}] [{pdb_id}] [SKIP] CIF file missing: {cif_path}")
            continue

        raw_p = row.get('Protein_chain')
        raw_d = row.get('DNA_chain')
        raw_r = row.get('RNA_chain')
        raw_p1 = row.get('Protein_chain_1') or row.get('Protein_sub1') or row.get('Subunit1') or row.get('p1')
        raw_p2 = row.get('Protein_chain_2') or row.get('Protein_sub2') or row.get('Subunit2') or row.get('p2') or row.get('q')

        p_chars = extract_chain_chars(raw_p)
        d_chars = extract_chain_chars(raw_d)
        r_chars = extract_chain_chars(raw_r)
        p1_chars = extract_chain_chars(raw_p1)
        p2_chars = extract_chain_chars(raw_p2)

        # 1. NA interface task
        if d_chars or r_chars:
            na_out_dir = out_dir_base / f"{pdb_id}_NA"
            tasks.append((prince_bin, row_num, total_rows, pdb_id, 'NA', p_chars, d_chars, r_chars, '', '', cif_path, na_out_dir))

        # 2. PP interface task
        if p1_chars and p2_chars:
            pp_out_dir = out_dir_base / f"{pdb_id}_PP"
            tasks.append((prince_bin, row_num, total_rows, pdb_id, 'PP', p_chars, d_chars, r_chars, p1_chars, p2_chars, cif_path, pp_out_dir))

    summary_records = []
    print(f"Executing {len(tasks)} batch PRince runs using {num_workers} parallel worker threads...")

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(run_batch_task, task) for task in tasks]
        for future in as_completed(futures):
            try:
                res = future.result()
                summary_records.append(res)
            except Exception as exc:
                print(f"[ERROR] Worker generated an exception: {exc}", file=sys.stderr)

    summary_csv = out_dir_base / 'prince_batch_summary.csv'
    fieldnames = ['PDB_ID', 'Interface_Type', 'Protein_chain', 'DNA_chain', 'RNA_chain', 'Subunit1_Chains', 'Subunit2_Chains', 'Status', 'Output_Directory', 'BSA_Summary']
    with open(summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_records)

    print(f"[COMPLETE] Standard Batch execution finished. Summary saved to {summary_csv}")


def run_unique_na_task(task_data):
    """Worker task for Unique Protein-NA interface calculations."""
    prince_bin, row_num, total_rows, pdb_id, raw_p, p_chars, d_chars, r_chars, tf_status, oligo_state, cif_path, out_dir = task_data
    na_type = "Protein-DNA-RNA" if (d_chars and r_chars) else ("Protein-DNA" if d_chars else ("Protein-RNA" if r_chars else "Protein-NA"))
    print(f"[{row_num}/{total_rows}] [{pdb_id}] Running Unique NA interface (Prot='{p_chars}', DNA='{d_chars}', RNA='{r_chars}')...")

    success, msg = run_prince_protein_nucleic_acid(prince_bin, cif_path, p_chars, r_chars, d_chars, out_dir)
    status = "SUCCESS" if success else "FAILED"
    res_file = out_dir / f"{pdb_id}.res"
    res_summary = parse_prince_res_file(res_file) if success else {}

    return {
        'PDB_ID': pdb_id,
        'Protein_chain': p_chars,
        'DNA_chain': d_chars,
        'RNA_chain': r_chars,
        'TF_Status': tf_status,
        'Oligomeric_State': oligo_state,
        'Interface_Type': na_type,
        'Status': status,
        'Output_Directory': str(out_dir),
        'BSA_Summary': res_summary.get('BSA_line', '')
    }


def execute_unique_na_mode(prince_bin, excel_path, cif_dir, output_base_dir, num_workers):
    """Executes Unique Protein-NA calculations based on sheet 'Unique_Interface_NA'."""
    sheet_name = 'Unique_Interface_NA'
    out_dir_base = output_base_dir / 'prince_results_unique_na'
    out_dir_base.mkdir(parents=True, exist_ok=True)

    print(f"\n" + "=" * 70)
    print(f"MODE: Unique Protein-Nucleic Acid Interface Calculation (Sheet: '{sheet_name}')")
    print("=" * 70)

    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except Exception as e:
        print(f"[ERROR] Failed to read sheet '{sheet_name}' from {excel_path}: {e}", file=sys.stderr)
        return

    total_rows = len(df)
    tasks = []
    skipped_count = 0

    for idx, row in df.iterrows():
        row_num = idx + 1
        pdb_id = str(row.get('PDB ID', '')).strip().upper()
        if not pdb_id or pdb_id.lower() in ('nan', 'none', 'null'):
            continue

        raw_p = row.get('Protein_chain')
        raw_d = row.get('DNA_chain')
        raw_r = row.get('RNA_chain')
        tf_status = row.get('TF_Status', '')
        oligo_state = row.get('Oligomeric state', '')

        if pd.isna(raw_p) or not str(raw_p).strip() or str(raw_p).strip().lower() in ('nan', 'none', 'null', 'na'):
            print(f"[{row_num}/{total_rows}] [{pdb_id}] [SKIP] Protein_chain is NaN or empty: '{raw_p}'")
            skipped_count += 1
            continue

        p_chars = extract_chain_chars(raw_p)
        d_chars = extract_chain_chars(raw_d)
        r_chars = extract_chain_chars(raw_r)

        cif_path = cif_dir / f"{pdb_id}.cif"
        if not cif_path.exists():
            print(f"[{row_num}/{total_rows}] [{pdb_id}] [SKIP] CIF file missing: {cif_path}")
            skipped_count += 1
            continue

        folder_tag = f"{pdb_id}_p{p_chars}"
        if d_chars:
            folder_tag += f"_d{d_chars}"
        if r_chars:
            folder_tag += f"_r{r_chars}"
        out_dir = out_dir_base / folder_tag

        tasks.append((prince_bin, row_num, total_rows, pdb_id, raw_p, p_chars, d_chars, r_chars, tf_status, oligo_state, cif_path, out_dir))

    summary_records = []
    print(f"Executing {len(tasks)} Unique NA PRince runs using {num_workers} parallel worker threads...")

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(run_unique_na_task, task) for task in tasks]
        for future in as_completed(futures):
            try:
                res = future.result()
                summary_records.append(res)
            except Exception as exc:
                print(f"[ERROR] Worker generated an exception: {exc}", file=sys.stderr)

    summary_csv = out_dir_base / 'prince_unique_na_batch_summary.csv'
    fieldnames = ['PDB_ID', 'Protein_chain', 'DNA_chain', 'RNA_chain', 'TF_Status', 'Oligomeric_State', 'Interface_Type', 'Status', 'Output_Directory', 'BSA_Summary']
    with open(summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_records)

    print(f"[COMPLETE] Unique NA calculation finished. Processed: {len(summary_records)} | Skipped: {skipped_count}")
    print(f"Summary saved to {summary_csv}")


def run_unique_pp_task(task_data):
    """Worker task for Unique Protein-Protein interface calculations."""
    prince_bin, row_num, total_rows, pdb_id, raw_p, target_p1, candidate_p2, cif_path, out_dir_base = task_data
    folder_tag = f"{pdb_id}_p{target_p1}_q{candidate_p2}"
    out_dir = out_dir_base / folder_tag
    print(f"[{row_num}/{total_rows}] [{pdb_id}] Running PP interface (Subunit1='{target_p1}', Subunit2='{candidate_p2}')...")

    success, msg = run_prince_protein_protein(prince_bin, cif_path, target_p1, candidate_p2, out_dir)
    status = "SUCCESS" if success else "FAILED"

    res_file = out_dir / f"{pdb_id}.res"
    res_summary = parse_prince_res_file(res_file) if success else {}

    cmplx_f = out_dir / f"{pdb_id}.int"
    c1_f = out_dir / f"{pdb_id}P.int"
    c2_f = out_dir / f"{pdb_id}R.int"
    if not c2_f.exists():
        c2_f = out_dir / f"{pdb_id}D.int"

    bsa_cmplx, cmplx_chains = get_interacting_chains_and_bsa(cmplx_f)
    bsa_sub1, sub1_chains = get_interacting_chains_and_bsa(c1_f)
    bsa_sub2, sub2_chains = get_interacting_chains_and_bsa(c2_f)

    forming_partners = sorted(list(set(sub2_chains)))
    has_interface = (bsa_cmplx > 0.0 or bsa_sub1 > 0.0 or bsa_sub2 > 0.0) and len(forming_partners) > 0
    forming_partner_str = ",".join(forming_partners) if forming_partners else ("NONE" if not has_interface else candidate_p2)

    return {
        'PDB_ID': pdb_id,
        'Protein_chain': str(raw_p).strip() if pd.notna(raw_p) else '',
        'Protein_1': target_p1,
        'Protein_2': candidate_p2,
        'Subunit1_Target': target_p1,
        'Subunit2_Candidates': candidate_p2,
        'Has_Interface': has_interface,
        'Forming_Partner_Chains': forming_partner_str,
        'BSA_Complex': bsa_cmplx,
        'BSA_Subunit1': bsa_sub1,
        'BSA_Subunit2': bsa_sub2,
        'Status': status,
        'Output_Directory': str(out_dir),
        'BSA_Summary': res_summary.get('BSA_line', '')
    }


def execute_unique_pp_mode(prince_bin, excel_path, cif_dir, output_base_dir, num_workers):
    """Executes Unique Protein-Protein calculations based on sheet 'Unique_Interface_PP'."""
    sheet_name = 'Unique_Interface_PP'
    out_dir_base = output_base_dir / 'prince_results_unique_pp'
    out_dir_base.mkdir(parents=True, exist_ok=True)

    print(f"\n" + "=" * 70)
    print(f"MODE: Unique Protein-Protein Interface Calculation (Sheet: '{sheet_name}')")
    print("=" * 70)

    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except Exception as e:
        print(f"[ERROR] Failed to read sheet '{sheet_name}' from {excel_path}: {e}", file=sys.stderr)
        return

    total_rows = len(df)
    tasks = []
    skipped_rows = 0

    for idx, row in df.iterrows():
        row_num = idx + 1
        pdb_id = str(row.get('PDB_ID', row.get('PDB ID', ''))).strip().upper()
        if not pdb_id or pdb_id.lower() in ('nan', 'none', 'null'):
            continue

        raw_p = row.get('Protein_chain')
        raw_p1 = row.get('Protein_1')
        raw_p2 = row.get('Protein_2')

        if pd.isna(raw_p) or not str(raw_p).strip() or str(raw_p).strip().lower() in ('nan', 'none', 'null', 'na'):
            print(f"[{row_num}/{total_rows}] [{pdb_id}] [SKIP] Protein_chain is NaN or empty: '{raw_p}'")
            skipped_rows += 1
            continue

        p_chars = extract_chain_chars(raw_p)
        p1_chars = extract_chain_chars(raw_p1)
        p2_chars = extract_chain_chars(raw_p2)

        cif_path = cif_dir / f"{pdb_id}.cif"
        if not cif_path.exists():
            print(f"[{row_num}/{total_rows}] [{pdb_id}] [SKIP] CIF file missing: {cif_path}")
            skipped_rows += 1
            continue

        run_pairs = []
        if p1_chars and p2_chars:
            run_pairs.append((p1_chars, p2_chars))
        elif raw_p and (':' in str(raw_p) or ',' in str(raw_p)):
            sep = ':' if ':' in str(raw_p) else ','
            parts = str(raw_p).split(sep)
            if len(parts) >= 2:
                sp1 = extract_chain_chars(parts[0])
                sp2 = extract_chain_chars(parts[1])
                if sp1 and sp2:
                    run_pairs.append((sp1, sp2))

        if not run_pairs:
            if len(p_chars) < 2:
                print(f"[{row_num}/{total_rows}] [{pdb_id}] [SKIP] Insufficient protein chains for PP interface: '{raw_p}' (P1='{p1_chars}', P2='{p2_chars}')")
                skipped_rows += 1
                continue

            for c in p_chars:
                target_chain = c
                candidate_partners = "".join([ch for ch in p_chars if ch != c])
                run_pairs.append((target_chain, candidate_partners))

        for target_p1, candidate_p2 in run_pairs:
            tasks.append((prince_bin, row_num, total_rows, pdb_id, raw_p, target_p1, candidate_p2, cif_path, out_dir_base))

    summary_records = []
    print(f"Executing {len(tasks)} Unique PP PRince runs using {num_workers} parallel worker threads...")

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(run_unique_pp_task, task) for task in tasks]
        for future in as_completed(futures):
            try:
                res = future.result()
                summary_records.append(res)
            except Exception as exc:
                print(f"[ERROR] Worker generated an exception: {exc}", file=sys.stderr)

    summary_csv = out_dir_base / 'prince_unique_pp_batch_summary.csv'
    fieldnames = [
        'PDB_ID', 'Protein_chain', 'Protein_1', 'Protein_2',
        'Subunit1_Target', 'Subunit2_Candidates', 'Has_Interface',
        'Forming_Partner_Chains', 'BSA_Complex', 'BSA_Subunit1', 'BSA_Subunit2',
        'Status', 'Output_Directory', 'BSA_Summary'
    ]

    with open(summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_records)

    print(f"[COMPLETE] Unique PP calculation finished. Processed runs: {len(summary_records)} | Skipped rows: {skipped_rows}")
    print(f"Summary saved to {summary_csv}")


# ==============================================================================
# Interface Feature Collection Module (BSA, FNP, FBU, LD)
# ==============================================================================

def calculate_interface_metrics(filepath):
    """
    Calculates interface metrics from a single PRince .int file:
      - BSA: Buried Surface Area (Å²)
      - FNP: Fraction of Non-Polar BSA (%)
      - FBU: Fraction of Buried Residues (%)
      - LD: Local Density (average contacts within 12 Å)
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

    return {'bsa': round(total_bsa, 2), 'fnp': fnp, 'fbu': fbu, 'ld': ld}


def load_excel_na_lookup(excel_path):
    """Loads metadata fallback lookup from nr_EM_TF_with_chain.xlsx."""
    lookup = {}
    excel_path = Path(excel_path)
    if not excel_path.exists():
        return lookup
    try:
        df_excel = pd.read_excel(excel_path, sheet_name='Interface_NA')
        for _, r in df_excel.iterrows():
            pid = str(r.get('PDB_ID', r.get('PDB ID', ''))).strip().upper()
            if pid and pid.lower() not in ('nan', 'none', 'null'):
                lookup[pid] = {
                    'Protein_chain': r.get('Protein_chain'),
                    'DNA_chain': r.get('DNA_chain'),
                    'RNA_chain': r.get('RNA_chain'),
                }
    except Exception as e:
        logging.warning(f"Could not parse Excel lookup from {excel_path}: {e}")
    return lookup


def parse_int_atom_chains(filepath):
    """Parses ATOM lines with BSA > 0 from PRince .int file to extract active chain IDs."""
    p_chains, d_chains, r_chains, other_chains = set(), set(), set(), set()
    if not Path(filepath).exists():
        return p_chains, d_chains, r_chains, other_chains

    AA_LIST = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'MSE', 'SEC', 'PYL'
    }
    RNA_LIST = {'A', 'U', 'C', 'G', 'RA', 'RU', 'RC', 'RG', 'I'}
    DNA_LIST = {'DA', 'DT', 'DC', 'DG', 'DI', 'DU', 'T'}

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


def determine_interface_formed(target_dir, pdb_id, meta):
    """Determines active interface chain string e.g. 'IJLM:OP' or 'NONE'."""
    complex_f = target_dir / f"{pdb_id}.int"
    if not complex_f.exists():
        c1 = str(meta.get('Subunit1_Chains', meta.get('Subunit1_Target', meta.get('Protein_1', '')))).strip()
        c2 = str(meta.get('Subunit2_Chains', meta.get('Subunit2_Candidates', meta.get('Protein_2', '')))).strip()
        if c1 and c2:
            alt_cmplx = target_dir / f"{pdb_id}_{c1}{c2}.int"
            if alt_cmplx.exists():
                complex_f = alt_cmplx

    if not complex_f.exists():
        return "NONE"

    p_ch, d_ch, r_ch, oth_ch = parse_int_atom_chains(complex_f)
    int_type = str(meta.get('Interface_Type', '')).strip()
    c1_meta = str(meta.get('Subunit1_Chains', meta.get('Subunit1_Target', meta.get('Protein_1', '')))).strip()
    c2_meta = str(meta.get('Subunit2_Chains', meta.get('Subunit2_Candidates', meta.get('Protein_2', '')))).strip()

    if int_type == 'Protein-Protein' or (c1_meta and c2_meta and not meta.get('DNA_chain') and not meta.get('RNA_chain')):
        s1_active = [c for c in c1_meta if c in p_ch or c in oth_ch]
        s2_active = [c for c in c2_meta if c in p_ch or c in oth_ch]

        if not s1_active and not s2_active and len(p_ch) >= 2:
            sorted_p = sorted(p_ch)
            s1_active = [sorted_p[0]]
            s2_active = sorted_p[1:]

        s1_str = "".join(sorted(set(s1_active)))
        s2_str = "".join(sorted(set(s2_active)))
        if s1_str and s2_str:
            return f"{s1_str}:{s2_str}"
        return "NONE"

    p_str = "".join(sorted(p_ch))
    na_ch = sorted(d_ch | r_ch)
    na_str = "".join(na_ch)

    if p_str and na_str:
        return f"{p_str}:{na_str}"
    elif p_str:
        return f"{p_str}:NONE"
    return "NONE"


def compute_na_interaction_chains(target_dir, pdb_id, p_chains_str, d_chains_str, r_chains_str):
    """Determines which protein chains interact with DNA, RNA, or both (dual)."""
    def extract_clean_chains(chain_str):
        if pd.isna(chain_str) or not chain_str:
            return []
        s = str(chain_str).strip()
        if not s or s.lower() in ('nan', 'none', 'null', 'na'):
            return []
        seen = set()
        res = []
        for c in s:
            if c not in seen and c not in (':', ',', ';', ' ', '\t'):
                seen.add(c)
                res.append(c)
        return res

    def parse_int_coords_by_chain(filepath):
        chain_coords = {}
        if not Path(filepath).exists():
            return chain_coords
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('ATOM'):
                    parts = line.split()
                    if len(parts) >= 10:
                        raw_ch = parts[4]
                        base_ch = raw_ch[0]
                        try:
                            x = float(parts[-6])
                            y = float(parts[-5])
                            z = float(parts[-4])
                            bsa = float(parts[-1])
                            if bsa > 0.0:
                                if base_ch not in chain_coords:
                                    chain_coords[base_ch] = []
                                chain_coords[base_ch].append([x, y, z])
                        except (ValueError, IndexError):
                            pass
        return {k: np.array(v) for k, v in chain_coords.items()}

    p_chains = extract_clean_chains(p_chains_str)
    d_chains = extract_clean_chains(d_chains_str)
    r_chains = extract_clean_chains(r_chains_str)

    if not p_chains:
        return ("NONE", "NONE", "NONE")

    complex_f = target_dir / f"{pdb_id}.int"
    d_f = target_dir / f"{pdb_id}D.int"
    r_f = target_dir / f"{pdb_id}R.int"

    cmplx_coords = parse_int_coords_by_chain(complex_f)
    d_coords_dict = parse_int_coords_by_chain(d_f)
    r_coords_dict = parse_int_coords_by_chain(r_f)

    dna_list = []
    for ch in d_chains:
        if ch in d_coords_dict:
            dna_list.append(d_coords_dict[ch])
        elif ch in cmplx_coords:
            dna_list.append(cmplx_coords[ch])
    dna_arr = np.vstack(dna_list) if dna_list else None

    rna_list = []
    for ch in r_chains:
        if ch in r_coords_dict:
            rna_list.append(r_coords_dict[ch])
        elif ch in cmplx_coords:
            rna_list.append(cmplx_coords[ch])
    rna_arr = np.vstack(rna_list) if rna_list else None

    dna_tree = cKDTree(dna_arr) if dna_arr is not None and len(dna_arr) > 0 else None
    rna_tree = cKDTree(rna_arr) if rna_arr is not None and len(rna_arr) > 0 else None

    dna_prot, rna_prot, dual_prot = [], [], []

    for p_ch in p_chains:
        if p_ch not in cmplx_coords:
            continue
        p_arr = cmplx_coords[p_ch]
        p_tree = cKDTree(p_arr)

        has_dna = False
        if dna_tree is not None:
            hits = p_tree.query_ball_tree(dna_tree, r=5.0)
            if any(hits):
                has_dna = True

        has_rna = False
        if rna_tree is not None:
            hits = p_tree.query_ball_tree(rna_tree, r=5.0)
            if any(hits):
                has_rna = True

        if has_dna:
            dna_prot.append(p_ch)
        if has_rna:
            rna_prot.append(p_ch)
        if has_dna and has_rna:
            dual_prot.append(p_ch)

    return (
        ",".join(sorted(dna_prot)) if dna_prot else "NONE",
        ",".join(sorted(rna_prot)) if rna_prot else "NONE",
        ",".join(sorted(dual_prot)) if dual_prot else "NONE"
    )


def process_prince_results_directory(results_dir, summary_csv_filename, excel_path):
    """Parses and calculates features for all subdirectories in a PRince results folder."""
    results_dir = Path(results_dir)
    if not results_dir.exists():
        logging.warning(f"Results directory does not exist: {results_dir}")
        return pd.DataFrame()

    summary_csv = results_dir / summary_csv_filename
    subfolder_meta = {}

    if summary_csv.exists():
        try:
            df_sum = pd.read_csv(summary_csv)
            for _, row in df_sum.iterrows():
                out_dir = Path(str(row.get('Output_Directory', ''))).name
                if out_dir:
                    subfolder_meta[out_dir] = dict(row)
        except Exception as e:
            logging.warning(f"Could not parse summary CSV {summary_csv}: {e}")

    excel_lookup = load_excel_na_lookup(excel_path)
    subdirs = [d for d in results_dir.iterdir() if d.is_dir()]
    logging.info(f"Collecting features for {len(subdirs)} subdirectories in {results_dir.name}")

    rows = []
    for target_dir in sorted(subdirs, key=lambda x: x.name):
        folder_name = target_dir.name
        meta = subfolder_meta.get(folder_name, {})

        pdb_id = meta.get('PDB_ID')
        if not pdb_id:
            pdb_id = folder_name.split('_')[0].upper()

        int_type = meta.get('Interface_Type')
        if not int_type or pd.isna(int_type):
            if 'Protein_1' in meta and 'Protein_2' in meta:
                int_type = 'Protein-Protein'
            elif meta.get('DNA_chain') or meta.get('RNA_chain'):
                int_type = 'Protein-NA'
            else:
                int_type = 'Unknown'

        entry = {
            'PDB_ID': pdb_id,
            'Folder_Name': folder_name,
            'Interface_Type': int_type,
            'Interface_Formed': determine_interface_formed(target_dir, pdb_id, meta),
        }

        for k in ['Protein_chain', 'Protein_1', 'Protein_2', 'DNA_chain', 'RNA_chain', 
                  'TF_Status', 'Oligomeric_State', 'Subunit1_Chains', 'Subunit2_Chains',
                  'Subunit1_Target', 'Subunit2_Candidates', 'Has_Interface', 'Forming_Partner_Chains',
                  'DNA_Interacting_Protein_Chains', 'RNA_Interacting_Protein_Chains', 'Dual_Interacting_Protein_Chains',
                  'Status']:
            if k in meta and pd.notna(meta[k]):
                entry[k] = meta[k]

        lookup_entry = excel_lookup.get(pdb_id, {})
        for chain_k in ['Protein_chain', 'DNA_chain', 'RNA_chain']:
            if chain_k not in entry or pd.isna(entry[chain_k]):
                if chain_k in lookup_entry and pd.notna(lookup_entry[chain_k]):
                    entry[chain_k] = lookup_entry[chain_k]

        p_str = entry.get('Protein_chain', '')
        d_str = entry.get('DNA_chain', '')
        r_str = entry.get('RNA_chain', '')
        if d_str or r_str:
            dna_prot, rna_prot, dual_prot = compute_na_interaction_chains(target_dir, pdb_id, p_str, d_str, r_str)
            entry['DNA_Interacting_Protein_Chains'] = dna_prot
            entry['RNA_Interacting_Protein_Chains'] = rna_prot
            entry['Dual_Interacting_Protein_Chains'] = dual_prot

        complex_f = target_dir / f"{pdb_id}.int"
        protein_f = target_dir / f"{pdb_id}P.int"
        sub2_f = target_dir / f"{pdb_id}R.int"
        if not sub2_f.exists():
            sub2_f = target_dir / f"{pdb_id}D.int"

        if not complex_f.exists():
            c1 = str(meta.get('Subunit1_Chains', meta.get('Subunit1_Target', meta.get('Protein_1', '')))).strip()
            c2 = str(meta.get('Subunit2_Chains', meta.get('Subunit2_Candidates', meta.get('Protein_2', '')))).strip()
            if c1 and c2:
                alt_cmplx = target_dir / f"{pdb_id}_{c1}{c2}.int"
                if alt_cmplx.exists():
                    complex_f = alt_cmplx
                alt_c1 = target_dir / f"{pdb_id}_{c1}.int"
                if alt_c1.exists():
                    protein_f = alt_c1
                alt_c2 = target_dir / f"{pdb_id}_{c2}.int"
                if alt_c2.exists():
                    sub2_f = alt_c2

        m_complex = calculate_interface_metrics(complex_f)
        m_protein = calculate_interface_metrics(protein_f)
        m_sub2 = calculate_interface_metrics(sub2_f) if sub2_f.exists() else None

        entry['BSA_complex'] = m_complex['bsa']
        entry['FNP_complex'] = m_complex['fnp']
        entry['FBU_complex'] = m_complex['fbu']
        entry['LD_complex'] = m_complex['ld']

        entry['BSA_subunit1'] = m_protein['bsa']
        entry['FNP_subunit1'] = m_protein['fnp']
        entry['FBU_subunit1'] = m_protein['fbu']
        entry['LD_subunit1'] = m_protein['ld']

        if m_sub2 is not None:
            entry['BSA_subunit2'] = m_sub2['bsa']
            entry['FNP_subunit2'] = m_sub2['fnp']
            entry['FBU_subunit2'] = m_sub2['fbu']
            entry['LD_subunit2'] = m_sub2['ld']
        else:
            entry['BSA_subunit2'] = 0.0
            entry['FNP_subunit2'] = 0.0
            entry['FBU_subunit2'] = 0.0
            entry['LD_subunit2'] = 0.0

        entry['Output_Directory'] = str(target_dir)
        rows.append(entry)

    return pd.DataFrame(rows)


def generate_summary_statistics(dfs):
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


def collect_features(active_modes, output_base_dir, excel_path, output_xlsx):
    """Collects features across all executed and existing directories and outputs multi-sheet Excel file."""
    print(f"\n" + "=" * 70)
    print(f"COLLECTING CALCULATED INTERFACE FEATURES")
    print("=" * 70)

    dfs = {}
    out_base = Path(output_base_dir)
    output_xlsx = Path(output_xlsx)

    # 1. Read existing sheets from Excel workbook if it already exists
    existing_sheets = {}
    if output_xlsx.exists():
        try:
            excel_file = pd.ExcelFile(output_xlsx)
            for sname in excel_file.sheet_names:
                if sname != 'Summary_Statistics':
                    existing_sheets[sname] = pd.read_excel(excel_file, sheet_name=sname)
        except Exception as e:
            logging.warning(f"Could not read existing Excel workbook {output_xlsx}: {e}")

    # 2. Check all potential result directories on disk (whether from current or previous runs)
    batch_dir = out_base / 'prince_results'
    if batch_dir.exists() and any(d.is_dir() for d in batch_dir.iterdir()):
        dfs['Batch_Interfaces'] = process_prince_results_directory(batch_dir, 'prince_batch_summary.csv', excel_path)
    elif 'Batch_Interfaces' in existing_sheets:
        dfs['Batch_Interfaces'] = existing_sheets['Batch_Interfaces']

    na_dir = out_base / 'prince_results_unique_na'
    if na_dir.exists() and any(d.is_dir() for d in na_dir.iterdir()):
        dfs['Unique_NA_Interfaces'] = process_prince_results_directory(na_dir, 'prince_unique_na_batch_summary.csv', excel_path)
    elif 'Unique_NA_Interfaces' in existing_sheets:
        dfs['Unique_NA_Interfaces'] = existing_sheets['Unique_NA_Interfaces']

    pp_dir = out_base / 'prince_results_unique_pp'
    if pp_dir.exists() and any(d.is_dir() for d in pp_dir.iterdir()):
        dfs['Unique_PP_Interfaces'] = process_prince_results_directory(pp_dir, 'prince_unique_pp_batch_summary.csv', excel_path)
    elif 'Unique_PP_Interfaces' in existing_sheets:
        dfs['Unique_PP_Interfaces'] = existing_sheets['Unique_PP_Interfaces']

    # Retain any custom existing sheets not covered above
    for sname, df_sheet in existing_sheets.items():
        if sname not in dfs:
            dfs[sname] = df_sheet

    # Filter out empty DataFrames
    dfs = {k: v for k, v in dfs.items() if not v.empty}

    if not dfs:
        print("[WARNING] No interface calculation directories or existing sheets found for feature collection.")
        return

    df_stats = generate_summary_statistics(dfs)

    output_xlsx.parent.mkdir(parents=True, exist_ok=True)

    logging.info(f"Writing Excel workbook to: {output_xlsx}")
    with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
        for sname, df in dfs.items():
            if not df.empty:
                df.to_excel(writer, sheet_name=sname, index=False)
                logging.info(f"Wrote {len(df)} rows to sheet '{sname}'")

        if not df_stats.empty:
            df_stats.to_excel(writer, sheet_name='Summary_Statistics', index=False)
            logging.info(f"Wrote {len(df_stats)} rows to sheet 'Summary_Statistics'")

    print("\n" + "=" * 70)
    print("PRINCE INTERFACE FEATURE COLLECTION COMPLETE")
    print("=" * 70)
    for sname, df in dfs.items():
        print(f"  {sname}: {len(df)} entries")
    print(f"  Output Excel Workbook: {output_xlsx.resolve()}")
    print("=" * 70 + "\n")


# ==============================================================================
# Main CLI Entry Point
# ==============================================================================

def main():
    default_workers = min(8, os.cpu_count() or 4)

    parser = argparse.ArgumentParser(
        description="Unified PRince interface execution and feature collection framework."
    )
    parser.add_argument(
        "-m", "--mode", type=str, choices=['batch', 'unique_na', 'unique_pp', 'all'], default='batch',
        help="Execution mode: 'batch' (default, standard batch run), 'unique_na' (unique NA interfaces), 'unique_pp' (unique PP interfaces), or 'all' (runs all modes)."
    )
    parser.add_argument(
        "--unique_na", action="store_true",
        help="Shortcut to run unique Protein-NA interface mode."
    )
    parser.add_argument(
        "--unique_pp", action="store_true",
        help="Shortcut to run unique Protein-Protein interface mode."
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Shortcut to run all interface modes."
    )
    parser.add_argument(
        "-w", "--num_workers", type=int, default=default_workers,
        help=f"Number of parallel CPU worker threads for PRince execution (default: {default_workers})."
    )
    parser.add_argument(
        "--excel_path", type=Path, default=EXCEL_PATH_DEFAULT,
        help=f"Path to input Excel spreadsheet (default: {EXCEL_PATH_DEFAULT})."
    )
    parser.add_argument(
        "--cif_dir", type=Path, default=CIF_DIR_DEFAULT,
        help=f"Path to directory containing structure CIF files (default: {CIF_DIR_DEFAULT})."
    )
    parser.add_argument(
        "--output_dir", type=Path, default=BASE_DIR_DEFAULT,
        help=f"Base output directory for results (default: {BASE_DIR_DEFAULT})."
    )
    parser.add_argument(
        "--prince_bin", type=Path, default=PRINCE_BIN_DEFAULT,
        help=f"Path to PRince binary executable (default: {PRINCE_BIN_DEFAULT})."
    )
    parser.add_argument(
        "--output_xlsx", type=Path, default=OUTPUT_XLSX_DEFAULT,
        help=f"Target Excel path for collected features (default: {OUTPUT_XLSX_DEFAULT})."
    )
    parser.add_argument(
        "--skip_collect", action="store_true",
        help="Skip automatic feature collection after PRince runs complete."
    )

    args = parser.parse_args()

    # Priority to explicit shortcut flags
    selected_mode = args.mode
    if args.all:
        selected_mode = 'all'
    elif args.unique_na:
        selected_mode = 'unique_na'
    elif args.unique_pp:
        selected_mode = 'unique_pp'

    prince_bin = Path(args.prince_bin)
    if not prince_bin.exists():
        print(f"[ERROR] PRince binary executable not found at: {prince_bin}", file=sys.stderr)
        sys.exit(1)

    excel_path = Path(args.excel_path)
    if not excel_path.exists():
        print(f"[ERROR] Input Excel file not found at: {excel_path}", file=sys.stderr)
        sys.exit(1)

    cif_dir = Path(args.cif_dir)
    output_dir = Path(args.output_dir)

    print("=" * 70)
    print("PRINCE UNIFIED INTERFACE FRAMEWORK")
    print("=" * 70)
    print(f"Selected Mode: {selected_mode}")
    print(f"PRince Binary: {prince_bin}")
    print(f"Input Excel  : {excel_path}")
    print(f"CIF Directory: {cif_dir}")
    print(f"Worker Threads: {args.num_workers}")
    print("=" * 70)

    modes_to_run = []
    if selected_mode == 'all':
        modes_to_run = ['batch', 'unique_na', 'unique_pp']
    else:
        modes_to_run = [selected_mode]

    for mode in modes_to_run:
        if mode == 'batch':
            execute_batch_mode(prince_bin, excel_path, cif_dir, output_dir, args.num_workers)
        elif mode == 'unique_na':
            execute_unique_na_mode(prince_bin, excel_path, cif_dir, output_dir, args.num_workers)
        elif mode == 'unique_pp':
            execute_unique_pp_mode(prince_bin, excel_path, cif_dir, output_dir, args.num_workers)

    if not args.skip_collect:
        collect_features(modes_to_run if selected_mode != 'all' else ['all'], output_dir, excel_path, args.output_xlsx)


if __name__ == '__main__':
    main()
