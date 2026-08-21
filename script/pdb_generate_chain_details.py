#!/home/labuser/anaconda3/bin/python
"""
pdb_generate_chain_details.py
------------------------------
Generates detailed chain mappings for biomolecular structures (X-ray, NMR, Cryo-EM).
Uses pathlib for robust PATH handling and cKDTree for spatial contact queries (default <= 5.0 Å).

Extracts:
  1. Experimental Method (X-ray, NMR, Cryo-EM)
  2. Raw chains present in the CIF structure (Raw_Protein_chain, Raw_DNA_chain, Raw_RNA_chain)
  3. Interacting chains determined via cKDTree spatial contact queries (Interacting_Protein_chain, Interacting_DNA_chain, Interacting_RNA_chain)
  4. Consolidated target chains (Protein_chain, DNA_chain, RNA_chain)
"""

import sys
import json
import csv
import argparse
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Add utils directory to sys.path
SCRIPT_DIR = Path(__file__).resolve().parent
UTILS_DIR = SCRIPT_DIR.parent / 'utils'
if UTILS_DIR.exists() and str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

try:
    from mmcif_clean_reader import parse_file
except ImportError:
    print(f"[ERROR] Could not import mmcif_clean_reader from {UTILS_DIR}", file=sys.stderr)
    sys.exit(1)

# Standard Residue Sets & Backbone Identifier Atoms
AA_LIST = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEC', 'PYL'
}
RNA_LIST = {'A', 'U', 'C', 'G', 'RA', 'RU', 'RC', 'RG', 'I'}
DNA_LIST = {'DA', 'DT', 'DC', 'DG', 'DI', 'DU', 'T'}

PROTEIN_BACKBONE_ATOMS = {"CA", "N", "C", "O"}
NUCLEIC_BACKBONE_ATOMS = {"P", "O5'", "C5'", "C4'", "C3'", "O3'", "C1'"}


def extract_exp_method(cif_path: Path) -> str:
    """Extracts experimental method (X-RAY, NMR, ELECTRON MICROSCOPY, etc.) from mmCIF file header."""
    if not cif_path.exists():
        return 'UNKNOWN'
    try:
        with cif_path.open('r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('_exptl.method'):
                    parts = line.strip().split(None, 1)
                    if len(parts) > 1:
                        return parts[1].strip("'\" ")
    except Exception as e:
        logging.warning(f"Failed to parse experimental method for {cif_path.name}: {e}")
    return 'UNKNOWN'


def analyze_cif_chains_kdtree(cif_path: Path, cutoff=5.0) -> dict:
    """
    Parses an mmCIF file (X-ray, NMR, or Cryo-EM) and uses cKDTree spatial queries
    to determine Raw and Interacting Protein, DNA, and RNA chains.
    """
    cif_path = Path(cif_path)
    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    cols, rows = parse_file(str(cif_path))
    chain_res = {}
    chain_atoms = {}
    chain_coords = {}

    chain_key = '_atom_site.auth_asym_id'
    res_key = '_atom_site.auth_comp_id'
    atom_key = '_atom_site.auth_atom_id'
    x_key, y_key, z_key = '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'
    model_key = '_atom_site.pdbx_PDB_model_num'

    # Filter to model 1 for NMR / multi-model structures
    first_model = None

    for row in rows:
        model_num = str(row.get(model_key, '1')).strip()
        if first_model is None:
            first_model = model_num
        elif model_num != first_model:
            continue

        ch = (row.get(chain_key) or row.get('_atom_site.label_asym_id') or '').strip()
        res = (row.get(res_key) or row.get('_atom_site.label_comp_id') or '').strip()
        atom = (row.get(atom_key) or row.get('_atom_site.label_atom_id') or '').strip()

        if not ch or not res:
            continue

        chain_res.setdefault(ch, set()).add(res)
        chain_atoms.setdefault(ch, set()).add(atom)

        try:
            xyz = (float(row[x_key]), float(row[y_key]), float(row[z_key]))
            chain_coords.setdefault(ch, []).append(xyz)
        except (ValueError, KeyError, TypeError):
            pass

    raw_protein = set()
    raw_dna = set()
    raw_rna = set()

    for ch, resnames in chain_res.items():
        atoms = chain_atoms.get(ch, set())
        n_aa = sum(1 for r in resnames if r in AA_LIST)
        n_dna = sum(1 for r in resnames if r in DNA_LIST)
        n_rna = sum(1 for r in resnames if r in RNA_LIST)

        if n_aa / max(1, len(resnames)) > 0.3 or any(r in AA_LIST for r in resnames):
            raw_protein.add(ch)
        elif n_dna > 0:
            raw_dna.add(ch)
        elif n_rna > 0:
            raw_rna.add(ch)
        else:
            if bool(atoms & NUCLEIC_BACKBONE_ATOMS) and not bool(atoms & PROTEIN_BACKBONE_ATOMS):
                raw_rna.add(ch)

    trees = {ch: cKDTree(np.array(coords)) for ch, coords in chain_coords.items() if len(coords) > 0}

    int_protein = set()
    int_dna = set()
    int_rna = set()

    # Query Protein vs DNA spatial contacts
    for p in raw_protein:
        if p not in trees:
            continue
        for d in raw_dna:
            if d not in trees:
                continue
            hits = trees[p].query_ball_tree(trees[d], r=cutoff)
            if any(hits):
                int_protein.add(p)
                int_dna.add(d)

        # Query Protein vs RNA spatial contacts
        for r in raw_rna:
            if r not in trees:
                continue
            hits = trees[p].query_ball_tree(trees[r], r=cutoff)
            if any(hits):
                int_protein.add(p)
                int_rna.add(r)

    def fmt_chains(chain_set):
        return "".join(sorted(chain_set)) if chain_set else ""

    p_interacting_str = fmt_chains(int_protein)
    d_interacting_str = fmt_chains(int_dna)
    r_interacting_str = fmt_chains(int_rna)

    p_raw_str = fmt_chains(raw_protein)
    d_raw_str = fmt_chains(raw_dna)
    r_raw_str = fmt_chains(raw_rna)

    # Consolidated target chains (Interacting if available, fallback to Raw)
    p_final_str = p_interacting_str if p_interacting_str else p_raw_str
    d_final_str = d_interacting_str if d_interacting_str else d_raw_str
    r_final_str = r_interacting_str if r_interacting_str else r_raw_str

    exp_method = extract_exp_method(cif_path)

    return {
        'PDB_ID': cif_path.stem.upper(),
        'Exp_Method': exp_method,
        'Raw_Protein_chain': p_raw_str,
        'Raw_DNA_chain': d_raw_str,
        'Raw_RNA_chain': r_raw_str,
        'Interacting_Protein_chain': p_interacting_str,
        'Interacting_DNA_chain': d_interacting_str,
        'Interacting_RNA_chain': r_interacting_str,
        'Protein_chain': p_final_str,
        'DNA_chain': d_final_str,
        'RNA_chain': r_final_str
    }


def main():
    default_test_dir = SCRIPT_DIR.parent / 'test'
    default_cif_dir = default_test_dir if default_test_dir.exists() else SCRIPT_DIR.parent / 'cif'

    parser = argparse.ArgumentParser(
        description="Generates raw and interacting chain details for X-ray, NMR, and Cryo-EM mmCIF structures using pathlib & cKDTree."
    )
    parser.add_argument(
        "-i", "--cif_dir", type=Path, default=default_cif_dir,
        help=f"Directory containing CIF structure files (default: {default_cif_dir})."
    )
    parser.add_argument(
        "-l", "--list_file", type=Path, default=None,
        help="Path to text file containing target PDB IDs."
    )
    parser.add_argument(
        "-o", "--output_dir", type=Path, default=None,
        help="Output directory to write results (default: same as cif_dir)."
    )
    parser.add_argument(
        "-c", "--cutoff", type=float, default=5.0,
        help="Spatial contact distance cutoff in Angstroms (default: 5.0 Å)."
    )

    args = parser.parse_args()

    cif_dir = Path(args.cif_dir)
    if not cif_dir.exists():
        print(f"[ERROR] CIF directory not found: {cif_dir}", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.output_dir) if args.output_dir else cif_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Determine CIF file list
    cif_files = []
    if args.list_file:
        list_path = Path(args.list_file)
        if list_path.exists():
            lines = list_path.read_text().splitlines()
            pdb_ids = [line.strip().upper() for line in lines if line.strip()]
            cif_files = [cif_dir / f"{pid}.cif" for pid in pdb_ids]
    
    if not cif_files:
        # Check for standard List files inside cif_dir
        for list_name in ['List_2.txt', 'List_3.txt', 'List_55_EM_TF.txt', 'List.txt']:
            check_list = cif_dir / list_name
            if check_list.exists():
                lines = check_list.read_text().splitlines()
                pdb_ids = [line.strip().upper() for line in lines if line.strip()]
                cif_files = [cif_dir / f"{pid}.cif" for pid in pdb_ids if (cif_dir / f"{pid}.cif").exists()]
                if cif_files:
                    break

    if not cif_files:
        cif_files = sorted(list(cif_dir.glob('*.cif')))

    if not cif_files:
        print(f"[ERROR] No CIF files found in {cif_dir}", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print("PDB CHAIN DETAILS GENERATOR (X-ray, NMR, Cryo-EM)")
    print("=" * 70)
    print(f"CIF Directory: {cif_dir}")
    print(f"Output Dir   : {out_dir}")
    print(f"Contact Cutoff: {args.cutoff} Å")
    print(f"Found {len(cif_files)} CIF files to process.")
    print("=" * 70 + "\n")

    results = []

    for cif_path in cif_files:
        if not cif_path.exists():
            print(f"[SKIP] File not found: {cif_path}")
            continue

        try:
            info = analyze_cif_chains_kdtree(cif_path, cutoff=args.cutoff)
            results.append(info)
            print(
                f"[{info['PDB_ID']}] Method: {info['Exp_Method']:20s} | "
                f"Raw (P:{info['Raw_Protein_chain']:10s} D:{info['Raw_DNA_chain']:5s} R:{info['Raw_RNA_chain']:5s}) | "
                f"Interacting (P:{info['Interacting_Protein_chain']:10s} D:{info['Interacting_DNA_chain']:5s} R:{info['Interacting_RNA_chain']:5s})"
            )
        except Exception as e:
            print(f"[ERROR] Failed to process {cif_path.name}: {e}", file=sys.stderr)

    if not results:
        print("[WARNING] No valid structures processed.")
        return

    # Write output files (CSV, JSON, Excel)
    output_csv = out_dir / 'test_nr_EM_TF_with_chain.csv' if 'test' in str(cif_dir) else out_dir / 'nr_EM_TF_with_chain.csv'
    output_json = out_dir / 'test_nr_EM_TF_with_chain.json' if 'test' in str(cif_dir) else out_dir / 'nr_EM_TF_with_chain.json'
    output_excel = out_dir / 'test_nr_EM_TF_with_chain.xlsx' if 'test' in str(cif_dir) else out_dir / 'nr_EM_TF_with_chain.xlsx'

    fieldnames = [
        'PDB_ID', 'Exp_Method',
        'Raw_Protein_chain', 'Raw_DNA_chain', 'Raw_RNA_chain',
        'Interacting_Protein_chain', 'Interacting_DNA_chain', 'Interacting_RNA_chain',
        'Protein_chain', 'DNA_chain', 'RNA_chain'
    ]

    # Write CSV
    with output_csv.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"\n[SUCCESS] Wrote CSV summary ({len(results)} entries) to: {output_csv}")

    # Write JSON
    with output_json.open('w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    print(f"[SUCCESS] Wrote JSON summary to: {output_json}")

    # Write / Update Excel sheet 'Interface_NA'
    try:
        df_res = pd.DataFrame(results)
        if output_excel.exists():
            with pd.ExcelWriter(output_excel, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
                df_res.to_excel(writer, sheet_name='Interface_NA', index=False)
        else:
            with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
                df_res.to_excel(writer, sheet_name='Interface_NA', index=False)
        print(f"[SUCCESS] Updated Excel workbook (Sheet: 'Interface_NA') at: {output_excel}")
    except Exception as e:
        print(f"[WARNING] Could not update Excel file {output_excel}: {e}")

    print("=" * 70 + "\n")


if __name__ == '__main__':
    main()
