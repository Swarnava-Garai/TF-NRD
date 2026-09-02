#!/usr/bin/env python3
"""
Interface_calculations.py
=========================
Universal grouped interface calculator using NACCESS.
Supports chain groups from text lists (e.g. New_list.txt, List_3.txt) or tabular files.
Generates individual and combined complex interface files (.int, .asa, .pdb) and JSON metrics.
Automatically converts CIF/mmCIF files to PDB on the fly for NACCESS.

Usage:
  # Run on test directory:
  python3 Interface_calculations.py -b test -i test/List_3.txt -o test/output_naccess -j test/output_naccess/bsa_summary.json
"""

import sys
import shutil
import logging
import argparse
import json
import time
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

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

UTILS_DIR = PROJECT_ROOT / "utils"
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

try:
    from mmcif_clean_reader import parse_file
except ImportError:
    try:
        from utils.mmcif_clean_reader import parse_file
    except ImportError:
        parse_file = None

DEFAULT_BASE_DIR = PROJECT_ROOT / "test"
DEFAULT_INPUT_FILE = DEFAULT_BASE_DIR / "List_3.txt"
DEFAULT_OUTPUT_DIR = DEFAULT_BASE_DIR / "output_naccess"
DEFAULT_JSON_OUTPUT = DEFAULT_OUTPUT_DIR / "bsa_interface_summary.json"
CUTOFF = 0.1


# ---------------- PDB & CIF FORMATTERS / PARSERS ---------------- #

def row_to_pdb_line(row: dict, atom_idx: int) -> str:
    """Formats an atom_site data row from mmcif_clean_reader into standard PDB ATOM record."""
    rec = row.get('_atom_site.group_PDB', 'ATOM')
    rec = f"{rec:<6}"[:6]

    atom_name = row.get('_atom_site.auth_atom_id') or row.get('_atom_site.label_atom_id', '')
    elem = row.get('_atom_site.type_symbol', '').strip()

    if len(atom_name) < 4 and len(elem) == 1:
        atom_field = f" {atom_name:<3}"
    else:
        atom_field = f"{atom_name:<4}"
    atom_field = atom_field[:4]

    resname = row.get('_atom_site.auth_comp_id') or row.get('_atom_site.label_comp_id', '')
    resname = f"{resname:>3}"[:3]

    chain = row.get('_atom_site.auth_asym_id') or row.get('_atom_site.label_asym_id', 'A')
    chain_char = chain[0] if chain else 'A'

    resseq_raw = row.get('_atom_site.auth_seq_id') or row.get('_atom_site.label_seq_id', '1')
    try:
        resseq = int(resseq_raw)
        resseq_str = f"{resseq:4d}"
    except (ValueError, TypeError):
        resseq_str = f"{resseq_raw:>4}"[:4]

    ins = row.get('_atom_site.pdbx_PDB_ins_code', '.')
    ins_char = ins if ins not in ('.', '?') else ' '

    try:
        x = float(row.get('_atom_site.Cartn_x', 0.0))
        y = float(row.get('_atom_site.Cartn_y', 0.0))
        z = float(row.get('_atom_site.Cartn_z', 0.0))
    except (ValueError, TypeError):
        x = y = z = 0.0

    try:
        occ = float(row.get('_atom_site.occupancy', 1.0))
    except (ValueError, TypeError):
        occ = 1.0

    try:
        bfac = float(row.get('_atom_site.B_iso_or_equiv', 0.0))
    except (ValueError, TypeError):
        bfac = 0.0

    elem_str = f"{elem:>2}"[:2]
    return f"{rec}{atom_idx % 100000:5d} {atom_field} {resname} {chain_char}{resseq_str}{ins_char}   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}          {elem_str}  \n"


def parse_pdbline(atom_line: str) -> List[str]:
    arr = [''] * 6
    arr[0] = atom_line[:4]
    arr[1] = atom_line[4:11].strip()
    arr[2] = atom_line[11:17].strip()
    arr[3] = atom_line[17:20].strip()
    arr[4] = atom_line[21] if len(atom_line) > 21 else ''  # chain ID
    arr[5] = atom_line[22:26].strip()
    return arr


def find_structure_file(pdb_id: str, base_dir: Path) -> Optional[Path]:
    """Finds CIF or PDB structure file for a given pdb_id."""
    candidates = [
        base_dir / f"{pdb_id}.cif",
        base_dir / f"{pdb_id.upper()}.cif",
        base_dir / f"{pdb_id.lower()}.cif",
        base_dir / f"{pdb_id}.mmcif",
        base_dir / f"{pdb_id.upper()}.mmcif",
        base_dir / f"{pdb_id.lower()}.mmcif",
        base_dir / f"{pdb_id}.pdb",
        base_dir / f"{pdb_id.upper()}.pdb",
        base_dir / f"{pdb_id.lower()}.pdb",
        PROJECT_ROOT / "test" / f"{pdb_id.upper()}.cif",
        PROJECT_ROOT / "input_data" / "RCSB_PDB" / f"{pdb_id.upper()}.cif",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def fetch_atomline(structure_path: Path, chain_group: str) -> str:
    """
    Extracts ATOM lines for the given chain_group from either a CIF or PDB file.
    """
    ext = structure_path.suffix.lower()
    if ext in ('.cif', '.mmcif'):
        if parse_file is None:
            raise ImportError("mmcif_clean_reader module could not be imported.")
        cols, rows = parse_file(str(structure_path))
        lines = []
        idx = 1
        for r in rows:
            chain = r.get('_atom_site.auth_asym_id') or r.get('_atom_site.label_asym_id', '')
            if chain and (chain in chain_group or any(c in chain_group for c in chain)):
                lines.append(row_to_pdb_line(r, idx))
                idx += 1
        if not lines:
            logger.error(f"Chain(s) '{chain_group}' not found in {structure_path.name}")
            return ""
        return "".join(lines)
    else:
        # Standard PDB
        lines = []
        found = False
        with open(structure_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    parsed = parse_pdbline(line)
                    if parsed[4] in chain_group:
                        lines.append(line)
                        found = True
        if not found:
            logger.error(f"Chain(s) '{chain_group}' not found in {structure_path.name}")
            return ""
        return "".join(lines)


# ---------------- INTERFACE CORE ---------------- #

def format_values(val: Any) -> str:
    val_str = str(val)
    return (' ' * (6 - len(val_str))) + val_str if len(val_str) < 6 else val_str


def generate_interface_atomfile(complex_asa: Path, subunit_asa: Path, out_int: Path, output_dir: Path, cutoff: float = CUTOFF):
    with open(complex_asa, "r", encoding="utf-8", errors="replace") as f:
        complex_lines = f.readlines()
    with open(subunit_asa, "r", encoding="utf-8", errors="replace") as f:
        subunit_lines = f.readlines()

    # Build fast O(1) lookup map from complex lines
    complex_map = {}
    for c in complex_lines:
        if c.strip():
            key = c[:55]
            try:
                complex_map[key] = float(c[55:63].strip())
            except (ValueError, IndexError):
                pass

    output = ""
    total_area = 0.0

    for s in subunit_lines:
        if not s.strip():
            continue

        try:
            sub_area = float(s[55:63].strip())
        except (ValueError, IndexError):
            continue

        total_area += sub_area
        key = s[:55]

        if key in complex_map:
            comp_area = complex_map[key]
            diff = sub_area - comp_area

            if diff > cutoff:
                output += (
                    s[:55]
                    + format_values(f"{sub_area:.0f}")
                    + format_values(f"{comp_area:.0f}")
                    + format_values(f"{diff:.0f}")
                    + "\n"
                )

    with open(out_int, "w", encoding="utf-8") as f:
        f.write(output)

    asa_final = output_dir / "ASA_FINAL"
    with open(asa_final, "a", encoding="utf-8") as f:
        f.write(f"{subunit_asa.name}\t{total_area:.0f}\n")


def calc_interface_area(int_file: Path) -> float:
    """Calculates total interface buried area by summing the last column of .int file."""
    area = 0.0
    with open(int_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            parts = line_str.split()
            if len(parts) >= 3:
                try:
                    area += float(parts[-1])
                except ValueError:
                    pass
    return area


def move_result_files(pdb_id: str, result_folder: Path, output_dir: Path):
    result_folder.mkdir(parents=True, exist_ok=True)
    for f in output_dir.iterdir():
        if f.is_file() and f.name.startswith(pdb_id + "_") and not f.name.startswith("BSA_FINAL") and not f.name.startswith("ASA_FINAL"):
            shutil.move(str(f), str(result_folder / f.name))


# ---------------- MAIN INTERFACE FUNCTION ---------------- #

def run_interface(pdb_id: str, group1: str, group2: str, base_dir: Path, output_dir: Path, cutoff: float = CUTOFF) -> Optional[Dict[str, Any]]:
    struct_file = find_structure_file(pdb_id, base_dir)
    if struct_file is None:
        logger.error(f"Structure file for '{pdb_id}' not found in {base_dir}")
        return None

    logger.info(f"Processing {pdb_id} ({struct_file.name}) -> {group1} vs {group2}")

    g1_lines = fetch_atomline(struct_file, group1)
    g2_lines = fetch_atomline(struct_file, group2)

    if not g1_lines or not g2_lines:
        logger.error(f"Failed to fetch atom lines for {pdb_id} ({group1} or {group2})")
        return None

    g1_pdb = output_dir / f"{pdb_id}_{group1}.pdb"
    g2_pdb = output_dir / f"{pdb_id}_{group2}.pdb"
    complex_pdb = output_dir / f"{pdb_id}_{group1}{group2}.pdb"

    g1_pdb.write_text(g1_lines, encoding="utf-8")
    g2_pdb.write_text(g2_lines, encoding="utf-8")
    complex_pdb.write_text(g1_lines + g2_lines, encoding="utf-8")

    # Run naccess directly inside output_dir
    for pdb_file in [g1_pdb, g2_pdb, complex_pdb]:
        subprocess.run(["naccess", pdb_file.name], cwd=output_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

    g1_asa = output_dir / f"{pdb_id}_{group1}.asa"
    g2_asa = output_dir / f"{pdb_id}_{group2}.asa"
    complex_asa = output_dir / f"{pdb_id}_{group1}{group2}.asa"

    if not all(x.exists() for x in [g1_asa, g2_asa, complex_asa]):
        logger.error(f"Missing NACCESS output for {pdb_id} ({group1} vs {group2})")
        return None

    g1_int = output_dir / f"{pdb_id}_{group1}.int"
    g2_int = output_dir / f"{pdb_id}_{group2}.int"
    complex_int = output_dir / f"{pdb_id}_{group1}{group2}.int"

    # ---- generate individual interface files ----
    generate_interface_atomfile(complex_asa, g1_asa, g1_int, output_dir, cutoff=cutoff)
    generate_interface_atomfile(complex_asa, g2_asa, g2_int, output_dir, cutoff=cutoff)

    # ---- generate combined complex interface file ----
    with open(complex_int, "w", encoding="utf-8") as out:
        if g1_int.exists():
            with open(g1_int, "r", encoding="utf-8") as f:
                out.write(f.read())
        if g2_int.exists():
            with open(g2_int, "r", encoding="utf-8") as f:
                out.write(f.read())

    # ---- area calculation ----
    area1 = calc_interface_area(g1_int) if g1_int.exists() else 0.0
    area2 = calc_interface_area(g2_int) if g2_int.exists() else 0.0
    total = area1 + area2

    bsa_final = output_dir / "BSA_FINAL"
    with open(bsa_final, "a", encoding="utf-8") as f:
        f.write(f"{pdb_id}\t{group1}\t{group2}\t{total:.0f}\t{area1:.0f}\t{area2:.0f}\n")

    result_folder = output_dir / f"{pdb_id}_{group1}_{group2}"
    move_result_files(pdb_id, result_folder, output_dir)

    return {
        "PDB_ID": pdb_id,
        "Group_1": group1,
        "Group_2": group2,
        "Total_BSA": round(total, 2),
        "BSA_Subunit1": round(area1, 2),
        "BSA_Subunit2": round(area2, 2),
        "Result_Folder": str(result_folder)
    }


def load_pairs_from_input(input_file: Path, base_dir: Path) -> List[Tuple[str, str, str]]:
    """Loads (pdb_id, group1, group2) pairs from text or tabular input."""
    pairs = []
    ext = input_file.suffix.lower()

    if ext in ('.csv', '.xlsx', '.xls', '.json'):
        import pandas as pd
        if ext == '.csv':
            df = pd.read_csv(input_file)
        elif ext == '.json':
            df = pd.read_json(input_file)
        else:
            df = pd.read_excel(input_file)

        for _, row in df.iterrows():
            pdb_id = str(row.get('PDB_ID', row.get('PDB ID', ''))).strip().upper()
            p_chain = str(row.get('Protein_chain', '')).replace('nan', '').strip()
            d_chain = str(row.get('DNA_chain', '')).replace('nan', '').strip()
            r_chain = str(row.get('RNA_chain', '')).replace('nan', '').strip()

            if p_chain and d_chain:
                pairs.append((pdb_id, p_chain, d_chain))
            elif p_chain and r_chain:
                pairs.append((pdb_id, p_chain, r_chain))
    else:
        # Text file mode (e.g. List_3.txt or New_list.txt)
        # Check if corresponding CSV exists in base_dir for single-column PDB IDs
        csv_lookup = {}
        for candidate_csv in [base_dir / "test_nr_EM_TF_with_chain.csv", PROJECT_ROOT / "test" / "test_nr_EM_TF_with_chain.csv"]:
            if candidate_csv.exists():
                import pandas as pd
                try:
                    df_c = pd.read_csv(candidate_csv)
                    for _, r in df_c.iterrows():
                        pid = str(r.get('PDB_ID', '')).strip().upper()
                        p = str(r.get('Protein_chain', '')).replace('nan', '').strip()
                        d = str(r.get('DNA_chain', '')).replace('nan', '').strip()
                        r_ch = str(r.get('RNA_chain', '')).replace('nan', '').strip()
                        csv_lookup[pid] = (p, d or r_ch)
                except Exception:
                    pass

        with open(input_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line_str = line.strip()
                if not line_str or line_str.startswith("#"):
                    continue

                parts = line_str.split()
                pdb_id = parts[0].upper()

                if len(parts) >= 2:
                    chain_def = parts[1].replace('"', ':')
                    chain_groups = [g.strip() for g in chain_def.split(':') if g.strip()]
                    if len(chain_groups) >= 2:
                        group1 = chain_groups[0]
                        for target_partner in chain_groups[1:]:
                            pairs.append((pdb_id, group1, target_partner))
                elif pdb_id in csv_lookup:
                    g1, g2 = csv_lookup[pdb_id]
                    if g1 and g2:
                        pairs.append((pdb_id, g1, g2))

    return pairs


def parse_args():
    parser = argparse.ArgumentParser(
        description="Universal grouped interface calculator using NACCESS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-b", "--base_dir",
        dest="base_dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Base directory containing structure CIF/PDB files."
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_file",
        type=Path,
        default=DEFAULT_INPUT_FILE,
        help="Path to input text/CSV/Excel file containing PDB IDs and chain groupings."
    )
    parser.add_argument(
        "-o", "--output_dir",
        dest="output_dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Target output directory for calculation results."
    )
    parser.add_argument(
        "-j", "--json",
        dest="json_output",
        type=Path,
        default=DEFAULT_JSON_OUTPUT,
        help="Path to export structured JSON metrics summary."
    )
    parser.add_argument(
        "-c", "--cutoff",
        dest="cutoff",
        type=float,
        default=CUTOFF,
        help="Interface area cutoff threshold in Å²."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    input_file = Path(args.input_file).resolve()
    output_dir = Path(args.output_dir).resolve()
    json_output = Path(args.json_output).resolve() if args.json_output else None

    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("NACCESS GROUPED INTERFACE CALCULATOR")
    logger.info("=" * 60)
    logger.info(f"Structure Dir : {base_dir}")
    logger.info(f"Input List    : {input_file}")
    logger.info(f"Output Dir    : {output_dir}")

    for fname in ["BSA_FINAL", "ASA_FINAL"]:
        fpath = output_dir / fname
        if fpath.exists():
            fpath.unlink()

    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)

    pairs = load_pairs_from_input(input_file, base_dir)
    logger.info(f"Identified {len(pairs)} chain group pairs to calculate.")

    results = []
    for pdb_id, group1, group2 in pairs:
        res = run_interface(pdb_id, group1, group2, base_dir, output_dir, cutoff=args.cutoff)
        if res:
            results.append(res)

    if json_output:
        json_output.parent.mkdir(parents=True, exist_ok=True)
        with open(json_output, "w", encoding="utf-8") as jf:
            json.dump(results, jf, indent=2)
        logger.info(f"Saved structured JSON summary ({len(results)} pairs) -> {json_output}")

    logger.info("### ALL CALCULATIONS COMPLETED ###")


if __name__ == "__main__":
    main()
