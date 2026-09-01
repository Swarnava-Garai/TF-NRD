#!/usr/bin/python3
# Universal grouped interface calculator using NACCESS
# Supports chain groups exactly as given in New_list.txt / List_1.txt
# Generates individual AND combined complex interface files
# Safe file handling (no folder nesting, no file mixing)
# Automatically converts CIF/mmCIF files to PDB on the fly for NACCESS

import sys
import shutil
import logging
import time
import subprocess
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Setup project directories and import CIF reader utility
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
UTILS_DIR = REPO_ROOT / "utils"
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

try:
    from mmcif_clean_reader import parse_file
except ImportError:
    from utils.mmcif_clean_reader import parse_file

BASE_DIR = Path("/home/labuser/Projects/PhD_projects/swarnava_TF_work/Interface/Revision/cif")
OUTPUT_DIR = BASE_DIR / "check"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = BASE_DIR / "New_list.txt"
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


def parse_pdbline(atom_line):
    arr = [''] * 6
    arr[0] = atom_line[:4]
    arr[1] = atom_line[4:11].strip()
    arr[2] = atom_line[11:17].strip()
    arr[3] = atom_line[17:20].strip()
    arr[4] = atom_line[21]  # chain ID
    arr[5] = atom_line[22:26].strip()
    return arr


def find_structure_file(pdb_id: str, base_dir: Path) -> Path:
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
        Path(f"{pdb_id}.cif"),
        Path(f"{pdb_id}.pdb"),
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def fetch_atomline(structure_path: Path, chain_group: str) -> str:
    """
    Extracts ATOM lines for the given chain_group from either a CIF or PDB file.
    If CIF, uses mmcif_clean_reader to parse and convert to standard PDB format.
    """
    ext = structure_path.suffix.lower()
    if ext in ('.cif', '.mmcif'):
        cols, rows = parse_file(str(structure_path))
        lines = []
        idx = 1
        for r in rows:
            chain = r.get('_atom_site.auth_asym_id') or r.get('_atom_site.label_asym_id', '')
            if chain and (chain in chain_group or any(c in chain_group for c in chain)):
                lines.append(row_to_pdb_line(r, idx))
                idx += 1
        if not lines:
            logging.error(f"Chain(s) '{chain_group}' not found in {structure_path}")
            return ""
        return "".join(lines)
    else:
        # Standard PDB
        lines = []
        found = False
        with open(structure_path) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    parsed = parse_pdbline(line)
                    if parsed[4] in chain_group:
                        lines.append(line)
                        found = True
        if not found:
            logging.error(f"Chain(s) '{chain_group}' not found in {structure_path}")
            return ""
        return "".join(lines)


# ---------------- INTERFACE CORE ---------------- #

def format_values(val):
    val = str(val)
    return (' ' * (6 - len(val))) + val if len(val) < 6 else val


def generate_interface_atomfile(complex_asa: Path, subunit_asa: Path, out_int: Path):
    complex_lines = open(complex_asa).readlines()
    subunit_lines = open(subunit_asa).readlines()

    output = ""
    total_area = 0.0

    for s in subunit_lines:
        if not s.strip():
            continue

        sub_area = float(s[55:63].strip())
        total_area += sub_area

        for c in complex_lines:
            if not c.strip():
                continue

            if s[:55] == c[:55]:
                comp_area = float(c[55:63].strip())
                diff = sub_area - comp_area

                if diff > CUTOFF:
                    output += (
                        s[:55]
                        + format_values(f"{sub_area:.0f}")
                        + format_values(f"{comp_area:.0f}")
                        + format_values(f"{diff:.0f}")
                        + "\n"
                    )
                break

    with open(out_int, "w") as f:
        f.write(output)

    asa_final = OUTPUT_DIR / "ASA_FINAL"
    with open(asa_final, "a") as f:
        f.write(f"{subunit_asa.name}\t{total_area:.0f}\n")


def calc_interface_area(int_file: Path) -> float:
    """Calculates total interface buried area by summing the last column of .int file."""
    area = 0.0
    with open(int_file) as f:
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


# ---------------- SAFE FILE HANDLING ---------------- #

def move_result_files(pdb_id: str, result_folder: Path):
    result_folder.mkdir(parents=True, exist_ok=True)
    for f in OUTPUT_DIR.iterdir():
        if f.is_file() and f.name.startswith(pdb_id + "_") and not f.name.startswith("BSA_FINAL") and not f.name.startswith("ASA_FINAL"):
            shutil.move(str(f), str(result_folder / f.name))


# ---------------- MAIN INTERFACE FUNCTION ---------------- #

def run_interface(pdb_id: str, group1: str, group2: str):
    struct_file = find_structure_file(pdb_id, BASE_DIR)
    if struct_file is None:
        logging.error(f"Structure file (CIF/PDB) for '{pdb_id}' not found in {BASE_DIR}")
        return

    logging.info(f"Processing {pdb_id} ({struct_file.name}) → {group1} vs {group2}")

    g1_lines = fetch_atomline(struct_file, group1)
    g2_lines = fetch_atomline(struct_file, group2)

    if not g1_lines or not g2_lines:
        logging.error(f"Failed to fetch atom lines for {pdb_id} ({group1} or {group2})")
        return

    g1_pdb = OUTPUT_DIR / f"{pdb_id}_{group1}.pdb"
    g2_pdb = OUTPUT_DIR / f"{pdb_id}_{group2}.pdb"
    complex_pdb = OUTPUT_DIR / f"{pdb_id}_{group1}{group2}.pdb"

    g1_pdb.write_text(g1_lines)
    g2_pdb.write_text(g2_lines)
    complex_pdb.write_text(g1_lines + g2_lines)

    # Run naccess directly inside OUTPUT_DIR
    for pdb_file in [g1_pdb, g2_pdb, complex_pdb]:
        subprocess.run(["naccess", pdb_file.name], cwd=OUTPUT_DIR, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

    g1_asa = OUTPUT_DIR / f"{pdb_id}_{group1}.asa"
    g2_asa = OUTPUT_DIR / f"{pdb_id}_{group2}.asa"
    complex_asa = OUTPUT_DIR / f"{pdb_id}_{group1}{group2}.asa"

    if not all(x.exists() for x in [g1_asa, g2_asa, complex_asa]):
        logging.error(f"Missing NACCESS output for {pdb_id} ({group1} vs {group2})")
        return

    g1_int = OUTPUT_DIR / f"{pdb_id}_{group1}.int"
    g2_int = OUTPUT_DIR / f"{pdb_id}_{group2}.int"
    complex_int = OUTPUT_DIR / f"{pdb_id}_{group1}{group2}.int"

    # ---- generate individual interface files ----
    generate_interface_atomfile(complex_asa, g1_asa, g1_int)
    generate_interface_atomfile(complex_asa, g2_asa, g2_int)

    # ---- generate combined complex interface file ----
    with open(complex_int, "w") as out:
        if g1_int.exists():
            with open(g1_int) as f:
                out.write(f.read())
        if g2_int.exists():
            with open(g2_int) as f:
                out.write(f.read())

    # ---- area calculation ----
    area1 = calc_interface_area(g1_int) if g1_int.exists() else 0.0
    area2 = calc_interface_area(g2_int) if g2_int.exists() else 0.0
    total = area1 + area2

    bsa_final = OUTPUT_DIR / "BSA_FINAL"
    with open(bsa_final, "a") as f:
        f.write(f"{pdb_id}\t{group1}\t{group2}\t{total:.0f}\t{area1:.0f}\t{area2:.0f}\n")

    result_folder = OUTPUT_DIR / f"{pdb_id}_{group1}_{group2}"
    move_result_files(pdb_id, result_folder)


# ---------------- DRIVER ---------------- #

def main():
    start_time = time.strftime("%d-%m-%Y %H:%M:%S", time.localtime())
    logging.info(f"\n{start_time}\n### START ###")

    for fname in ["BSA_FINAL", "ASA_FINAL"]:
        fpath = OUTPUT_DIR / fname
        if fpath.exists():
            fpath.unlink()

    if not INPUT_FILE.exists():
        logging.error(f"Input file not found: {INPUT_FILE}")
        return

    with open(INPUT_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue
            pdb_id = parts[0]
            chain_def = parts[1].replace('"', ':')
            chain_groups = [g.strip() for g in chain_def.split(':') if g.strip()]
            if len(chain_groups) < 2:
                continue
            group1 = chain_groups[0]

            # For each target partner (e.g., Protein vs DNA, then Protein vs RNA)
            for target_partner in chain_groups[1:]:
                run_interface(pdb_id, group1, target_partner)

    logging.info("### ALL CALCULATIONS COMPLETED ###")


if __name__ == "__main__":
    main()
