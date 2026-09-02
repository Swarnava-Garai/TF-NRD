#!/usr/bin/env python3
"""
tf_interface_common_uniqueness.py
----------------------------------
Unified Common & Unique Interface Analysis Suite for Biomolecular Complex Interfaces.

Merges and extends:
  1. common_atom.py: Extracts common/unique atoms per interface pair and exports common atom list text files.
  2. final_compare_difference_similarity.py: Computes interface atom counts, line counts, percentage similarity, and Jaccard indices.

Supports:
  - Input list mode (e.g. generated_int_paths.txt)
  - Auto-discovery mode (finds and pairs .int interface files within PDB folders)
  - Exports summary TSV, CSV, Excel, and JSON.
"""

import os
import sys
import csv
import json
import argparse
import logging
from pathlib import Path
from typing import List, Tuple, Set, Dict, Optional, Any
import pandas as pd

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

# Default project paths
DEFAULT_BASE_DIR = PROJECT_ROOT / "results" / "Interface"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Supplementary"
DEFAULT_INPUT_LIST = DEFAULT_BASE_DIR / "generated_int_paths.txt"


# ==============================================================================
# Interface Atom Parsing & Comparison Functions
# ==============================================================================

def parse_atoms(file_lines: List[str], mode: str = "atom_resnum") -> Set[Any]:
    """
    Parses ATOM lines from PRince .int file into a set of atom identifier tuples.
    - 'atom_resnum' (default): (atom_name, res_num)
    - 'full': (chain_id, res_name, res_num, atom_name)
    """
    atoms = set()
    for line in file_lines:
        if line.startswith("ATOM"):
            if len(line) >= 26:
                atom_name = line[12:16].strip()
                res_num = line[22:26].strip()
                res_name = line[17:20].strip() if len(line) >= 20 else ""
                chain_id = line[21].strip() if len(line) >= 22 else ""

                if mode == "full":
                    atoms.add((chain_id, res_name, res_num, atom_name))
                else:
                    atoms.add((atom_name, res_num))
    return atoms


def count_atom_lines(file_lines: List[str]) -> int:
    """Counts non-empty ATOM lines in file content."""
    return len([l for l in file_lines if l.startswith("ATOM") or l.strip()])


def compare_int_files(file1_lines: List[str], file2_lines: List[str], mode: str = "atom_resnum") -> Tuple[Dict[str, Any], Set[Any]]:
    """
    Compares two .int interface file contents and returns set metrics + common atoms set.
    """
    atoms1 = parse_atoms(file1_lines, mode=mode)
    atoms2 = parse_atoms(file2_lines, mode=mode)

    common_atoms = atoms1 & atoms2
    unique1 = atoms1 - atoms2
    unique2 = atoms2 - atoms1
    total_unique_atoms = atoms1 | atoms2

    len1 = count_atom_lines(file1_lines)
    len2 = count_atom_lines(file2_lines)
    max_lines = max(len1, len2)

    common_count = len(common_atoms)
    similarity = round((common_count / max_lines * 100), 4) if max_lines > 0 else 0.0
    jaccard = round((common_count / len(total_unique_atoms) * 100), 4) if total_unique_atoms else 0.0

    min_count = min(len(atoms1), len(atoms2)) if (atoms1 and atoms2) else 0
    overlap_coeff = round((common_count / min_count * 100), 4) if min_count > 0 else 0.0

    metrics = {
        "Total_Common_Atoms": common_count,
        "Unique_to_File1": len(unique1),
        "Unique_to_File2": len(unique2),
        "Total_Unique_Atoms": len(total_unique_atoms),
        "Lines_File1": len1,
        "Lines_File2": len2,
        "Max_Lines": max_lines,
        "Percentage_Similarity": similarity,
        "Jaccard_Index": jaccard,
        "Overlap_Coefficient": overlap_coeff
    }

    return metrics, common_atoms


def save_common_atoms(common_atoms: Set[Any], file1_rel: str, file2_rel: str, output_dir: Path) -> Path:
    """
    Saves common atoms list to a formatted text file in output_dir / 'common_atom_lists'.
    """
    out_dir = Path(output_dir) / "common_atom_lists"
    out_dir.mkdir(parents=True, exist_ok=True)

    safe1 = file1_rel.replace("/", "_").replace("\\", "_").replace(".", "_")
    safe2 = file2_rel.replace("/", "_").replace("\\", "_").replace(".", "_")
    out_filename = f"common_atoms_{safe1}_vs_{safe2}.txt"
    out_path = out_dir / out_filename

    with out_path.open("w", encoding="utf-8") as out:
        out.write("Atom_Name\tResidue_Number\n")
        for atom in sorted(common_atoms):
            if isinstance(atom, tuple) and len(atom) == 2:
                out.write(f"{atom[0]}\t{atom[1]}\n")
            elif isinstance(atom, tuple) and len(atom) == 4:
                out.write(f"{atom[3]}\t{atom[2]}\t{atom[1]}\t{atom[0]}\n")
            else:
                out.write(f"{atom}\n")

    return out_path


# ==============================================================================
# Auto-Discovery & Pairwise Processing Engines
# ==============================================================================

def auto_discover_int_pairs(base_dir: Path) -> List[Tuple[Path, Path, str, str]]:
    """
    Auto-discovers .int file pairs within PDB subdirectories under base_dir.
    Returns: list of (full_path1, full_path2, rel_path1, rel_path2)
    """
    base_dir = Path(base_dir).resolve()
    pairs = []

    # Find all .int files under base_dir
    int_files = sorted(list(base_dir.glob("**/*.int")))

    # Group by PDB ID (first token before '_')
    pdb_groups: Dict[str, List[Path]] = {}
    for f in int_files:
        pdb_id = f.name.split(".")[0].split("_")[0].upper()
        pdb_groups.setdefault(pdb_id, []).append(f)

    for pdb_id, files in pdb_groups.items():
        if len(files) >= 2:
            # Pair consecutive or all combination files
            for i in range(len(files) - 1):
                f1 = files[i]
                f2 = files[i + 1]
                try:
                    rel1 = str(f1.relative_to(base_dir))
                    rel2 = str(f2.relative_to(base_dir))
                except ValueError:
                    rel1 = f1.name
                    rel2 = f2.name
                pairs.append((f1, f2, rel1, rel2))

    return pairs


def process_interface_pairs(
    pairs_list: List[Tuple[Path, Path, str, str]],
    output_dir: Path,
    mode: str = "atom_resnum",
    save_atom_lists: bool = True
) -> pd.DataFrame:
    """
    Processes a list of interface file pairs, computes set metrics, and returns summary DataFrame.
    """
    summary_records = []

    for idx, (f1_path, f2_path, f1_rel, f2_rel) in enumerate(pairs_list, 1):
        if not f1_path.exists():
            logger.warning(f"File 1 not found: {f1_path} (skipping pair)")
            continue
        if not f2_path.exists():
            logger.warning(f"File 2 not found: {f2_path} (skipping pair)")
            continue

        try:
            lines1 = [l.strip() for l in f1_path.read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()]
            lines2 = [l.strip() for l in f2_path.read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()]

            metrics, common_atoms = compare_int_files(lines1, lines2, mode=mode)

            atom_list_file = ""
            if save_atom_lists and common_atoms:
                saved_path = save_common_atoms(common_atoms, f1_rel, f2_rel, output_dir)
                atom_list_file = str(saved_path.name)

            record = {
                "Pair_Index": idx,
                "File1_Path": f1_rel,
                "File2_Path": f2_rel,
                "File1_Name": f1_path.name,
                "File2_Name": f2_path.name,
                **metrics,
                "Common_Atoms_File": atom_list_file
            }
            summary_records.append(record)

        except Exception as e:
            logger.error(f"Error processing pair ({f1_path.name}, {f2_path.name}): {e}", exc_info=True)

    return pd.DataFrame(summary_records)


# ==============================================================================
# CLI Entrypoint
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Unified Common & Unique Interface Analysis Suite for Biomolecular Complex Interfaces."
    )
    parser.add_argument(
        "-l", "--input_list", type=Path, default=None,
        help="Path to text file containing paired relative .int paths."
    )
    parser.add_argument(
        "-b", "--base_dir", type=Path, default=DEFAULT_BASE_DIR,
        help=f"Base directory containing interface output folders (default: {DEFAULT_BASE_DIR})."
    )
    parser.add_argument(
        "-o", "--output_dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help=f"Target output directory for summary files (default: {DEFAULT_OUTPUT_DIR})."
    )
    parser.add_argument(
        "--mode", type=str, choices=["atom_resnum", "full"], default="atom_resnum",
        help="Atom matching mode: 'atom_resnum' ((atom_name, res_num), default) or 'full' ((chain_id, res_name, res_num, atom_name))."
    )
    parser.add_argument(
        "--auto_discover", action="store_true",
        help="Automatically discover and pair .int files within PDB folders."
    )
    parser.add_argument(
        "--skip_atom_lists", action="store_true",
        help="Skip exporting individual common atom text files."
    )

    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("TF-NRD COMMON & UNIQUE INTERFACE ANALYSIS SUITE")
    logger.info("=" * 60)
    logger.info(f"Base Directory : {base_dir}")
    logger.info(f"Output Dir     : {out_dir}")
    logger.info(f"Atom Mode      : {args.mode}")

    pairs_to_process = []

    # 1. Input list mode if provided or found
    input_list = Path(args.input_list) if args.input_list else None
    if not input_list and not args.auto_discover:
        if DEFAULT_INPUT_LIST.exists():
            input_list = DEFAULT_INPUT_LIST
        elif (base_dir / "generated_int_paths.txt").exists():
            input_list = base_dir / "generated_int_paths.txt"

    if input_list and input_list.exists() and not args.auto_discover:
        logger.info(f"Reading paired interface paths from list: {input_list}")
        raw_paths = [line.strip() for line in input_list.read_text(encoding="utf-8").splitlines() if line.strip()]
        i = 0
        while i < len(raw_paths):
            p1_rel = raw_paths[i]
            if i + 1 < len(raw_paths):
                p2_rel = raw_paths[i + 1]
                full1 = base_dir / p1_rel if not Path(p1_rel).is_absolute() else Path(p1_rel)
                full2 = base_dir / p2_rel if not Path(p2_rel).is_absolute() else Path(p2_rel)

                pairs_to_process.append((full1, full2, p1_rel, p2_rel))
            i += 2
    else:
        logger.info("Auto-discovering .int interface file pairs across directory tree...")
        pairs_to_process = auto_discover_int_pairs(base_dir)

    logger.info(f"Found {len(pairs_to_process)} interface pairs to compare.")

    if not pairs_to_process:
        logger.warning("No valid interface pairs found for comparison.")
        return

    # Process interface pairs
    df_summary = process_interface_pairs(
        pairs_to_process,
        output_dir=out_dir,
        mode=args.mode,
        save_atom_lists=not args.skip_atom_lists
    )

    # Save summary outputs
    tsv_out = out_dir / "atom_line_common_summary.tsv"
    csv_out = out_dir / "atom_line_common_summary.csv"
    json_out = out_dir / "atom_line_common_summary.json"
    xlsx_out = out_dir / "atom_line_common_summary.xlsx"

    df_summary.to_csv(tsv_out, sep="\t", index=False)
    df_summary.to_csv(csv_out, index=False)
    df_summary.to_json(json_out, orient="records", indent=2)

    logger.info(f"Exported TSV summary ({len(df_summary)} pairs) -> {tsv_out}")
    logger.info(f"Exported CSV summary -> {csv_out}")
    logger.info(f"Exported JSON summary -> {json_out}")

    try:
        with pd.ExcelWriter(xlsx_out, engine="openpyxl") as writer:
            df_summary.to_excel(writer, sheet_name="Common_Unique_Atoms", index=False)
        logger.info(f"Exported Excel summary workbook -> {xlsx_out}")
    except Exception as e:
        logger.warning(f"Could not export standalone Excel {xlsx_out}: {e}")

    # Update Table S4 in Supplementary Excel workbooks if present
    for supp_path in [out_dir / "Supplementary_Tables.xlsx", out_dir / "Supplementary.xlsx", PROJECT_ROOT / "results" / "Supplementary.xlsx"]:
        if supp_path.exists():
            try:
                with pd.ExcelWriter(supp_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
                    df_summary.to_excel(writer, sheet_name="Table S4", index=False)
                logger.info(f"Updated sheet 'Table S4' in Supplementary workbook: {supp_path}")
            except Exception as e:
                logger.warning(f"Could not update sheet 'Table S4' in {supp_path}: {e}")

    logger.info("COMMON & UNIQUE INTERFACE ANALYSIS COMPLETE")


if __name__ == "__main__":
    main()
