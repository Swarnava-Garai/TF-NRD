#!/usr/bin/env python3
"""
cross_validate_interface.py
===========================
Unified interface cross-validation and quality control framework for biomolecular complexes.

Combines three essential structural checks for interface (.int) datasets:
  1. Insertion Code Detection (formerly check_insertion_codes.sh):
     Scans ATOM lines for PDB insertion codes (column 27) that could cause silent
     residue-number merging during atom-matching.
  2. Residue Numbering Consistency (formerly check_residue_numbering.py):
     Verifies if residue numbers map to identical amino acid types across all chain
     copies within each PDB entry.
  3. Chain Correspondence & Register Shift Diagnosis (formerly check_chain_correspondence.py):
     Diagnoses residue mismatches between chain pairs to distinguish register shifts
     (offset numbering) from genuine biological mismatches (hetero-oligomers / distinct subunits).

Usage:
  # Run all validation checks on an interface directory:
  python3 cross_validate_interface.py /path/to/interface_root

  # Run specific checks:
  python3 cross_validate_interface.py /path/to/interface_root --check insertion
  python3 cross_validate_interface.py /path/to/interface_root --check numbering
  python3 cross_validate_interface.py /path/to/interface_root --check correspondence

  # Save report to file:
  python3 cross_validate_interface.py /path/to/interface_root -o validation_report.txt
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional

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

# Nucleotide residue codes (DNA + RNA) to exclude from protein sequence matching
NUCLEOTIDES = {
    "DA", "DC", "DG", "DT", "DU", "DI",
    "A", "C", "G", "U", "I", "T",
    "RA", "RC", "RG", "RU"
}


# ==============================================================================
# Helper Functions & File Discovery
# ==============================================================================

def get_pdb_id_from_dirname(dirname: str) -> str:
    """Extracts PDB ID as the first underscore-delimited token from directory name."""
    return os.path.basename(dirname).split("_")[0].upper()


def discover_interface_groups(base_dir: Path) -> Dict[str, List[Path]]:
    """
    Discovers all subdirectories containing .int files and groups them by PDB ID.
    If .int files exist directly in base_dir, groups by PDB ID prefix of filenames.
    """
    base_dir = Path(base_dir).resolve()
    groups = defaultdict(list)

    # 1. Search subdirectories
    subdirs = [p for p in base_dir.iterdir() if p.is_dir()]
    for sdir in sorted(subdirs):
        int_files = list(sdir.glob("*.int"))
        if int_files:
            pdb_id = get_pdb_id_from_dirname(sdir.name)
            groups[pdb_id].append(sdir)

    # 2. Check if .int files are located directly in base_dir
    direct_ints = list(base_dir.glob("*.int"))
    if direct_ints and not groups:
        for f in direct_ints:
            pdb_id = f.name.split("_")[0].split(".")[0].upper()
            groups[pdb_id].append(base_dir)

    return groups


def get_all_int_files(base_dir: Path) -> List[Path]:
    """Recursively collects all .int files under base_dir."""
    return sorted(list(Path(base_dir).glob("**/*.int")))


# ==============================================================================
# Check 1: Insertion Code Detection
# ==============================================================================

def check_insertion_codes(base_dir: Path) -> Tuple[int, Dict[Path, List[str]]]:
    """
    Scans all .int files for PDB insertion codes in column 27 (1-indexed).

    Returns:
        Tuple of (total_flagged_files_count, flagged_files_dict).
    """
    int_files = get_all_int_files(base_dir)
    flagged = {}

    for fpath in int_files:
        codes = []
        with open(fpath, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("ATOM") and len(line) >= 27:
                    # Column 27 (1-indexed) is index 26
                    ins_code = line[26].strip()
                    if ins_code:
                        codes.append(ins_code)
        if codes:
            flagged[fpath] = codes

    return len(flagged), flagged


# ==============================================================================
# Check 2: Residue Numbering Consistency
# ==============================================================================

def parse_atoms_from_int(filepath: Path):
    """Yields (chain, res_num, res_name) for each ATOM record in a .int file."""
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM") or len(line) < 26:
                continue
            res_name = line[17:20].strip()
            chain = line[21].strip() if len(line) > 21 else "A"
            res_num = line[22:26].strip()
            if res_num and res_name:
                yield chain, res_num, res_name


def check_residue_numbering(groups: Dict[str, List[Path]]) -> Tuple[int, int, Dict[str, Dict]]:
    """
    Cross-checks that residue numbers map to the same amino acid type across all chain copies.

    Returns:
        Tuple of (total_flagged_pdbs, total_flagged_residues, details_dict).
    """
    total_flagged_pdbs = 0
    total_flagged_residues = 0
    flagged_details = {}

    for pdb_id, dirs in sorted(groups.items()):
        # res_num -> res_name -> set of (chain, filepath)
        res_map = defaultdict(lambda: defaultdict(set))

        for d in dirs:
            int_files = list(d.glob("*.int")) if d.is_dir() else [d]
            for fpath in int_files:
                for chain, res_num, res_name in parse_atoms_from_int(fpath):
                    if res_name in NUCLEOTIDES:
                        continue  # Skip nucleic acids (complementary strands differ naturally)
                    res_map[res_num][res_name].add((chain, str(fpath)))

        inconsistent = {
            res_num: names
            for res_num, names in res_map.items()
            if len(names) > 1
        }

        if inconsistent:
            total_flagged_pdbs += 1
            total_flagged_residues += len(inconsistent)
            flagged_details[pdb_id] = inconsistent

    return total_flagged_pdbs, total_flagged_residues, flagged_details


# ==============================================================================
# Check 3: Chain Correspondence & Register Shift Diagnosis
# ==============================================================================

def parse_chain_residues_for_alignment(filepath: Path) -> Dict[str, Dict[int, str]]:
    """Returns {chain: {res_num (int): res_name}} for protein residues in a .int file."""
    chains = defaultdict(dict)
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM") or len(line) < 26:
                continue
            res_name = line[17:20].strip()
            chain = line[21].strip() if len(line) > 21 else "A"
            res_num_str = line[22:26].strip()

            if res_name in NUCLEOTIDES or not res_num_str:
                continue
            try:
                res_num = int(res_num_str)
                chains[chain][res_num] = res_name
            except ValueError:
                continue
    return chains


def compute_best_offset_match(
    res_a: Dict[int, str],
    res_b: Dict[int, str],
    max_shift: int = 30
) -> Tuple[int, float, int, float]:
    """
    Evaluates sequence identity across offsets [-max_shift..+max_shift].

    Returns:
        (best_offset, best_fraction, overlap_at_0, frac_at_0)
    """
    nums_a = set(res_a.keys())
    nums_b = set(res_b.keys())
    overlap0 = nums_a & nums_b
    match0 = sum(1 for n in overlap0 if res_a[n] == res_b[n])
    frac0 = (match0 / len(overlap0)) if overlap0 else 0.0

    best_offset = 0
    best_frac = frac0

    for shift in range(-max_shift, max_shift + 1):
        if shift == 0:
            continue
        overlap = 0
        match = 0
        for n, name in res_a.items():
            shifted = n + shift
            if shifted in res_b:
                overlap += 1
                if res_b[shifted] == name:
                    match += 1
        if overlap == 0:
            continue
        frac = match / overlap
        min_required = max(3, int(0.5 * len(overlap0))) if overlap0 else 3
        if overlap >= min_required and frac > best_frac:
            best_offset = shift
            best_frac = frac

    return best_offset, best_frac, len(overlap0), frac0


def check_chain_correspondence(
    groups: Dict[str, List[Path]],
    max_shift: int = 30
) -> Dict[str, List[Dict]]:
    """
    Evaluates all protein chain pairs within each PDB entry for register shifts vs genuine mismatches.

    Returns:
        Dict mapping pdb_id -> list of evaluated pair records.
    """
    results = {}

    for pdb_id, dirs in sorted(groups.items()):
        chain_residues = defaultdict(dict)
        for d in dirs:
            int_files = list(d.glob("*.int")) if d.is_dir() else [d]
            for fpath in int_files:
                for chain, resmap in parse_chain_residues_for_alignment(fpath).items():
                    chain_residues[chain].update(resmap)

        chains = sorted(chain_residues.keys())
        if len(chains) < 2:
            continue

        pair_diagnostics = []
        for i in range(len(chains)):
            for j in range(i + 1, len(chains)):
                ca, cb = chains[i], chains[j]
                offset, best_frac, overlap0, frac0 = compute_best_offset_match(
                    chain_residues[ca], chain_residues[cb], max_shift=max_shift
                )

                if overlap0 == 0 or frac0 >= 0.999:
                    continue  # Perfectly consistent or no overlap

                if overlap0 < 5:
                    verdict = "INSUFFICIENT DATA (too few shared residues)"
                elif offset != 0 and best_frac >= 0.9:
                    verdict = "REGISTER SHIFT"
                elif best_frac < 0.5:
                    verdict = "GENUINE MISMATCH (likely different subunits)"
                else:
                    verdict = "PARTIAL / UNCLEAR"

                pair_diagnostics.append({
                    "chain_a": ca,
                    "chain_b": cb,
                    "overlap0": overlap0,
                    "match0": frac0,
                    "best_offset": offset,
                    "best_match": best_frac,
                    "verdict": verdict
                })

        if pair_diagnostics:
            results[pdb_id] = pair_diagnostics

    return results


# ==============================================================================
# Comprehensive Execution & Reporting
# ==============================================================================

def run_cross_validation(
    base_dir: Path,
    checks_to_run: Set[str],
    max_shift: int = 30,
    output_report: Optional[Path] = None
) -> None:
    """Executes selected validation modules and prints/saves structured report."""
    base_dir = Path(base_dir).resolve()
    report_lines = []

    def out(line: str = ""):
        print(line)
        report_lines.append(line)

    out("=" * 80)
    out(" BIOMOLECULAR COMPLEX INTERFACE CROSS-VALIDATION REPORT")
    out("=" * 80)
    out(f"Target Directory : {base_dir}")

    groups = discover_interface_groups(base_dir)
    all_int_files = get_all_int_files(base_dir)

    out(f"Total .int Files : {len(all_int_files)}")
    out(f"PDB ID Groups    : {len(groups)}")
    out("=" * 80)

    # 1. Insertion Code Check
    if "all" in checks_to_run or "insertion" in checks_to_run:
        out("\n[CHECK 1] PDB Insertion Codes (Column 27)")
        out("-" * 80)
        n_flagged, ins_flagged = check_insertion_codes(base_dir)
        if n_flagged == 0:
            out("✓ PASS: No insertion codes found in any .int file.")
            out("  Residue sequence numbers are clean (no silent merging hazard).")
        else:
            out(f"⚠ WARNING: Found insertion codes in {n_flagged} file(s):")
            for fpath, codes in sorted(ins_flagged.items()):
                unique_codes = sorted(list(set(codes)))
                rel_path = fpath.relative_to(base_dir) if fpath.is_relative_to(base_dir) else fpath
                out(f"  - {rel_path} : {len(codes)} line(s) with code(s) {unique_codes}")
            out("  Action: Ensure atom-matching keys include insertion code (atom, res_num, ins_code).")

    # 2. Residue Numbering Consistency Check
    if "all" in checks_to_run or "numbering" in checks_to_run:
        out("\n[CHECK 2] Residue Numbering Consistency Across Chain Copies")
        out("-" * 80)
        n_pdbs, n_res, num_flagged = check_residue_numbering(groups)
        if n_pdbs == 0:
            out("✓ PASS: Residue numbering is 100% consistent across all protein chain copies.")
            out("  Every residue number refers to the same amino acid type across copies.")
        else:
            out(f"⚠ WARNING: {n_res} inconsistent residue number(s) across {n_pdbs} PDB ID group(s):")
            for pdb_id, inc_map in sorted(num_flagged.items()):
                out(f"  [FLAGGED] PDB {pdb_id}: {len(inc_map)} inconsistent residue(s)")
                for rnum, names in sorted(inc_map.items(), key=lambda x: int(x[0]) if x[0].isdigit() else 0):
                    details = []
                    for rname, occs in sorted(names.items()):
                        occ_str = ", ".join(f"{c}:{Path(fp).parent.name}/{Path(fp).name}" for c, fp in sorted(occs))
                        details.append(f"{rname} ({occ_str})")
                    out(f"    - res {rnum:>4s}: {'; '.join(details)}")
            out("  Action: Review flagged entries; residue-number matching may compare non-homologous positions.")

    # 3. Chain Correspondence & Register Shift Diagnosis
    if "all" in checks_to_run or "correspondence" in checks_to_run:
        out("\n[CHECK 3] Chain Correspondence & Register Shift Diagnosis")
        out("-" * 80)
        corr_results = check_chain_correspondence(groups, max_shift=max_shift)
        if not corr_results:
            out("✓ PASS: All multi-chain protein pairs agree at offset 0 (no register shifts or mismatches).")
        else:
            out(f"⚠ DIAGNOSTICS: Identified chain-pair numbering discrepancies in {len(corr_results)} PDB ID group(s):")
            for pdb_id, pairs in sorted(corr_results.items()):
                out(f"  PDB {pdb_id}:")
                for p in pairs:
                    ca, cb = p["chain_a"], p["chain_b"]
                    ov0, m0 = p["overlap0"], p["match0"]
                    bo, bm = p["best_offset"], p["best_match"]
                    verd = p["verdict"]
                    out(f"    - Chains {ca} vs {cb}: overlap={ov0}, match@0={m0:.1%} | "
                        f"best_offset={bo:+d} -> match={bm:.1%} [{verd}]")

    out("\n" + "=" * 80)
    out(" SUMMARY & RECOMMENDATIONS")
    out("=" * 80)
    out("1. If all checks PASS: The interface identity matching key (atom_name, res_num) is safe.")
    out("2. If REGISTER SHIFT is flagged: Align numbering by the detected offset or via sequence alignment.")
    out("3. If GENUINE MISMATCH is flagged: Chains are distinct proteins; do not treat them as identical copies.")
    out("=" * 80)

    # Save to report file if requested
    if output_report:
        output_report = Path(output_report).resolve()
        output_report.parent.mkdir(parents=True, exist_ok=True)
        with open(output_report, "w", encoding="utf-8") as f:
            f.write("\n".join(report_lines) + "\n")
        print(f"\n✓ Saved validation report to: {output_report}")


# ==============================================================================
# CLI Entrypoint
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Unified Interface Cross-Validation Pipeline for Biomolecular Complex Datasets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "target_dir",
        nargs="?",
        default=".",
        help="Path to interface directory containing .int files or subdirectories."
    )
    parser.add_argument(
        "-d", "--dir",
        dest="dir_flag",
        help="Explicit path to interface directory."
    )
    parser.add_argument(
        "-c", "--check",
        choices=["all", "insertion", "numbering", "correspondence"],
        default="all",
        help="Specify which validation check to execute."
    )
    parser.add_argument(
        "--max-shift",
        type=int,
        default=30,
        help="Maximum residue register offset search range [-N..+N]."
    )
    parser.add_argument(
        "-o", "--output-report",
        dest="output_report",
        help="Optional path to save text report."
    )

    return parser.parse_args()


def main():
    args = parse_args()
    target_path = Path(args.dir_flag or args.target_dir).resolve()

    if not target_path.exists() or not target_path.is_dir():
        print(f"Error: Target path '{target_path}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    checks = {args.check}
    run_cross_validation(
        base_dir=target_path,
        checks_to_run=checks,
        max_shift=args.max_shift,
        output_report=args.output_report
    )


if __name__ == "__main__":
    main()
