#!/usr/bin/env python3
"""
blast2nr.py
===========
Filters sequence databases to generate Non-Redundant (NR) datasets from BLAST tabular output.

Formula:
  Score = 0.01 * Identity (%) * Coverage (%)
  If Score >= Cutoff, subject sequence is marked as redundant to the query sequence.

Key Features:
  - Efficient O(1) set-based redundancy clustering and FASTA filtering.
  - Robust pathlib-based path resolution (resolves TF-NRD project structure).
  - Defaults configured to TF-NRD/input_data/blast and TF-NRD/results/blast.
  - Dynamic FASTA naming embedding the exact count of selected non-redundant sequences.
  - CLI supporting both positional syntax (legacy compatibility) and standard flags.
  - Comprehensive logging and reporting of non-redundant and rejected sequences.

Usage:
  # 1. Run with default project paths (input_data/blast -> results/blast):
  python3 blast2nr.py

  # 2. Run with custom cutoff threshold:
  python3 blast2nr.py -c 25.0

  # 3. Custom files and output directory:
  python3 blast2nr.py -b blast.txt -f input.fasta -c 35 -o results/blast/nr_output.fa -l results/blast/classified.txt

  # 4. Legacy positional syntax:
  python3 blast2nr.py <blast_output> <cutoff> <fasta_file>
"""

import re
import sys
import logging
import argparse
from pathlib import Path
from typing import Set, Tuple, List, Optional
from Bio import SeqIO

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
PROJECT_ROOT = SCRIPT_DIR.parent
LOG_FILE = SCRIPT_PATH.with_suffix(".log")

# Configure logger to stdout and script-named log file
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

DEFAULT_INPUT_DIR = PROJECT_ROOT / "input_data" / "blast"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "blast"

DEFAULT_BLAST_FILE = DEFAULT_INPUT_DIR / "blastoutput_fasta_prot.txt"
DEFAULT_FASTA_FILE = DEFAULT_INPUT_DIR / "List_protein_seq_498_complexes.fa"
DEFAULT_CUTOFF = 35.0


def parse_blast_for_non_redundant(
    blast_file_path: Path,
    cutoff: float,
    query_col: int = 0,
    subject_col: int = 1,
    identity_col: int = 2,
    coverage_col: int = 12,
    truncate_id_len: Optional[int] = None
    ) -> Tuple[List[str], List[str]]:
    """
    Parses BLAST tabular data to identify and segregate redundant sequence IDs.

    Args:
        blast_file_path: Path to the tabular BLAST output file.
        cutoff: Float cutoff threshold (e.g. 35.0).
        query_col: 0-indexed column for query ID (default 0).
        subject_col: 0-indexed column for subject ID (default 1).
        identity_col: 0-indexed column for sequence identity percentage (default 2).
        coverage_col: 0-indexed column for alignment coverage percentage (default 12).
        truncate_id_len: Optional integer to slice IDs (e.g. 6). Defaults to None (full ID).

    Returns:
        Tuple of (rejected_list, selected_list).
    """
    blast_file_path = Path(blast_file_path).resolve()
    if not blast_file_path.exists():
        raise FileNotFoundError(f"BLAST input file not found: {blast_file_path}")

    max_col_needed = max(query_col, subject_col, identity_col, coverage_col)
    rejected_set: Set[str] = set()
    all_seen: Set[str] = set()

    with open(blast_file_path, "r", encoding="utf-8", errors="replace") as infile:
        for line in infile:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue

            columns = line_str.split("\t")
            if len(columns) <= max_col_needed:
                continue

            q_id = columns[query_col].strip()
            s_id = columns[subject_col].strip()

            if truncate_id_len:
                q_id = q_id[:truncate_id_len]
                s_id = s_id[:truncate_id_len]

            all_seen.add(q_id)
            all_seen.add(s_id)

            # Skip self-hits
            if q_id == s_id:
                continue

            # If query itself has already been rejected, it cannot serve as representative
            if q_id in rejected_set:
                continue

            # Compute length-weighted identity-coverage score: (0.01 * Identity * Coverage)
            try:
                identity = float(columns[identity_col].strip())
                coverage = float(columns[coverage_col].strip())
                score = 0.01 * identity * coverage
            except (ValueError, IndexError):
                continue

            # If score breaches cutoff, mark subject as redundant
            if score >= cutoff:
                rejected_set.add(s_id)

    # Non-redundant representative set
    selected_set = all_seen - rejected_set

    # Return sorted lists for deterministic ordering
    return sorted(list(rejected_set)), sorted(list(selected_set))


def filter_fasta_sequences(
    fasta_input_path: Path,
    selected_ids: Set[str],
    output_fasta_path: Path
    ) -> int:
    """
    Filters FASTA records, exporting only non-redundant representatives.

    Args:
        fasta_input_path: Path to the input FASTA file.
        selected_ids: Set of non-redundant sequence IDs to retain.
        output_fasta_path: Path to write the output FASTA file.

    Returns:
        Number of sequences written.
    """
    fasta_input_path = Path(fasta_input_path).resolve()
    output_fasta_path = Path(output_fasta_path).resolve()
    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    written_count = 0
    with open(output_fasta_path, "w", encoding="utf-8") as out_f:
        for record in SeqIO.parse(str(fasta_input_path), "fasta"):
            rec_id = record.id.strip()
            if rec_id in selected_ids:
                out_f.write(f">{record.description}\n{str(record.seq)}\n")
                written_count += 1

    return written_count


def write_classification_report(
    output_list_path: Path,
    selected_list: List[str],
    rejected_list: List[str]
    ) -> None:
    """Writes classified PDB / sequence IDs report."""
    output_list_path = Path(output_list_path).resolve()
    output_list_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_list_path, "w", encoding="utf-8") as out_f:
        out_f.write(f"# Selected Non-Redundant IDs ({len(selected_list)} total)\n")
        for s_id in selected_list:
            out_f.write(f"{s_id}\n")

        out_f.write(f"\n# Rejected Redundant IDs ({len(rejected_list)} total)\n")
        for r_id in rejected_list:
            out_f.write(f"{r_id}\n")


def process_blast2nr(
    blast_path: Path = DEFAULT_BLAST_FILE,
    fasta_path: Path = DEFAULT_FASTA_FILE,
    cutoff: float = DEFAULT_CUTOFF,
    output_fasta_path: Optional[Path] = None,
    output_list_path: Optional[Path] = None,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    query_col: int = 0,
    subject_col: int = 1,
    identity_col: int = 2,
    coverage_col: int = 12,
    truncate_id_len: Optional[int] = None
    ) -> Tuple[List[str], List[str], Path, Path]:
    """
    Executes end-to-end BLAST non-redundancy filtering.

    Returns:
        Tuple of (selected_list, reject_list, final_out_fasta, final_out_list).
    """
    blast_path = Path(blast_path).resolve()
    fasta_path = Path(fasta_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not blast_path.exists():
        raise FileNotFoundError(f"BLAST output file not found: {blast_path}")
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA input file not found: {fasta_path}")

    logger.info(f"Processing BLAST data from: {blast_path}")
    logger.info(f"Cutoff Threshold: {cutoff}%")
    logger.info(f"Input FASTA: {fasta_path}")

    # 1. Process BLAST tabular data
    reject_list, selected_list = parse_blast_for_non_redundant(
        blast_file_path=blast_path,
        cutoff=cutoff,
        query_col=query_col,
        subject_col=subject_col,
        identity_col=identity_col,
        coverage_col=coverage_col,
        truncate_id_len=truncate_id_len
    )

    selected_count = len(selected_list)
    logger.info(f"Classified {selected_count} Non-Redundant and {len(reject_list)} Redundant sequence IDs.")

    # 2. Determine output filenames
    if output_fasta_path:
        final_out_fasta = Path(output_fasta_path).resolve()
    else:
        stem = fasta_path.stem
        if re.search(r"\d+", stem):
            new_stem = re.sub(r"\d+", str(selected_count), stem)
            if not new_stem.startswith("nr_"):
                new_stem = f"nr_{new_stem}"
            out_fasta_name = f"{new_stem}{fasta_path.suffix}"
        else:
            out_fasta_name = f"nr_{selected_count}_{fasta_path.name}"
        final_out_fasta = output_dir / out_fasta_name

    if output_list_path:
        final_out_list = Path(output_list_path).resolve()
    else:
        final_out_list = output_dir / "blastOutputClassifiedList"

    # 3. Write classification reports
    write_classification_report(final_out_list, selected_list, reject_list)
    logger.info(f"Saved classification summary -> {final_out_list}")

    # 4. Filter FASTA records
    selected_set = set(selected_list)
    written_count = filter_fasta_sequences(fasta_path, selected_set, final_out_fasta)
    logger.info(f"Saved {written_count} non-redundant sequences -> {final_out_fasta}")

    return selected_list, reject_list, final_out_fasta, final_out_list


def make_nrFromBlastOutput(blast_file_path, cutoff_value):
    """Legacy backward compatibility wrapper."""
    return parse_blast_for_non_redundant(Path(blast_file_path), float(cutoff_value))


def parse_args():
    """Parses command-line arguments supporting both positional and optional flags."""
    parser = argparse.ArgumentParser(
        description="Filter redundant sequences from BLAST tabular output and FASTA database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Positional or flagged arguments
    parser.add_argument("blast_pos", nargs="?", default=None, help="Path to BLAST tabular output file")
    parser.add_argument("cutoff_pos", nargs="?", type=float, default=None, help="Score cutoff threshold (e.g. 35)")
    parser.add_argument("fasta_pos", nargs="?", default=None, help="Path to input FASTA file")

    parser.add_argument("-b", "--blast", dest="blast_flag", default=str(DEFAULT_BLAST_FILE), help="Path to BLAST tabular output file")
    parser.add_argument("-c", "--cutoff", dest="cutoff_flag", type=float, default=DEFAULT_CUTOFF, help="Cutoff percentage threshold")
    parser.add_argument("-f", "--fasta", dest="fasta_flag", default=str(DEFAULT_FASTA_FILE), help="Path to input FASTA file")
    parser.add_argument("-d", "--output-dir", dest="output_dir", default=str(DEFAULT_OUTPUT_DIR), help="Output directory for results")
    parser.add_argument("-o", "--output-fasta", dest="output_fasta", default=None, help="Path to output non-redundant FASTA")
    parser.add_argument("-l", "--output-list", dest="output_list", default=None, help="Path to output classified ID list")

    parser.add_argument("--query-col", type=int, default=0, help="0-based column index for Query ID")
    parser.add_argument("--subject-col", type=int, default=1, help="0-based column index for Subject ID")
    parser.add_argument("--identity-col", type=int, default=2, help="0-based column index for Identity percentage")
    parser.add_argument("--coverage-col", type=int, default=12, help="0-based column index for Coverage percentage")
    parser.add_argument("--truncate-len", type=int, default=None, help="Optional length to slice ID strings (e.g. 6)")

    return parser.parse_args()


def main():
    args = parse_args()

    # Determine paths prioritizing positional if provided, otherwise flags / defaults
    blast_input = args.blast_pos if args.blast_pos else args.blast_flag
    cutoff_val = args.cutoff_pos if args.cutoff_pos is not None else args.cutoff_flag
    fasta_input = args.fasta_pos if args.fasta_pos else args.fasta_flag

    blast_path = Path(blast_input)
    fasta_path = Path(fasta_input)
    output_dir = Path(args.output_dir)

    try:
        process_blast2nr(
            blast_path=blast_path,
            fasta_path=fasta_path,
            cutoff=cutoff_val,
            output_fasta_path=Path(args.output_fasta) if args.output_fasta else None,
            output_list_path=Path(args.output_list) if args.output_list else None,
            output_dir=output_dir,
            query_col=args.query_col,
            subject_col=args.subject_col,
            identity_col=args.identity_col,
            coverage_col=args.coverage_col,
            truncate_id_len=args.truncate_len
        )
    except Exception as e:
        logger.error(f"Error during BLAST non-redundancy processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
