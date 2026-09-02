#!/usr/bin/env python3
"""
mapper_TF_pdb_Pfam_Rfam.py
==========================
Maps target transcription factor PDB entries against domain and RNA family databases:
  1. Pfam Database ('pdb_pfam_mapping.txt'): Protein domain mapping. User needs to download the mapping file (https://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt).
  2. Rfam Database ('Rfam.pdb'): RNA family / structural motif mapping. User needs to download the Rfam file (https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.pdb.gz).

Features:
  - Common input IDs file ('final_TF_ids_377.txt') used across both Pfam and Rfam mappings.
  - Fast single-pass O(N) streaming parsing for large (70MB+) mapping files.
  - Supports both Pfam (protein domain) and Rfam (RNA motif) mapping modes.
  - Generates standard TSV, formatted readable text, Excel (.xlsx), and JSON outputs.
  - Configurable CLI with argparse.

Usage:
  # 1. Run default Pfam protein domain mapping (377 TFs):
  python3 mapper_TF_pdb_Pfam_Rfam.py

  # 2. Run Rfam RNA domain mapping:
  python3 mapper_TF_pdb_Pfam_Rfam.py --mode rfam

  # 3. Run both Pfam and Rfam mapping pipelines:
  python3 mapper_TF_pdb_Pfam_Rfam.py --mode all

  # 4. Custom input IDs, database, and outputs:
  python3 mapper_TF_pdb_Pfam_Rfam.py --mode rfam -i input_data/domains/final_TF_ids_377.txt -m input_data/domains/Rfam.pdb
"""

import sys
import json
import logging
import argparse
from pathlib import Path
from typing import List, Set, Dict, Optional, Tuple, Any

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

# Common input PDB IDs file for both Pfam and Rfam mapping
DEFAULT_PFAM_IDS_FILE = PROJECT_ROOT / "input_data" / "domains" / "final_TF_ids_377.txt"

# Pfam defaults
DEFAULT_PFAM_MAPPING_FILE = PROJECT_ROOT / "input_data" / "domains" / "pdb_pfam_mapping.txt"
DEFAULT_PFAM_OUTPUT_TSV = PROJECT_ROOT / "results" / "domains" / "377_TF_domain_mapped.tsv"
DEFAULT_PFAM_OUTPUT_XLSX = PROJECT_ROOT / "results" / "domains" / "377_TF_domain_mapped.xlsx"
DEFAULT_PFAM_OUTPUT_JSON = PROJECT_ROOT / "results" / "domains" / "377_TF_domain_mapped.json"
DEFAULT_PFAM_PRETTY_OUTPUT = PROJECT_ROOT / "results" / "domains" / "377_TF_domain_mapped_readable.txt"

# Rfam defaults
DEFAULT_RFAM_MAPPING_FILE = PROJECT_ROOT / "input_data" / "domains" / "Rfam.pdb"
DEFAULT_RFAM_OUTPUT_TSV = PROJECT_ROOT / "results" / "domains" / "55_EM_TF_rfam_mapped.tsv"
DEFAULT_RFAM_OUTPUT_XLSX = PROJECT_ROOT / "results" / "domains" / "55_EM_TF_rfam_mapped.xlsx"
DEFAULT_RFAM_OUTPUT_JSON = PROJECT_ROOT / "results" / "domains" / "55_EM_TF_rfam_mapped.json"
DEFAULT_RFAM_PRETTY_OUTPUT = PROJECT_ROOT / "results" / "domains" / "55_EM_TF_rfam_mapped_readable.txt"

# Schema definitions
PFAM_COLUMNS = [
    "PDB",
    "CHAIN",
    "PDB_START",
    "PDB_END",
    "PFAM_ACCESSION",
    "PFAM_NAME",
    "AUTH_PDBRES_START",
    "AUTH_PDBRES_START_INS_CODE",
    "AUTH_PDBRES_END",
    "AUTH_PDBRES_END_INS_CODE",
    "UNIPROT_ACCESSION",
    "UNP_START",
    "UNP_END",
    "UNP_STR_START",
    "UNP_STR_END"
]

RFAM_COLUMNS = [
    "RFAM_ACCESSION",
    "PDB",
    "CHAIN",
    "PDB_START",
    "PDB_END",
    "BIT_SCORE",
    "E_VALUE",
    "RFAM_START",
    "RFAM_END",
    "HEX_FLAG"
]


# ==============================================================================
# Helper Functions
# ==============================================================================

def load_target_pdb_ids(ids_path: Path) -> Tuple[List[str], Set[str]]:
    """
    Loads target PDB IDs from file.
    Returns (ordered_ids_list, uppercase_ids_set).
    """
    ids_path = Path(ids_path).resolve()
    if not ids_path.exists():
        raise FileNotFoundError(f"PDB IDs file not found: {ids_path}")

    ordered_ids = []
    seen = set()

    with open(ids_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            # Extract first whitespace-delimited token (PDB ID)
            pdb_id = line_str.split()[0].upper()
            if pdb_id not in seen:
                seen.add(pdb_id)
                ordered_ids.append(pdb_id)

    return ordered_ids, seen


def write_formatted_outputs(
    records: List[List[str]],
    columns: List[str],
    output_tsv: Optional[Path] = None,
    output_pretty: Optional[Path] = None,
    output_xlsx: Optional[Path] = None,
    output_json: Optional[Path] = None
    ) -> None:
    """Writes TSV, pretty-aligned text, Excel spreadsheet, and JSON."""
    # 1. TSV Output
    if output_tsv:
        output_tsv = Path(output_tsv).resolve()
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(output_tsv, "w", encoding="utf-8") as f:
            f.write("\t".join(columns) + "\n")
            for row in records:
                f.write("\t".join(row) + "\n")
        logger.info(f"Saved TSV output -> {output_tsv}")

    # 2. Column-Aligned Pretty Output
    if output_pretty:
        output_pretty = Path(output_pretty).resolve()
        output_pretty.parent.mkdir(parents=True, exist_ok=True)

        widths = [len(c) for c in columns]
        for row in records:
            for i, val in enumerate(row):
                if i < len(widths):
                    widths[i] = max(widths[i], len(str(val)))

        header_line = "  ".join(f"{columns[i]:<{widths[i]}}" for i in range(len(columns)))
        separator_line = "  ".join("-" * widths[i] for i in range(len(columns)))

        with open(output_pretty, "w", encoding="utf-8") as f:
            f.write(header_line + "\n")
            f.write(separator_line + "\n")
            for row in records:
                row_line = "  ".join(f"{row[i]:<{widths[i]}}" if i < len(row) else "" for i in range(len(columns)))
                f.write(row_line + "\n")
        logger.info(f"Saved aligned readable text -> {output_pretty}")

    # 3. Excel Output
    if output_xlsx:
        import pandas as pd
        output_xlsx = Path(output_xlsx).resolve()
        output_xlsx.parent.mkdir(parents=True, exist_ok=True)
        df_out = pd.DataFrame(records, columns=columns)
        df_out.to_excel(output_xlsx, index=False)
        logger.info(f"Saved Excel workbook -> {output_xlsx}")

    # 4. JSON Output
    if output_json:
        output_json = Path(output_json).resolve()
        output_json.parent.mkdir(parents=True, exist_ok=True)
        dict_records = [
            {columns[i]: (row[i] if i < len(row) else "") for i in range(len(columns))}
            for row in records
        ]
        with open(output_json, "w", encoding="utf-8") as jf:
            json.dump(dict_records, jf, indent=2, ensure_ascii=False)
        logger.info(f"Saved structured JSON output -> {output_json}")


# ==============================================================================
# Module 1: Pfam Protein Domain Mapping
# ==============================================================================

def parse_pfam_mapping_line(line: str) -> Optional[List[str]]:
    """Parses a line from pdb_pfam_mapping.txt into 15 standard columns."""
    line_str = line.strip()
    if not line_str or line_str.startswith("#"):
        return None

    if "\t" in line_str:
        fields = line_str.split("\t")
        out = [fields[i].strip() if i < len(fields) else "" for i in range(15)]
    else:
        fields = line_str.split()
        if not fields:
            return None
        out = [""] * 15
        for i in range(min(7, len(fields))):
            out[i] = fields[i].strip()

        if len(fields) > 7:
            if fields[7].isdigit():
                out[7] = ""
                out[8] = fields[7]
                out[9] = ""
                out[10:15] = [fields[i].strip() if i < len(fields) else "" for i in range(8, 13)]
            else:
                out[7] = fields[7]
                out[8] = fields[8] if len(fields) > 8 else ""
                if len(fields) > 9 and len(fields[9]) == 1 and fields[9].isalnum() and not fields[9].isdigit():
                    out[9] = fields[9]
                    out[10:15] = [fields[i].strip() if i < len(fields) else "" for i in range(10, 15)]
                else:
                    out[9] = ""
                    out[10:15] = [fields[i].strip() if i < len(fields) else "" for i in range(9, 14)]

    out[0] = out[0].upper()
    return out


def map_pdb_to_pfam(
    ids_path: Path = DEFAULT_PFAM_IDS_FILE,
    mapping_path: Path = DEFAULT_PFAM_MAPPING_FILE,
    output_tsv: Optional[Path] = DEFAULT_PFAM_OUTPUT_TSV,
    output_pretty: Optional[Path] = DEFAULT_PFAM_PRETTY_OUTPUT,
    output_xlsx: Optional[Path] = DEFAULT_PFAM_OUTPUT_XLSX,
    output_json: Optional[Path] = DEFAULT_PFAM_OUTPUT_JSON
    ) -> List[List[str]]:
    """
    Streams the PDB-Pfam database in a single fast pass and exports matching protein domain records.
    """
    ids_path = Path(ids_path).resolve()
    mapping_path = Path(mapping_path).resolve()

    if not mapping_path.exists():
        raise FileNotFoundError(f"PDB-Pfam mapping file not found: {mapping_path}")

    ordered_ids, target_set = load_target_pdb_ids(ids_path)
    logger.info(f"[Pfam Mode] Loaded {len(target_set)} target PDB IDs from: {ids_path}")
    logger.info(f"[Pfam Mode] Scanning PDB-Pfam database: {mapping_path} ({mapping_path.stat().st_size / (1024*1024):.1f} MB)...")

    matched_records = []
    matched_pdbs = set()

    with open(mapping_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue

            first_token = line_str.split(None, 1)[0].split("\t", 1)[0].upper()
            if first_token in target_set:
                parsed = parse_pfam_mapping_line(line_str)
                if parsed and parsed[0] in target_set:
                    matched_records.append(parsed)
                    matched_pdbs.add(parsed[0])

    logger.info(f"[Pfam Mode] Mapping complete: Found {len(matched_records)} domain records for {len(matched_pdbs)} / {len(target_set)} PDB IDs.")

    write_formatted_outputs(
        records=matched_records,
        columns=PFAM_COLUMNS,
        output_tsv=output_tsv,
        output_pretty=output_pretty,
        output_xlsx=output_xlsx,
        output_json=output_json
    )

    return matched_records


# ==============================================================================
# Module 2: Rfam RNA Motif / Family Mapping
# ==============================================================================

def map_pdb_to_rfam(
    ids_path: Path = DEFAULT_PFAM_IDS_FILE,
    rfam_path: Path = DEFAULT_RFAM_MAPPING_FILE,
    output_tsv: Optional[Path] = DEFAULT_RFAM_OUTPUT_TSV,
    output_pretty: Optional[Path] = DEFAULT_RFAM_PRETTY_OUTPUT,
    output_xlsx: Optional[Path] = DEFAULT_RFAM_OUTPUT_XLSX,
    output_json: Optional[Path] = DEFAULT_RFAM_OUTPUT_JSON
    ) -> List[List[str]]:
    """
    Streams the Rfam.pdb database and exports matching RNA family / motif records for target PDB IDs.
    """
    ids_path = Path(ids_path).resolve()
    rfam_path = Path(rfam_path).resolve()

    if not rfam_path.exists():
        raise FileNotFoundError(f"Rfam mapping file not found: {rfam_path}")

    ordered_ids, target_set = load_target_pdb_ids(ids_path)
    logger.info(f"[Rfam Mode] Loaded {len(target_set)} target PDB IDs from: {ids_path}")
    logger.info(f"[Rfam Mode] Scanning Rfam database: {rfam_path} ({rfam_path.stat().st_size / (1024*1024):.2f} MB)...")

    matched_records = []
    matched_pdbs = set()

    with open(rfam_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue

            cols = line_str.split("\t") if "\t" in line_str else line_str.split()
            if len(cols) >= 2:
                # Column 1 contains the PDB ID
                pdb_id = cols[1].strip().upper()
                if pdb_id in target_set:
                    row = [cols[i].strip() if i < len(cols) else "" for i in range(len(RFAM_COLUMNS))]
                    row[1] = pdb_id  # Uppercase PDB ID
                    matched_records.append(row)
                    matched_pdbs.add(pdb_id)

    logger.info(f"[Rfam Mode] Mapping complete: Found {len(matched_records)} RNA motif records for {len(matched_pdbs)} / {len(target_set)} PDB IDs.")

    write_formatted_outputs(
        records=matched_records,
        columns=RFAM_COLUMNS,
        output_tsv=output_tsv,
        output_pretty=output_pretty,
        output_xlsx=output_xlsx,
        output_json=output_json
    )

    return matched_records


# ==============================================================================
# CLI Entrypoint
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Unified PDB Domain & RNA Motif Mapping Pipeline (Pfam & Rfam).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--mode",
        choices=["pfam", "rfam", "all"],
        default="pfam",
        help="Select mapping database: 'pfam' (Protein domains), 'rfam' (RNA motifs), or 'all'."
    )
    parser.add_argument(
        "-i", "--input-ids",
        dest="ids_file",
        default=str(DEFAULT_PFAM_IDS_FILE),
        help="Path to file containing target PDB IDs (used for both Pfam and Rfam)."
    )
    parser.add_argument(
        "-m", "--mapping-file",
        dest="mapping_file",
        default=None,
        help="Path to mapping database file (pdb_pfam_mapping.txt or Rfam.pdb)."
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_tsv",
        default=None,
        help="Path to write output TSV file."
    )
    parser.add_argument(
        "-p", "--pretty-output",
        dest="pretty_output",
        default=None,
        help="Path to write column-aligned readable text file."
    )
    parser.add_argument(
        "-x", "--xlsx",
        dest="output_xlsx",
        default=None,
        help="Optional path to write Excel (.xlsx) file."
    )
    parser.add_argument(
        "-j", "--json",
        dest="output_json",
        default=None,
        help="Optional path to write structured JSON file."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    mode = args.mode
    ids_file = Path(args.ids_file) if args.ids_file else DEFAULT_PFAM_IDS_FILE

    try:
        if mode in ("pfam", "all"):
            map_file = Path(args.mapping_file) if args.mapping_file else DEFAULT_PFAM_MAPPING_FILE
            out_tsv = Path(args.output_tsv) if args.output_tsv else DEFAULT_PFAM_OUTPUT_TSV
            out_pretty = Path(args.pretty_output) if args.pretty_output else DEFAULT_PFAM_PRETTY_OUTPUT
            out_xlsx = Path(args.output_xlsx) if args.output_xlsx else DEFAULT_PFAM_OUTPUT_XLSX
            out_json = Path(args.output_json) if args.output_json else DEFAULT_PFAM_OUTPUT_JSON

            map_pdb_to_pfam(
                ids_path=ids_file,
                mapping_path=map_file,
                output_tsv=out_tsv,
                output_pretty=out_pretty,
                output_xlsx=out_xlsx,
                output_json=out_json
            )

        if mode in ("rfam", "all"):
            map_file = Path(args.mapping_file) if (args.mapping_file and mode == "rfam") else DEFAULT_RFAM_MAPPING_FILE
            out_tsv = Path(args.output_tsv) if (args.output_tsv and mode == "rfam") else DEFAULT_RFAM_OUTPUT_TSV
            out_pretty = Path(args.pretty_output) if (args.pretty_output and mode == "rfam") else DEFAULT_RFAM_PRETTY_OUTPUT
            out_xlsx = Path(args.output_xlsx) if (args.output_xlsx and mode == "rfam") else DEFAULT_RFAM_OUTPUT_XLSX
            out_json = Path(args.output_json) if (args.output_json and mode == "rfam") else DEFAULT_RFAM_OUTPUT_JSON

            map_pdb_to_rfam(
                ids_path=ids_file,
                rfam_path=map_file,
                output_tsv=out_tsv,
                output_pretty=out_pretty,
                output_xlsx=out_xlsx,
                output_json=out_json
            )

    except Exception as e:
        logger.error(f"Error during domain mapping: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
