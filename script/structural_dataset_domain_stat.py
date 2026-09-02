#!/usr/bin/env python3
"""
structural_dataset_domain_stat.py
=================================
Compiles and standardizes Pfam domain statistics and coordinate mappings from structural PDB-Pfam datasets.

Features:
  - Aggregates unique PDB counts and occurrence frequencies per Pfam domain.
  - Extracts standard and author domain start/end positions per PDB chain.
  - Direct in-memory processing without intermediate files.
  - Generates a multi-sheet Excel workbook ('Domain_Statistics', 'Domain_Positions', 'Full_Dataset')
    and structured JSON.
  - Supports CLI configuration with argparse.

Usage:
  # 1. Default execution using repository dataset:
  python3 structural_dataset_domain_stat.py

  # 2. Custom input mapping file and output destination:
  python3 structural_dataset_domain_stat.py -i input_data/domains/377_TF_domain_mapped.tsv -o results/domains/domain_stat.xlsx -j results/domains/domain_stat.json

  # 3. Also export domain summary statistics to CSV:
  python3 structural_dataset_domain_stat.py --csv results/domains/domain_stats.csv
"""

import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Tuple, Optional, Any
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

DEFAULT_INPUT_FILE = PROJECT_ROOT / "input_data" / "domains" / "377_TF_domain_mapped.tsv"
DEFAULT_OUTPUT_FILE = PROJECT_ROOT / "results" / "domains" / "domain_stat_structure_dataset_TF_NRD_final_377.xlsx"
DEFAULT_JSON_OUTPUT_FILE = PROJECT_ROOT / "results" / "domains" / "domain_stat_structure_dataset_TF_NRD_final_377.json"


def load_dataset(input_path: Path) -> pd.DataFrame:
    """Loads dataset supporting TSV, CSV, and Excel formats."""
    ext = input_path.suffix.lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        return pd.read_excel(input_path, dtype=str)
    elif ext == ".tsv":
        return pd.read_csv(input_path, sep="\t", dtype=str)
    elif ext == ".csv":
        return pd.read_csv(input_path, dtype=str)
    elif ext == ".txt":
        try:
            return pd.read_csv(input_path, sep="\t", dtype=str)
        except Exception:
            return pd.read_csv(input_path, dtype=str)
    else:
        raise ValueError(
            f"Unsupported file format '{ext}' for '{input_path.name}'. "
            "Please provide a TSV, CSV, or Excel spreadsheet."
        )


def process_structural_domain_stats(
    input_path: Path = DEFAULT_INPUT_FILE,
    output_path: Optional[Path] = DEFAULT_OUTPUT_FILE,
    json_output_path: Optional[Path] = DEFAULT_JSON_OUTPUT_FILE,
    csv_output_path: Optional[Path] = None
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parses PDB-Pfam mapped data and generates domain statistics and coordinate summaries.

    Args:
        input_path: Path to mapped PDB-Pfam input file (TSV/CSV/Excel).
        output_path: Path to save final multi-sheet Excel workbook.
        json_output_path: Path to save structured JSON output.
        csv_output_path: Optional path to export domain statistics CSV.

    Returns:
        Tuple of (domain_stats_df, domain_positions_df, full_df).
    """
    input_path = Path(input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    logger.info(f"Loading PDB-Pfam mapped dataset from: {input_path}")
    df = load_dataset(input_path)
    logger.info(f"Loaded {len(df)} mapped domain records.")

    # Clean whitespace in column names and values
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        df[col] = df[col].astype(str).str.strip().replace("nan", "")

    required_cols = ["PDB", "CHAIN", "PFAM_ACCESSION", "PFAM_NAME"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Required columns missing from input dataset: {missing}")

    # 1. Generate Domain Statistics Summary
    logger.info("Aggregating domain statistics...")
    stats_df = (
        df[df["PFAM_ACCESSION"] != ""]
        .groupby(["PFAM_ACCESSION", "PFAM_NAME"], as_index=False)
        .agg(
            Number_of_PDBs=("PDB", "nunique"),
            Total_Domain_Instances=("PDB", "count"),
            PDB_IDs=("PDB", lambda x: ", ".join(sorted(x[x != ""].unique())))
        )
        .sort_values(by=["Number_of_PDBs", "Total_Domain_Instances"], ascending=[False, False])
        .reset_index(drop=True)
    )

    logger.info(f"Identified {len(stats_df)} unique Pfam domains across {df['PDB'].nunique()} PDB entries.")
    logger.info("Top 5 most frequent Pfam domains in dataset:")
    for idx, row in stats_df.head(5).iterrows():
        logger.info(f"  {idx + 1}. {row['PFAM_ACCESSION']} ({row['PFAM_NAME']}) : {row['Number_of_PDBs']} PDBs ({row['Total_Domain_Instances']} instances)")

    # 2. Extract Standard Domain Positions
    position_cols = [
        "PDB", "CHAIN", "PFAM_ACCESSION", "PFAM_NAME",
        "PDB_START", "PDB_END", "AUTH_PDBRES_START", "AUTH_PDBRES_END",
        "UNIPROT_ACCESSION", "UNP_START", "UNP_END"
    ]
    available_pos_cols = [c for c in position_cols if c in df.columns]
    positions_df = df[available_pos_cols].copy()

    # Create explicit Start/End Position columns (fallback to AUTH if PDB is empty)
    if "PDB_START" in positions_df.columns:
        positions_df["Start_Position"] = positions_df["PDB_START"].where(
            positions_df["PDB_START"] != "",
            positions_df.get("AUTH_PDBRES_START", "")
        )
    if "PDB_END" in positions_df.columns:
        positions_df["End_Position"] = positions_df["PDB_END"].where(
            positions_df["PDB_END"] != "",
            positions_df.get("AUTH_PDBRES_END", "")
        )

    # 3. Save to Multi-Sheet Excel Workbook
    if output_path:
        output_path = Path(output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
            stats_df.to_excel(writer, sheet_name="Domain_Statistics", index=False)
            positions_df.to_excel(writer, sheet_name="Domain_Positions", index=False)
            df.to_excel(writer, sheet_name="Full_Dataset", index=False)
        logger.info(f"Saved multi-sheet Excel workbook -> {output_path}")

    # 4. Save structured JSON output
    if json_output_path:
        json_output_path = Path(json_output_path).resolve()
        json_output_path.parent.mkdir(parents=True, exist_ok=True)
        json_data = {
            "domain_statistics": stats_df.to_dict(orient="records"),
            "domain_positions": positions_df.to_dict(orient="records")
        }
        with open(json_output_path, "w", encoding="utf-8") as jf:
            json.dump(json_data, jf, indent=2, ensure_ascii=False)
        logger.info(f"Saved structured JSON output -> {json_output_path}")

    # 5. Export summary CSV if requested
    if csv_output_path:
        csv_output_path = Path(csv_output_path).resolve()
        csv_output_path.parent.mkdir(parents=True, exist_ok=True)
        stats_df.to_csv(csv_output_path, index=False)
        logger.info(f"Saved domain statistics CSV -> {csv_output_path}")

    return stats_df, positions_df, df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compile Pfam domain statistics and coordinate summaries from structural PDB datasets into Excel and JSON.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_file",
        default=str(DEFAULT_INPUT_FILE),
        help="Path to mapped PDB-Pfam dataset file (TSV/CSV/Excel)."
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        default=str(DEFAULT_OUTPUT_FILE),
        help="Path to save the final multi-sheet Excel workbook."
    )
    parser.add_argument(
        "-j", "--json",
        dest="json_output",
        default=str(DEFAULT_JSON_OUTPUT_FILE),
        help="Path to save structured JSON output."
    )
    parser.add_argument(
        "--csv",
        dest="csv_file",
        default=None,
        help="Optional path to export domain statistics summary table to CSV."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_file)
    output_path = Path(args.output_file) if args.output_file else None
    json_path = Path(args.json_output) if args.json_output else None
    csv_path = Path(args.csv_file) if args.csv_file else None

    try:
        process_structural_domain_stats(
            input_path=input_path,
            output_path=output_path,
            json_output_path=json_path,
            csv_output_path=csv_path
        )
    except Exception as e:
        logger.error(f"Error during structural domain statistics processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
