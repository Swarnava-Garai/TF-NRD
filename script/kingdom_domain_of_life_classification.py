#!/usr/bin/env python3
"""
kingdom_domain_of_life_classification.py
========================================
Extracts and standardizes the Domain of Life (Eukaryota, Bacteria, Archaea, Virus)
from UniProt taxonomic lineage annotations without generating intermediate files.

Pipeline:
  1. Loads raw UniProt taxonomic lineage dataset from Excel.
  2. Extracts Domain of Life using regex pattern r'([A-Za-z0-9_\\- ]+)\\s*\\(domain\\)'.
  3. Assigns 'Virus' for records lacking cellular domain annotations (e.g. viral factors).
  4. Groups by UniProt accession ('From' / 'Entry') and aggregates metadata.
  5. Cleans and deduplicates semicolon/comma-separated strings in aggregated cells.
  6. Exports directly to the final cleaned Excel workbook and structured JSON file.

Usage:
  # 1. Default execution using repository dataset:
  python3 kingdom_domain_of_life_classification.py

  # 2. Custom input and output paths:
  python3 kingdom_domain_of_life_classification.py -i raw_data.xlsx -o results/domain_of_life/cleaned.xlsx -j results/domain_of_life/cleaned.json
"""

import sys
import re
import json
import logging
import argparse
from pathlib import Path
from typing import Optional, List, Dict, Any
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

DEFAULT_INPUT_FILE = PROJECT_ROOT / "input_data" / "domain_of_life" / "Final_domain_of_life_structure_dataset_rawdata_377.xlsx"
DEFAULT_OUTPUT_FILE = PROJECT_ROOT / "results" / "domain_of_life" / "classify_TFs_domain_organism_cleaned.xlsx"
DEFAULT_JSON_OUTPUT_FILE = PROJECT_ROOT / "results" / "domain_of_life" / "classify_TFs_domain_organism_cleaned.json"

# Regex pattern for extracting domain from taxonomic lineage string
DOMAIN_PATTERN = re.compile(r"([A-Za-z0-9_\- ]+)\s*\(domain\)", flags=re.IGNORECASE)


def extract_domain_from_lineage(lineage_text: Optional[str]) -> str:
    """
    Extracts domain name from UniProt taxonomic lineage text.
    Returns 'Virus' if no cellular domain is found.
    """
    if pd.isna(lineage_text) or not str(lineage_text).strip():
        return "Virus"

    match = DOMAIN_PATTERN.search(str(lineage_text))
    if match:
        domain = match.group(1).strip()
        return domain if domain else "Virus"
    return "Virus"


def deduplicate_cell_values(val, default_empty: str = "") -> str:
    """
    Deduplicates items in a comma- or semicolon-separated cell value
    while preserving original encounter order.
    """
    if pd.isna(val) or str(val).strip() == "":
        return default_empty

    parts = [p.strip() for p in str(val).replace(";", ",").split(",") if p.strip()]
    seen = set()
    unique = []
    for p in parts:
        if p not in seen:
            seen.add(p)
            unique.append(p)

    return ", ".join(unique) if unique else default_empty


def load_dataset(input_path: Path) -> pd.DataFrame:
    """
    Loads dataset supporting Excel (.xlsx, .xls) and Tabular text (.tsv, .csv) formats.
    """
    ext = input_path.suffix.lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        return pd.read_excel(input_path)
    elif ext == ".tsv":
        return pd.read_csv(input_path, sep="\t")
    elif ext == ".csv":
        return pd.read_csv(input_path)
    elif ext == ".txt":
        try:
            return pd.read_csv(input_path, sep="\t")
        except Exception:
            return pd.read_csv(input_path)
    else:
        raise ValueError(
            f"Unsupported file format '{ext}' for '{input_path.name}'. "
            "Please provide an Excel spreadsheet (.xlsx, .xls) or tabular dataset (.tsv, .csv)."
        )


def classify_domains_of_life(
    input_path: Path = DEFAULT_INPUT_FILE,
    output_path: Optional[Path] = DEFAULT_OUTPUT_FILE,
    json_output_path: Optional[Path] = DEFAULT_JSON_OUTPUT_FILE,
    group_by_col: str = "From"
) -> pd.DataFrame:
    """
    Executes domain of life extraction and aggregation in-memory.

    Args:
        input_path: Path to raw input Excel or TSV/CSV file.
        output_path: Path to save the final cleaned Excel file.
        json_output_path: Path to save the structured JSON file.
        group_by_col: Column to group by (default 'From').

    Returns:
        Cleaned pd.DataFrame.
    """
    input_path = Path(input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    logger.info(f"Loading raw taxonomic dataset from: {input_path}")
    df = load_dataset(input_path)
    logger.info(f"Loaded {len(df)} total rows.")

    # 1. Detect taxonomic lineage column
    tax_col = None
    for col in df.columns:
        if "taxonomic" in str(col).lower():
            tax_col = col
            break

    if tax_col is None:
        raise ValueError("Could not find a column containing 'Taxonomic lineage' in input spreadsheet.")

    # 2. Extract Domain of Life
    logger.info(f"Extracting Domain of Life from column '{tax_col}'...")
    df["Domain_Name"] = df[tax_col].apply(extract_domain_from_lineage)

    # 3. Aggregate by target group column (e.g. 'From')
    if group_by_col not in df.columns:
        if "Entry" in df.columns:
            group_by_col = "Entry"
            logger.warning(f"Target group column not found, falling back to '{group_by_col}'.")
        else:
            raise ValueError(f"Grouping column '{group_by_col}' not found in input spreadsheet.")

    logger.info(f"Grouping and aggregating records by '{group_by_col}'...")

    agg_targets = [
        "Entry", "Protein names", "Organism", "Gene Names",
        tax_col, "Organism (ID)", "Domain_Name"
    ]
    available_targets = [c for c in agg_targets if c in df.columns]

    agg_dict = {
        col: (lambda x, c=col: deduplicate_cell_values(", ".join(x.dropna().astype(str)), default_empty="Virus" if c == "Domain_Name" else ""))
        for col in available_targets
    }

    merged_df = df.groupby(group_by_col, as_index=False).agg(agg_dict)

    logger.info(f"Aggregated into {len(merged_df)} unique entries.")

    # Domain summary
    domain_counts = merged_df["Domain_Name"].value_counts().to_dict()
    logger.info("Domain distribution summary:")
    for dom, count in domain_counts.items():
        logger.info(f"  - {dom}: {count}")

    # 4. Save directly to output Excel file
    if output_path:
        output_path = Path(output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged_df.to_excel(output_path, index=False)
        logger.info(f"Saved final cleaned dataset -> {output_path}")

    # 5. Save structured JSON output
    if json_output_path:
        json_output_path = Path(json_output_path).resolve()
        json_output_path.parent.mkdir(parents=True, exist_ok=True)
        records = merged_df.to_dict(orient="records")
        with open(json_output_path, "w", encoding="utf-8") as jf:
            json.dump(records, jf, indent=2, ensure_ascii=False)
        logger.info(f"Saved structured JSON output -> {json_output_path}")

    return merged_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify and standardize Domain of Life from UniProt taxonomic lineage into Excel and JSON.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_file",
        default=str(DEFAULT_INPUT_FILE),
        help="Path to raw UniProt dataset Excel file."
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        default=str(DEFAULT_OUTPUT_FILE),
        help="Path to save the final cleaned Excel file."
    )
    parser.add_argument(
        "-j", "--json",
        dest="json_output",
        default=str(DEFAULT_JSON_OUTPUT_FILE),
        help="Path to save structured JSON output."
    )
    parser.add_argument(
        "-g", "--group-by",
        dest="group_by",
        default="From",
        help="Column name to group and aggregate by."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_file)
    output_path = Path(args.output_file) if args.output_file else None
    json_path = Path(args.json_output) if args.json_output else None

    try:
        classify_domains_of_life(
            input_path=input_path,
            output_path=output_path,
            json_output_path=json_path,
            group_by_col=args.group_by
        )
    except Exception as e:
        logger.error(f"Error during domain classification: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
