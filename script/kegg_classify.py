#!/usr/bin/env python3
"""
kegg_classify.py
================
Classifies KEGG pathway and disease annotations into hierarchical functional and disease classes using the KEGG BRITE hierarchy.

Features:
  - Fetches and parses the official KEGG BRITE functional hierarchy (br:br08901) for all biological and human disease pathways.
  - Supports KEGG DISEASE entries (e.g., H00224) via KEGG REST flat-file CATEGORY parsing.
  - Automatically classifies pathway IDs (e.g., hsa05203, hsa04613) into Major Class (e.g., Human Diseases, Organismal Systems)
    and Subclass (e.g., Cancer: overview, Signal transduction, Infectious disease: viral).
  - Enriches input datasets directly in-memory without creating intermediate files.
  - Exports a structured multi-sheet Excel workbook ('Classified_Pathways', 'Major_Class_Summary', 'Subclass_Summary'),
    structured JSON, and CSV.
  - Logs execution details to script-named log file ('kegg_classify.log') and console.

Usage:
  # 1. Default execution using repository dataset:
  python3 kegg_classify.py

  # 2. Custom input and output paths:
  python3 kegg_classify.py -i input_data/kegg/TF_NRD_KEGG_pathways_combined_final.xlsx -o results/kegg/classified.xlsx -j results/kegg/classified.json

  # 3. Custom ID column specification:
  python3 kegg_classify.py -i custom_pathways.csv --col KEGG_ID
"""

import sys
import re
import json
import time
import logging
import argparse
import urllib.request
import urllib.error
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, Optional, List, Any
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

KEGG_BASE = "https://rest.kegg.jp"
DEFAULT_INPUT_FILE = PROJECT_ROOT / "input_data" / "kegg" / "TF_NRD_KEGG_pathways_combined_final.xlsx"
DEFAULT_OUTPUT_XLSX = PROJECT_ROOT / "results" / "kegg" / "TF_NRD_KEGG_pathways_classified.xlsx"
DEFAULT_OUTPUT_JSON = PROJECT_ROOT / "results" / "kegg" / "TF_NRD_KEGG_pathways_classified.json"
DEFAULT_OUTPUT_CSV = PROJECT_ROOT / "results" / "kegg" / "TF_NRD_KEGG_pathways_classified.csv"


# ==============================================================================
# Helper Functions & Dataset Loaders
# ==============================================================================

def load_dataset(input_path: Path) -> pd.DataFrame:
    """Loads dataset supporting Excel (.xlsx, .xls) and Tabular text (.tsv, .csv) formats."""
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


def http_get(url: str, max_retries: int = 4, backoff: float = 0.8) -> str:
    """Fetches text from URL using urllib with retry backoff."""
    req = urllib.request.Request(url, headers={"User-Agent": "TF-NRD-KEGG-Classifier/1.0"})
    for i in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=20) as response:
                if response.status == 200:
                    return response.read().decode("utf-8", errors="replace")
        except Exception as e:
            if i == max_retries - 1:
                raise RuntimeError(f"Failed to GET after {max_retries} retries from {url}: {e}")
            time.sleep(backoff * (2 ** i))
    raise RuntimeError(f"Failed to GET from {url}")


# ==============================================================================
# BRITE Functional Hierarchy Parsing
# ==============================================================================

def fetch_and_parse_brite_hierarchy() -> Dict[str, Tuple[str, str, str]]:
    """
    Fetches the KEGG BRITE functional hierarchy (br:br08901) and returns a mapping:
    key -> (Major_Class, Subclass, Pathway_Name)
    Keys include 5-digit number (e.g. '05203'), 'hsa05203', and 'map05203'.
    """
    logger.info("Fetching KEGG BRITE pathway functional hierarchy (br:br08901)...")
    url = f"{KEGG_BASE}/get/br:br08901"
    brite_text = http_get(url)

    mapping: Dict[str, Tuple[str, str, str]] = {}
    major = "Unknown"
    sub = "Unknown"

    for line in brite_text.splitlines():
        line_str = line.strip()
        if not line_str:
            continue
        if line.startswith("A") and not line.startswith("A<b>"):
            major = line[1:].strip()
            sub = "Unknown"
        elif line.startswith("B") and not line.startswith("B<b>"):
            sub = line[1:].strip()
        elif line.startswith("C"):
            parts = line[1:].strip().split(None, 1)
            if parts:
                num = parts[0].strip()
                pname = parts[1].strip() if len(parts) > 1 else ""
                mapping[num] = (major, sub, pname)
                mapping[f"map{num}"] = (major, sub, pname)
                mapping[f"hsa{num}"] = (major, sub, pname)

    logger.info(f"Loaded {len(mapping)} pathway BRITE hierarchy index entries.")
    return mapping


def parse_kegg_flat(text: str) -> Dict[str, List[str]]:
    """Parses a KEGG flat file into field: [values] mapping."""
    data = defaultdict(list)
    current_key = None
    for raw in text.splitlines():
        if not raw.strip():
            continue
        if len(raw) >= 12 and raw[:12].strip().isalpha() and raw[12] == " ":
            current_key = raw[:12].strip()
            value = raw[12:].strip()
            if value:
                data[current_key].append(value)
        else:
            if current_key is not None:
                data[current_key].append(raw.strip())
    return data


def get_disease_category_online(h_id: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Queries KEGG REST API for disease entry Hxxxxx and extracts name and category."""
    url = f"{KEGG_BASE}/get/{h_id}"
    try:
        txt = http_get(url)
    except Exception:
        return None, None, None

    data = parse_kegg_flat(txt)
    name = None
    if "NAME" in data and data["NAME"]:
        name = re.sub(r"\s*;$", "", data["NAME"][0]).strip()

    major = None
    subclass = None
    if "CATEGORY" in data and data["CATEGORY"]:
        first_cat = data["CATEGORY"][0]
        parts = [p.strip() for p in first_cat.split(";")]
        if parts:
            major = parts[0] or None
        if len(parts) > 1:
            subclass = parts[1] or None
    return name, major, subclass


# ==============================================================================
# Pipeline Processing Engine
# ==============================================================================

def process_kegg_classification(
    input_path: Path = DEFAULT_INPUT_FILE,
    output_xlsx: Optional[Path] = DEFAULT_OUTPUT_XLSX,
    output_json: Optional[Path] = DEFAULT_OUTPUT_JSON,
    output_csv: Optional[Path] = DEFAULT_OUTPUT_CSV,
    id_col: Optional[str] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Executes the KEGG pathway / disease classification pipeline.

    Args:
        input_path: Path to input dataset file.
        output_xlsx: Optional path to save multi-sheet Excel workbook.
        output_json: Optional path to save structured JSON.
        output_csv: Optional path to save CSV output.
        id_col: Optional column name for KEGG IDs.

    Returns:
        Tuple of (classified_df, major_summary_df, subclass_summary_df).
    """
    input_path = Path(input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    logger.info(f"Loading KEGG pathway dataset from: {input_path}")
    df = load_dataset(input_path)
    logger.info(f"Loaded {len(df)} records.")

    # 1. Detect KEGG ID column
    if id_col and id_col in df.columns:
        target_col = id_col
    else:
        candidates = ["KEGG_ID", "KEGG ID", "kegg_id", "KEGG", "ID", "hsa_id", "Entry", "From"]
        target_col = next((c for c in candidates if c in df.columns), df.columns[0])

    logger.info(f"Using column '{target_col}' as KEGG identifier.")

    # 2. Fetch BRITE Hierarchy
    try:
        brite_map = fetch_and_parse_brite_hierarchy()
    except Exception as e:
        logger.warning(f"Could not fetch live BRITE hierarchy from KEGG: {e}. Proceeding with fallback parsing.")
        brite_map = {}

    # 3. Classify Each Entry
    logger.info("Classifying pathways into Major Classes and Subclasses...")
    major_classes = []
    subclasses = []
    category_sources = []

    for idx, raw_id in enumerate(df[target_col]):
        kegg_str = str(raw_id).strip()
        # Clean ID
        kegg_clean = re.sub(r"[;,\s].*", "", kegg_str)

        major = "Unknown"
        sub = "Unknown"
        source = "Unclassified"

        # Check BRITE hierarchy first
        if kegg_clean in brite_map:
            major, sub, _ = brite_map[kegg_clean]
            source = "KEGG_BRITE_Pathway"
        elif kegg_clean.lower() in brite_map:
            major, sub, _ = brite_map[kegg_clean.lower()]
            source = "KEGG_BRITE_Pathway"
        elif re.sub(r"^[a-zA-Z]+", "", kegg_clean) in brite_map:
            num_only = re.sub(r"^[a-zA-Z]+", "", kegg_clean)
            major, sub, _ = brite_map[num_only]
            source = "KEGG_BRITE_Pathway"
        elif re.fullmatch(r"H\d{5}", kegg_clean, re.IGNORECASE):
            # Query KEGG DISEASE entry
            d_name, d_major, d_sub = get_disease_category_online(kegg_clean.upper())
            major = d_major or "Human Diseases"
            sub = d_sub or "Unspecified Disease"
            source = "KEGG_DISEASE_Category"

        major_classes.append(major)
        subclasses.append(sub)
        category_sources.append(source)

    df_classified = df.copy()
    df_classified["Major_Class"] = major_classes
    df_classified["Subclass"] = subclasses
    df_classified["Category_Source"] = category_sources

    # 4. Generate Summaries
    hits_col = next((c for c in ["Number_of_Hits", "Number of hit", "Hits", "Count"] if c in df_classified.columns), None)

    # Major Class Summary
    if hits_col:
        major_summary = (
            df_classified.groupby("Major_Class", as_index=False)
            .agg(
                Pathway_Count=(target_col, "count"),
                Total_Hits=(hits_col, "sum")
            )
            .sort_values(by=["Total_Hits", "Pathway_Count"], ascending=[False, False])
            .reset_index(drop=True)
        )
    else:
        major_summary = (
            df_classified.groupby("Major_Class", as_index=False)
            .agg(Pathway_Count=(target_col, "count"))
            .sort_values(by="Pathway_Count", ascending=False)
            .reset_index(drop=True)
        )

    tot_pathways = len(df_classified)
    major_summary["Pathway_Percentage"] = (major_summary["Pathway_Count"] / tot_pathways * 100).round(2)

    # Subclass Summary
    if hits_col:
        subclass_summary = (
            df_classified.groupby(["Major_Class", "Subclass"], as_index=False)
            .agg(
                Pathway_Count=(target_col, "count"),
                Total_Hits=(hits_col, "sum")
            )
            .sort_values(by=["Total_Hits", "Pathway_Count"], ascending=[False, False])
            .reset_index(drop=True)
        )
    else:
        subclass_summary = (
            df_classified.groupby(["Major_Class", "Subclass"], as_index=False)
            .agg(Pathway_Count=(target_col, "count"))
            .sort_values(by="Pathway_Count", ascending=False)
            .reset_index(drop=True)
        )

    subclass_summary["Pathway_Percentage"] = (subclass_summary["Pathway_Count"] / tot_pathways * 100).round(2)

    logger.info("Major Class Distribution:")
    for _, r in major_summary.iterrows():
        hits_str = f", {r['Total_Hits']} total hits" if hits_col else ""
        logger.info(f"  - {r['Major_Class']} : {r['Pathway_Count']} pathways ({r['Pathway_Percentage']}%{hits_str})")

    # 5. Export Multi-Sheet Excel Workbook
    if output_xlsx:
        output_xlsx = Path(output_xlsx).resolve()
        output_xlsx.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
            df_classified.to_excel(writer, sheet_name="Classified_Pathways", index=False)
            major_summary.to_excel(writer, sheet_name="Major_Class_Summary", index=False)
            subclass_summary.to_excel(writer, sheet_name="Subclass_Summary", index=False)
        logger.info(f"Saved final multi-sheet Excel workbook -> {output_xlsx}")

    # 6. Export Structured JSON
    if output_json:
        output_json = Path(output_json).resolve()
        output_json.parent.mkdir(parents=True, exist_ok=True)
        json_data = {
            "major_class_summary": major_summary.to_dict(orient="records"),
            "subclass_summary": subclass_summary.to_dict(orient="records"),
            "classified_pathways": df_classified.to_dict(orient="records")
        }
        with open(output_json, "w", encoding="utf-8") as jf:
            json.dump(json_data, jf, indent=2, ensure_ascii=False)
        logger.info(f"Saved structured JSON output -> {output_json}")

    # 7. Export CSV
    if output_csv:
        output_csv = Path(output_csv).resolve()
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        df_classified.to_csv(output_csv, index=False)
        logger.info(f"Saved classified CSV output -> {output_csv}")

    return df_classified, major_summary, subclass_summary


# ==============================================================================
# CLI Entrypoint
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify KEGG disease & biological pathway annotations into hierarchical classes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--in", "--input",
        dest="input_file",
        default=str(DEFAULT_INPUT_FILE),
        help="Path to input Excel or TSV/CSV dataset containing KEGG pathway IDs."
    )
    parser.add_argument(
        "-o", "--out", "--output",
        dest="output_xlsx",
        default=str(DEFAULT_OUTPUT_XLSX),
        help="Path to save the final multi-sheet Excel workbook."
    )
    parser.add_argument(
        "-j", "--json",
        dest="output_json",
        default=str(DEFAULT_OUTPUT_JSON),
        help="Path to save structured JSON output."
    )
    parser.add_argument(
        "--csv",
        dest="output_csv",
        default=str(DEFAULT_OUTPUT_CSV),
        help="Path to save classified CSV output."
    )
    parser.add_argument(
        "--col",
        dest="id_col",
        default=None,
        help="Column name containing KEGG IDs (auto-detected if omitted)."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_file)
    output_xlsx = Path(args.output_xlsx) if args.output_xlsx else None
    output_json = Path(args.output_json) if args.output_json else None
    output_csv = Path(args.output_csv) if args.output_csv else None

    try:
        process_kegg_classification(
            input_path=input_path,
            output_xlsx=output_xlsx,
            output_json=output_json,
            output_csv=output_csv,
            id_col=args.id_col
        )
    except Exception as e:
        logger.error(f"Error during KEGG classification: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
