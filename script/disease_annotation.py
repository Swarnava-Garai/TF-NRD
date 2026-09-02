#!/usr/bin/env python3
"""
disease_annotation.py
=====================
Cleans and standardizes UniProt disease involvement annotations in-memory from an input Excel file.

Processes raw UniProt disease annotations without creating intermediate files:
  1. Filters entries with valid 'Involvement in disease' annotations.
  2. Aggregates metadata per unique UniProt 'Entry'.
  3. Cleans ID fields ('GeneID', 'KEGG') removing redundant delimiter formatting.
  4. Filters out entries containing only non-specific note records ('DISEASE: Note').
  5. Extracts and normalizes disease names via regex, removing citation tags and brackets.
  6. Exports the final cleaned dataset into a formatted Excel workbook ('ByRow' and 'AllNames' sheets)
     and structured JSON.

Usage:
  # 1. Default execution using repository dataset:
  python3 disease_annotation.py

  # 2. Custom input and output paths:
  python3 disease_annotation.py -i raw_data.xlsx -o results/disease/cleaned.xlsx -j results/disease/cleaned.json

  # 3. Also export plain text disease names list:
  python3 disease_annotation.py -t results/disease/disease_names.txt
"""

import sys
import re
import json
import logging
import argparse
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any
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

DEFAULT_INPUT_FILE = PROJECT_ROOT / "input_data" / "disease" / "TF_377_plus_nr_sequence_TF_involvement_in_disease_raw_data.xlsx"
DEFAULT_OUTPUT_FILE = PROJECT_ROOT / "results" / "disease" / "TF_377_plus_nr_disease_annotation.xlsx"
DEFAULT_JSON_OUTPUT_FILE = PROJECT_ROOT / "results" / "disease" / "TF_377_plus_nr_disease_annotation.json"

# Regex pattern for extracting disease names from UniProt text
DISEASE_NAME_PATTERN = re.compile(r"DISEASE:\s*([^:\[\.\n]+)", flags=re.IGNORECASE)


def extract_clean_disease_names(text: str) -> List[str]:
    """
    Extracts and cleans disease names from a raw UniProt disease annotation string.

    Args:
        text: Raw text string from 'Involvement in disease'.

    Returns:
        List of cleaned disease names.
    """
    if pd.isna(text):
        return []

    raw_matches = DISEASE_NAME_PATTERN.findall(str(text))
    cleaned_names = []

    for name in raw_matches:
        # Normalize whitespace and strip common punctuation
        name = re.sub(r"\s+", " ", name).strip(" ;,\t\r\n")
        # Strip surrounding quotes and bracket artifacts
        name = re.sub(r"^['\"]|['\"]$", "", name).strip()
        name = name.replace("[", "").replace("]", "").strip()

        # Skip non-disease notes (e.g. Note=...)
        if name and not name.lower().startswith("note="):
            cleaned_names.append(name)

    return cleaned_names


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


def process_disease_annotations(
    input_path: Path = DEFAULT_INPUT_FILE,
    output_path: Optional[Path] = DEFAULT_OUTPUT_FILE,
    json_output_path: Optional[Path] = DEFAULT_JSON_OUTPUT_FILE,
    text_output_path: Optional[Path] = None
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Executes the end-to-end in-memory disease annotation cleaning pipeline.

    Args:
        input_path: Path to the input Excel spreadsheet or TSV/CSV file.
        output_path: Path to save the final cleaned Excel file.
        json_output_path: Optional path to save structured JSON output.
        text_output_path: Optional path to export disease names as a plain list.

    Returns:
        Tuple of (by_row_df, all_names_exploded_df).
    """
    input_path = Path(input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    logger.info(f"Loading raw disease dataset from: {input_path}")
    raw_df = load_dataset(input_path)
    logger.info(f"Loaded {len(raw_df)} total records.")

    # 1. Filter rows containing disease involvement
    disease_col = "Involvement in disease"
    if disease_col not in raw_df.columns:
        raise ValueError(f"Required column '{disease_col}' not found in input spreadsheet.")

    df_disease = raw_df[raw_df[disease_col].notna()].copy()
    logger.info(f"Identified {len(df_disease)} entries with disease annotations.")

    # 2. Aggregate / group by unique UniProt Entry
    agg_dict = {
        "Entry Name": lambda x: ", ".join(x.dropna().astype(str)),
        "Protein names": lambda x: ", ".join(x.dropna().astype(str)),
        "Gene Names": lambda x: ", ".join(x.dropna().astype(str)),
        "GeneID": lambda x: ", ".join(x.dropna().astype(str)),
        "Organism": lambda x: ", ".join(x.dropna().astype(str)),
        disease_col: lambda x: ", ".join(x.dropna().astype(str)),
        "KEGG": lambda x: ", ".join(x.dropna().astype(str)),
        "PDB": lambda x: ", ".join(x.dropna().astype(str))
    }

    # Filter to aggregate only columns present in the input dataframe
    available_agg = {k: v for k, v in agg_dict.items() if k in df_disease.columns}
    merged_df = df_disease.groupby("Entry", as_index=False).agg(available_agg)

    # 3. Clean delimiters from identifier columns
    if "GeneID" in merged_df.columns:
        merged_df["GeneID"] = merged_df["GeneID"].astype(str).str.replace(";", "").str.strip()
    if "KEGG" in merged_df.columns:
        merged_df["KEGG"] = merged_df["KEGG"].astype(str).str.replace(";", "").str.strip()

    # 4. Filter out entries that only contain general disease notes ('DISEASE: Note')
    df_filtered = merged_df[
        ~merged_df[disease_col].astype(str).str.startswith("DISEASE: Note")
    ].copy().reset_index(drop=True)

    # 5. Extract and normalize disease names
    df_filtered["Disease_Names_List"] = df_filtered[disease_col].apply(extract_clean_disease_names)
    df_filtered["Disease_Names"] = df_filtered["Disease_Names_List"].apply(lambda lst: "; ".join(lst))

    # Reorder columns with Disease_Names at the end
    desired_cols = [
        "Entry", "Entry Name", "Protein names", "Gene Names", "GeneID",
        "Organism", disease_col, "KEGG", "PDB", "Disease_Names"
    ]
    final_cols = [c for c in desired_cols if c in df_filtered.columns]
    by_row_df = df_filtered[final_cols].copy()

    # 6. Explode into individual disease records (AllNames)
    meta_cols = [c for c in ["Entry", "Entry Name", "Gene Names", "Protein names"] if c in df_filtered.columns]
    exploded_df = (
        df_filtered[meta_cols + ["Disease_Names_List"]]
        .explode("Disease_Names_List")
        .dropna(subset=["Disease_Names_List"])
        .rename(columns={"Disease_Names_List": "Disease_Name"})
        .reset_index(drop=True)
    )

    logger.info(f"Cleaned {len(by_row_df)} protein entries associated with {len(exploded_df)} total disease annotations ({exploded_df['Disease_Name'].nunique()} unique diseases).")

    # 7. Write directly to final Excel workbook
    if output_path:
        output_path = Path(output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
            by_row_df.to_excel(writer, sheet_name="ByRow", index=False)
            exploded_df.to_excel(writer, sheet_name="AllNames", index=False)
        logger.info(f"Saved final multi-sheet Excel file -> {output_path}")

    # 8. Write structured JSON output
    if json_output_path:
        json_output_path = Path(json_output_path).resolve()
        json_output_path.parent.mkdir(parents=True, exist_ok=True)
        json_records = []
        for _, row in df_filtered.iterrows():
            rec = {
                "entry": row.get("Entry", ""),
                "entry_name": row.get("Entry Name", ""),
                "protein_names": row.get("Protein names", ""),
                "gene_names": row.get("Gene Names", ""),
                "gene_id": row.get("GeneID", ""),
                "organism": row.get("Organism", ""),
                "kegg": row.get("KEGG", ""),
                "pdb": row.get("PDB", ""),
                "disease_names": row.get("Disease_Names_List", []),
                "raw_disease_text": row.get(disease_col, "")
            }
            json_records.append(rec)

        with open(json_output_path, "w", encoding="utf-8") as jf:
            json.dump(json_records, jf, indent=2, ensure_ascii=False)
        logger.info(f"Saved structured JSON output -> {json_output_path}")

    # 9. Export plain text list if requested
    if text_output_path:
        text_output_path = Path(text_output_path).resolve()
        text_output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(text_output_path, "w", encoding="utf-8") as f:
            for name in exploded_df["Disease_Name"]:
                f.write(f"{name}\n")
        logger.info(f"Saved disease names list -> {text_output_path}")

    return by_row_df, exploded_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean and standardize UniProt disease involvement annotations directly into Excel and JSON.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_file",
        default=str(DEFAULT_INPUT_FILE),
        help="Path to input raw UniProt disease annotations Excel file."
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        default=str(DEFAULT_OUTPUT_FILE),
        help="Path to save the final cleaned Excel workbook."
    )
    parser.add_argument(
        "-j", "--json",
        dest="json_output",
        default=str(DEFAULT_JSON_OUTPUT_FILE),
        help="Path to save structured JSON output."
    )
    parser.add_argument(
        "-t", "--text-output",
        dest="text_output",
        default=None,
        help="Optional path to save cleaned disease names as plain text list."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_file)
    output_path = Path(args.output_file) if args.output_file else None
    json_path = Path(args.json_output) if args.json_output else None
    text_path = Path(args.text_output) if args.text_output else None

    try:
        process_disease_annotations(
            input_path=input_path,
            output_path=output_path,
            json_output_path=json_path,
            text_output_path=text_path
        )
    except Exception as e:
        logger.error(f"Error processing disease annotations: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
