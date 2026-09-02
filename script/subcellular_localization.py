#!/usr/bin/env python3
"""
subcellular_localization.py
===========================
Cleans, standardizes, and classifies UniProt subcellular localization annotations for transcription factor complexes.

Pipeline:
  1. Parses raw UniProt 'Subcellular location [CC]' annotations.
  2. Strips evidence tags ({ECO:...}), isoform brackets ([...]), and freeform note sections (Note=...).
  3. Groups by query identifier ('From' / PDB ID) and deduplicates compartment terms.
  4. Classifies each complex into compartmental categories:
     - 'Only_Nucleus': Exclusively nuclear.
     - 'Only_Cytoplasm': Exclusively cytoplasmic.
     - 'Nucleus_Cytoplasm_Both': Shuttles between or localized to both Nucleus and Cytoplasm.
     - 'Nucleus_with_Other_not_Cytoplasm': Nuclear plus other organelle(s) excluding cytoplasm.
     - 'Cytoplasm_with_Other_not_Nucleus': Cytoplasmic plus other organelle(s) excluding nucleus.
     - 'Others': Membrane, extracellular, or organelle-only without nucleus/cytoplasm.
  5. Exports a structured multi-sheet Excel workbook ('Summary', 'Category_PDBs', 'Detailed_Per_PDB')
     and structured JSON.

Usage:
  # 1. Default execution using repository dataset:
  python3 subcellular_localization.py

  # 2. Custom input and output paths:
  python3 subcellular_localization.py -i input_data.xlsx -o results/subcellular_location/detailed.xlsx -j results/subcellular_location/detailed.json

  # 3. Also export summary statistics to CSV:
  python3 subcellular_localization.py --csv results/subcellular_location/summary.csv
"""

import sys
import re
import json
import logging
import argparse
from pathlib import Path
from typing import Tuple, List, Dict, Optional, Any
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

DEFAULT_INPUT_FILE = PROJECT_ROOT / "input_data" / "subcellular_location" / "Final_TF_377_to_uniprot.xlsx"
DEFAULT_OUTPUT_FILE = PROJECT_ROOT / "results" / "subcellular_location" / "Structure_subcellular_location_detailed.xlsx"
DEFAULT_JSON_OUTPUT_FILE = PROJECT_ROOT / "results" / "subcellular_location" / "Structure_subcellular_location_detailed.json"


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


def clean_subcellular_text(text: Optional[str]) -> str:
    """
    Strips evidence tags, isoform brackets, and note sections from UniProt location text.
    Deduplicates compartment tokens while preserving order.
    """
    if pd.isna(text) or not str(text).strip():
        return ""

    s = str(text).replace("SUBCELLULAR LOCATION:", "")
    # Remove everything starting at Note=
    s = re.split(r"Note=", s, flags=re.IGNORECASE)[0]
    # Remove isoform/variant brackets [ ... ]
    s = re.sub(r"\[.*?\]", "", s)
    # Remove evidence tags {ECO:...}
    s = re.sub(r"\{ECO:.*?\}", "", s)

    # Tokenize by common delimiters
    parts = re.split(r"[.;,/|\n]+", s)
    parts = [p.strip() for p in parts if p.strip()]

    # Deduplicate case-insensitively
    seen = set()
    unique_tokens = []
    for p in parts:
        lower_p = p.lower()
        if lower_p not in seen:
            seen.add(lower_p)
            unique_tokens.append(p)

    return ", ".join(unique_tokens)


def classify_subcellular_category(cleaned_location_str: str) -> str:
    """
    Categorizes localization into nuclear/cytoplasmic/other buckets.
    """
    if not cleaned_location_str:
        return "Others"

    tokens = [t.strip().lower() for t in re.split(r"[.,;]", cleaned_location_str) if t.strip()]
    has_nucleus = any("nucleus" in t for t in tokens)
    has_cytoplasm = any("cytoplasm" in t for t in tokens)

    if has_nucleus and not has_cytoplasm:
        return "Only_Nucleus" if len(tokens) == 1 else "Nucleus_with_Other_not_Cytoplasm"
    elif has_cytoplasm and not has_nucleus:
        return "Only_Cytoplasm" if len(tokens) == 1 else "Cytoplasm_with_Other_not_Nucleus"
    elif has_nucleus and has_cytoplasm:
        return "Nucleus_Cytoplasm_Both"
    else:
        return "Others"


def process_subcellular_localization(
    input_path: Path = DEFAULT_INPUT_FILE,
    output_path: Optional[Path] = DEFAULT_OUTPUT_FILE,
    json_output_path: Optional[Path] = DEFAULT_JSON_OUTPUT_FILE,
    csv_output_path: Optional[Path] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Executes the full subcellular localization parsing and categorization pipeline.

    Args:
        input_path: Path to input dataset file.
        output_path: Path to save final Excel workbook.
        json_output_path: Path to save structured JSON output.
        csv_output_path: Optional path to export summary CSV.

    Returns:
        Tuple of (detailed_df, category_columns_df, summary_df).
    """
    input_path = Path(input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    logger.info(f"Loading UniProt dataset from: {input_path}")
    df = load_dataset(input_path)
    logger.info(f"Loaded {len(df)} records.")

    # 1. Detect target columns
    from_col = "From" if "From" in df.columns else "PDB"
    if from_col not in df.columns:
        from_col = df.columns[0]

    subcell_col = next((c for c in df.columns if "subcellular" in str(c).lower()), None)
    if subcell_col is None:
        raise ValueError("Could not find a column containing 'Subcellular location' in input dataset.")

    # Filter rows with non-null ID and location
    df_valid = df[df[from_col].notna() & df[subcell_col].notna()].copy()
    df_valid["Cleaned_Location"] = df_valid[subcell_col].apply(clean_subcellular_text)
    df_valid = df_valid[df_valid["Cleaned_Location"] != ""].copy()

    logger.info(f"Identified {len(df_valid)} valid subcellular localization annotations.")

    # 2. Group by query ID (From / PDB)
    agg_dict = {
        subcell_col: lambda x: clean_subcellular_text(", ".join(x.astype(str)))
    }
    for meta_col in ["Entry", "Organism", "Protein names", "Gene Names"]:
        if meta_col in df_valid.columns:
            agg_dict[meta_col] = lambda x: ", ".join(sorted(x.dropna().astype(str).unique()))

    detailed_df = df_valid.groupby(from_col, as_index=False).agg(agg_dict)
    detailed_df = detailed_df.rename(columns={subcell_col: "Subcellular_Location_Cleaned"})

    # 3. Categorize
    detailed_df["Category"] = detailed_df["Subcellular_Location_Cleaned"].apply(classify_subcellular_category)
    detailed_df["Component_Count"] = detailed_df["Subcellular_Location_Cleaned"].apply(
        lambda s: len([p for p in s.split(",") if p.strip()])
    )

    # 4. Generate Category Columns DataFrame (PDBs under each category)
    categories = [
        "Only_Nucleus",
        "Only_Cytoplasm",
        "Nucleus_Cytoplasm_Both",
        "Nucleus_with_Other_not_Cytoplasm",
        "Cytoplasm_with_Other_not_Nucleus",
        "Others"
    ]
    cat_dict = {}
    for cat in categories:
        cat_pdbs = detailed_df[detailed_df["Category"] == cat][from_col].tolist()
        cat_dict[cat] = pd.Series(cat_pdbs)

    category_columns_df = pd.DataFrame(cat_dict)

    # 5. Generate Summary Table
    total_entries = len(detailed_df)
    summary_records = []
    for cat in categories:
        count = (detailed_df["Category"] == cat).sum()
        pct = (count / total_entries * 100) if total_entries > 0 else 0.0
        summary_records.append({
            "Category": cat,
            "Count": count,
            "Percentage": round(pct, 2)
        })

    summary_df = pd.DataFrame(summary_records)

    logger.info(f"Categorized {total_entries} unique PDB / UniProt complexes:")
    for _, r in summary_df.iterrows():
        logger.info(f"  - {r['Category']} : {r['Count']} ({r['Percentage']}%)")

    # 6. Save to Multi-Sheet Excel Workbook
    if output_path:
        output_path = Path(output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
            summary_df.to_excel(writer, sheet_name="Summary", index=False)
            category_columns_df.to_excel(writer, sheet_name="Category_PDBs", index=False)
            detailed_df.to_excel(writer, sheet_name="Detailed_Per_PDB", index=False)
        logger.info(f"Saved multi-sheet Excel workbook -> {output_path}")

    # 7. Save structured JSON output
    if json_output_path:
        json_output_path = Path(json_output_path).resolve()
        json_output_path.parent.mkdir(parents=True, exist_ok=True)
        json_data = {
            "summary": summary_df.to_dict(orient="records"),
            "category_pdbs": {cat: detailed_df[detailed_df["Category"] == cat][from_col].tolist() for cat in categories},
            "detailed_per_entry": detailed_df.to_dict(orient="records")
        }
        with open(json_output_path, "w", encoding="utf-8") as jf:
            json.dump(json_data, jf, indent=2, ensure_ascii=False)
        logger.info(f"Saved structured JSON output -> {json_output_path}")

    # 8. Export summary CSV if requested
    if csv_output_path:
        csv_output_path = Path(csv_output_path).resolve()
        csv_output_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(csv_output_path, index=False)
        logger.info(f"Saved summary statistics CSV -> {csv_output_path}")

    return detailed_df, category_columns_df, summary_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean and classify UniProt subcellular localization annotations into Excel and JSON.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_file",
        default=str(DEFAULT_INPUT_FILE),
        help="Path to input dataset file (Excel/TSV/CSV)."
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
        help="Optional path to export summary statistics table to CSV."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_file)
    output_path = Path(args.output_file) if args.output_file else None
    json_path = Path(args.json_output) if args.json_output else None
    csv_path = Path(args.csv_file) if args.csv_file else None

    try:
        process_subcellular_localization(
            input_path=input_path,
            output_path=output_path,
            json_output_path=json_path,
            csv_output_path=csv_path
        )
    except Exception as e:
        logger.error(f"Error during subcellular localization processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
