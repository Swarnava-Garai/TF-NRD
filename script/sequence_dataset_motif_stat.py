#!/usr/bin/env python3
"""
sequence_dataset_motif_stat.py
==============================
Parses, cleans, and compiles statistical summaries of sequence motifs from UniProt annotations.

Features:
  - Extracts all MOTIF occurrences (ranges and /note descriptions) per entry, including multi-motif rows.
  - Generates comprehensive statistics: unique entry counts, total occurrence counts, and associated UniProt IDs.
  - Direct in-memory processing without creating intermediate files.
  - Exports a structured multi-sheet Excel workbook ('Motif_Statistics', 'Individual_Motifs', 'Annotated_Dataset')
    and structured JSON.

Usage:
  # 1. Default execution using repository dataset:
  python3 sequence_dataset_motif_stat.py

  # 2. Custom input and output paths:
  python3 sequence_dataset_motif_stat.py -i input_motifs.xlsx -o results/motif/motif_stat.xlsx -j results/motif/motif_stat.json

  # 3. Also export motif statistics to CSV:
  python3 sequence_dataset_motif_stat.py --csv results/motif/motif_stats.csv
"""

import sys
import re
import json
import logging
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
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

DEFAULT_INPUT_FILE = PROJECT_ROOT / "input_data" / "motif" / "nr_sequence_dataset_motif_details.xlsx"
DEFAULT_OUTPUT_FILE = PROJECT_ROOT / "results" / "motif" / "nr_sequence_dataset_motif_stat.xlsx"
DEFAULT_JSON_OUTPUT_FILE = PROJECT_ROOT / "results" / "motif" / "nr_sequence_dataset_motif_stat.json"


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


def extract_all_motifs(text: Optional[str]) -> List[Dict[str, str]]:
    """
    Extracts all motif instances from UniProt Motif text, handling single or multiple MOTIF blocks.

    Returns:
        List of dicts: [{'Motif_Range': '...', 'Motif_Name': '...'}, ...]
    """
    if pd.isna(text) or not str(text).strip():
        return []

    # Split into individual MOTIF blocks if multiple are present
    blocks = re.split(r"(?=MOTIF\s+)", str(text))
    extracted = []

    for b in blocks:
        b_str = b.strip()
        if not b_str.startswith("MOTIF"):
            continue

        r_match = re.search(r"MOTIF\s+([\d.]+)", b_str)
        n_match = re.search(r'/note="([^"]+)"', b_str)

        m_range = r_match.group(1) if r_match else ""
        m_name = n_match.group(1).strip() if n_match else ""

        if m_name or m_range:
            extracted.append({
                "Motif_Range": m_range,
                "Motif_Name": m_name
            })

    return extracted


def process_motif_statistics(
    input_path: Path = DEFAULT_INPUT_FILE,
    output_path: Optional[Path] = DEFAULT_OUTPUT_FILE,
    json_output_path: Optional[Path] = DEFAULT_JSON_OUTPUT_FILE,
    csv_output_path: Optional[Path] = None
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parses motif annotations and generates comprehensive statistical summaries.

    Args:
        input_path: Path to input dataset file.
        output_path: Path to save final Excel workbook.
        json_output_path: Path to save structured JSON output.
        csv_output_path: Optional path to export statistics summary CSV.

    Returns:
        Tuple of (annotated_df, exploded_df, stats_summary_df).
    """
    input_path = Path(input_path).resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    logger.info(f"Loading sequence motif dataset from: {input_path}")
    df = load_dataset(input_path)
    logger.info(f"Loaded {len(df)} records.")

    motif_col = "Motif"
    if motif_col not in df.columns:
        raise ValueError(f"Required column '{motif_col}' not found in input spreadsheet.")

    # 1. Parse and extract all motifs per row
    logger.info("Extracting motif ranges and names...")
    df["Extracted_Motifs"] = df[motif_col].apply(extract_all_motifs)
    df["Motif_Range"] = df["Extracted_Motifs"].apply(lambda lst: "; ".join([m["Motif_Range"] for m in lst if m["Motif_Range"]]))
    df["Motif_Name"] = df["Extracted_Motifs"].apply(lambda lst: "; ".join([m["Motif_Name"] for m in lst if m["Motif_Name"]]))

    # Filter cleaned annotated dataframe
    annotated_df = df.drop(columns=["Extracted_Motifs", "Unnamed: 0"], errors="ignore").copy()

    # 2. Explode individual motif occurrences for statistical analysis
    has_motifs = df[df["Extracted_Motifs"].apply(len) > 0].copy()
    meta_cols = [c for c in ["Entry", "Entry Name", "Protein names", "Gene Names", "Organism"] if c in df.columns]

    exploded_df = has_motifs[meta_cols + ["Extracted_Motifs"]].explode("Extracted_Motifs").reset_index(drop=True)
    exploded_df["Motif_Range"] = exploded_df["Extracted_Motifs"].apply(lambda d: d.get("Motif_Range", "") if isinstance(d, dict) else "")
    exploded_df["Motif_Name"] = exploded_df["Extracted_Motifs"].apply(lambda d: d.get("Motif_Name", "") if isinstance(d, dict) else "")
    exploded_df = exploded_df.drop(columns=["Extracted_Motifs"])

    # 3. Aggregate statistics per unique motif name
    stats_df = (
        exploded_df[exploded_df["Motif_Name"] != ""]
        .groupby("Motif_Name", as_index=False)
        .agg(
            Number_of_Entries=("Entry", "nunique"),
            Total_Occurrences=("Entry", "count"),
            Uniprot_IDs=("Entry", lambda x: ", ".join(sorted(x.unique())))
        )
        .sort_values(by=["Number_of_Entries", "Total_Occurrences"], ascending=[False, False])
        .reset_index(drop=True)
    )

    logger.info(f"Identified {len(stats_df)} unique motifs across {exploded_df['Entry'].nunique()} entries ({len(exploded_df)} total motif instances).")
    logger.info("Top 5 most frequent motifs:")
    for idx, row in stats_df.head(5).iterrows():
        logger.info(f"  {idx + 1}. {row['Motif_Name']} : {row['Number_of_Entries']} entries ({row['Total_Occurrences']} occurrences)")

    # 4. Save to Multi-Sheet Excel Workbook
    if output_path:
        output_path = Path(output_path).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
            stats_df.to_excel(writer, sheet_name="Motif_Statistics", index=False)
            exploded_df.to_excel(writer, sheet_name="Individual_Motifs", index=False)
            annotated_df.to_excel(writer, sheet_name="Annotated_Dataset", index=False)
        logger.info(f"Saved final multi-sheet Excel workbook -> {output_path}")

    # 5. Save structured JSON output
    if json_output_path:
        json_output_path = Path(json_output_path).resolve()
        json_output_path.parent.mkdir(parents=True, exist_ok=True)
        json_data = {
            "motif_statistics": stats_df.to_dict(orient="records"),
            "individual_motifs": exploded_df.to_dict(orient="records")
        }
        with open(json_output_path, "w", encoding="utf-8") as jf:
            json.dump(json_data, jf, indent=2, ensure_ascii=False)
        logger.info(f"Saved structured JSON output -> {json_output_path}")

    # 6. Export summary CSV if requested
    if csv_output_path:
        csv_output_path = Path(csv_output_path).resolve()
        csv_output_path.parent.mkdir(parents=True, exist_ok=True)
        stats_df.to_csv(csv_output_path, index=False)
        logger.info(f"Saved motif statistics CSV -> {csv_output_path}")

    return annotated_df, exploded_df, stats_df


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract sequence motifs and compile statistical summaries into Excel, JSON, and CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_file",
        default=str(DEFAULT_INPUT_FILE),
        help="Path to input Excel or TSV/CSV dataset containing 'Motif' annotations."
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
        help="Optional path to export motif statistics summary table to CSV."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_file)
    output_path = Path(args.output_file) if args.output_file else None
    json_path = Path(args.json_output) if args.json_output else None
    csv_path = Path(args.csv_file) if args.csv_file else None

    try:
        process_motif_statistics(
            input_path=input_path,
            output_path=output_path,
            json_output_path=json_path,
            csv_output_path=csv_path
        )
    except Exception as e:
        logger.error(f"Error during motif statistics processing: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
