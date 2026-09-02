"""
generate_supplementary_json.py
--------------------------------
Generates JSON files for all tables in the manuscript Supplementary Excel workbook:
  - Input:  PhD_projects/TF-NRD/Supplementary/Supplementary.xlsx
  - Output: PhD_projects/TF-NRD/Supplementary/<Sheet_Name>.json

Parses all 17 supplementary sheets ('Table S1.A' to 'Table S15'), cleans column
headers, handles multi-level headers where applicable, converts rows into clean JSON
records, and exports them directly into the Supplementary directory.
"""

from pathlib import Path
import argparse
import logging
import json
import sys
import pandas as pd

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
PROJECT_DIR = SCRIPT_DIR.parent
SUPPLEMENTARY_DIR = PROJECT_DIR / "Supplementary"
DEFAULT_EXCEL = SUPPLEMENTARY_DIR / "Supplementary.xlsx"

# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(SCRIPT_PATH.stem)


def convert_excel_sheets_to_json(excel_file: Path = DEFAULT_EXCEL, output_dir: Path = SUPPLEMENTARY_DIR) -> dict:
    """
    Parses all sheets in Supplementary.xlsx and writes JSON files named after each sheet into output_dir.

    Returns:
        dict: Mapping of sheet names to generated JSON file paths.
    """
    excel_file = Path(excel_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not excel_file.exists():
        logger.error(f"Supplementary Excel file not found: {excel_file}")
        return {}

    logger.info(f"Opening Supplementary Excel workbook: {excel_file}")
    xls = pd.ExcelFile(excel_file)
    generated_files = {}

    for sheet in xls.sheet_names:
        logger.info(f"Processing sheet: '{sheet}'...")
        df_raw = pd.read_excel(xls, sheet_name=sheet, header=None)

        # Custom header handling for multi-level or titled sheets
        if sheet == 'Table S1.B':
            headers = [
                'PDB_ID', 'Protein_Name', 'Source Organism', 'Protein Length',
                'BSA_Protein', 'BSA_DNA', 'BSA_RNA', 'BSA_Complex',
                'FNP_Protein', 'FNP_DNA', 'FNP_RNA'
            ]
            df_data = df_raw.iloc[2:].copy()
        elif sheet == 'Table S5':
            headers = [
                'Experimental Methods',
                'Number of PDB IDs (TF-DNA)',
                'Number of PDB IDs (TF-RNA)',
                'Number of PDB IDs (TF-DNA-RNA)'
            ]
            df_data = df_raw.iloc[3:].copy()
        else:
            header_idx = 1 if sheet in [
                'Table S1.A', 'Table S1.C', 'Table S2', 'Table S3', 'Table S6',
                'Table S7', 'Table S8', 'Table S9', 'Table S10', 'Table S11',
                'Table S12', 'Table 13', 'Table S14', 'Table S15'
            ] else 0
            headers = [str(h).strip() for h in df_raw.iloc[header_idx].values]
            df_data = df_raw.iloc[header_idx + 1:].copy()

        # Deduplicate headers if needed
        seen = {}
        unique_headers = []
        for h in headers:
            if h == 'nan' or h == '':
                h = 'unnamed'
            if h in seen:
                seen[h] += 1
                unique_headers.append(f"{h}_{seen[h]}")
            else:
                seen[h] = 0
                unique_headers.append(h)

        df_data.columns = unique_headers
        df_data = df_data.dropna(how='all')
        df_data = df_data.loc[:, ~df_data.columns.str.startswith('unnamed')]

        records = []
        for _, row in df_data.iterrows():
            cleaned_row = {}
            for col in df_data.columns:
                val = row[col]
                if pd.isna(val) or val is None or str(val).strip() == '':
                    cleaned_row[col] = None
                else:
                    if isinstance(val, (int, float, str, bool)):
                        cleaned_row[col] = val
                    else:
                        cleaned_row[col] = str(val)
            if any(v is not None for v in cleaned_row.values()):
                records.append(cleaned_row)

        json_path = output_dir / f"{sheet}.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(records, f, indent=2, ensure_ascii=False)

        logger.info(f"Successfully generated JSON: {json_path.name} ({len(records)} records)")
        generated_files[sheet] = json_path

    logger.info(f"Completed JSON conversion for all {len(generated_files)} sheets into {output_dir}")
    return generated_files


def main():
    parser = argparse.ArgumentParser(description="Convert Supplementary XLSX sheets to JSON files.")
    parser.add_argument(
        "-i", "--input-file",
        default=str(DEFAULT_EXCEL),
        help="Path to Supplementary.xlsx workbook"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=str(SUPPLEMENTARY_DIR),
        help="Directory to save JSON files"
    )

    args = parser.parse_args()
    convert_excel_sheets_to_json(excel_file=Path(args.input_file), output_dir=Path(args.output_dir))


if __name__ == "__main__":
    main()
