"""
generate_webpages.py
--------------------
Generates modern, interactive, responsive HTML webpages for the TF-NRD database
from the master supplementary dataset:
  - Input:  PhD_projects/TF-NRD/Supplementary/Supplementary.xlsx
  - Target Pages:
      1. tfnrd.html         -> Database Landing / Home Page
      2. TFNRDv1.0.html     -> Table S1.A: Biological assembly of TFs & interfaces
      3. TFNRDv1.0_PNA.html -> Table S1.B: Unique TF-NA interface dataset
      4. TFNRDv1.0_PP.html  -> Table S1.C: Protein-Protein interfaces of TFs
      5. css/style.css      -> Modern, responsive, unified stylesheet

Features:
  - High-performance client-side instant search across all columns
  - Interactive category filtering (All, DNA-binding, RNA-binding, DNA-RNA-binding)
  - Interactive multi-column sorting (numeric & alphabetic)
  - Data export (CSV, JSON, Copy to Clipboard)
  - Click-to-inspect detail modal for PDB complexes with direct RCSB PDB 3D links
  - Mobile-responsive layout, glassmorphism navigation, modern typography
"""

from pathlib import Path
import json
import logging
import sys
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
SUPPLEMENTARY_DIR = PROJECT_DIR / "Supplementary"
EXCEL_PATH = SUPPLEMENTARY_DIR / "Supplementary.xlsx"

# Output target directories
TARGET_DIRS = [
    PROJECT_DIR / "TFNRDv1.0",
    PROJECT_DIR
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("generate_webpages")


def parse_table_s1a(xls: pd.ExcelFile) -> list:
    """Parses Table S1.A: Biological assembly of the TFs and their interface."""
    logger.info("Parsing Table S1.A...")
    df = pd.read_excel(xls, sheet_name='Table S1.A', header=None)
    records = []
    current_cat = 'DNA-binding'

    for idx, row in df.iterrows():
        val0 = str(row[0]).strip() if pd.notna(row[0]) else ''
        if val0.startswith('#'):
            current_cat = val0.replace('#', '').replace('Transcription Factors', '').strip()
            continue
        if idx <= 1 or val0 == '' or val0.startswith('PDB') or val0.startswith('Table'):
            continue

        def s(v):
            if pd.isna(v) or str(v).strip() in ['', 'nan', 'None']:
                return '-'
            val_str = str(v).strip()
            # Clean floating point trailing zeros if appropriate
            try:
                flt = float(val_str)
                if flt.is_integer():
                    return str(int(flt))
                return f"{flt:.2f}"
            except ValueError:
                return val_str

        rec = {
            "pdb_id": str(row[0]).strip(),
            "protein_name": s(row[1]),
            "source_organism": s(row[2]),
            "protein_length": s(row[3]),
            "pdb_chain_id": s(row[4]),
            "bio_assembly": s(row[5]),
            "exp_method": s(row[6]),
            "protein_entities": s(row[7]),
            "rna_entities": s(row[8]),
            "dna_entities": s(row[9]),
            "resolution": s(row[10]),
            "oligomeric_state": s(row[11]),
            "bsa_complex": s(row[12]),
            "bsa_protein": s(row[13]),
            "bsa_dna": s(row[14]),
            "bsa_rna": s(row[15]),
            "category": current_cat
        }
        records.append(rec)

    logger.info(f"Loaded {len(records)} records from Table S1.A")
    return records


def parse_table_s1b(xls: pd.ExcelFile) -> list:
    """Parses Table S1.B: Unique interface of the TFs."""
    logger.info("Parsing Table S1.B...")
    df = pd.read_excel(xls, sheet_name='Table S1.B', header=None)
    records = []
    current_cat = 'DNA-binding'

    for idx, row in df.iterrows():
        val0 = str(row[0]).strip() if pd.notna(row[0]) else ''
        if val0.startswith('#'):
            current_cat = val0.replace('#', '').replace('Transcription Factors', '').strip()
            continue
        if idx <= 1 or val0 == '' or val0.startswith('PDB') or val0.startswith('Table'):
            continue

        def s(v):
            if pd.isna(v) or str(v).strip() in ['', 'nan', 'None']:
                return '-'
            val_str = str(v).strip()
            try:
                flt = float(val_str)
                if flt.is_integer():
                    return str(int(flt))
                return f"{flt:.2f}"
            except ValueError:
                return val_str

        rec = {
            "pdb_id": str(row[0]).strip(),
            "protein_name": s(row[1]),
            "source_organism": s(row[2]),
            "protein_length": s(row[3]),
            "chain_protein": s(row[4]),
            "chain_dna": s(row[5]),
            "chain_rna": s(row[6]),
            "bsa_complex": s(row[7]),
            "bsa_protein": s(row[8]),
            "bsa_dna": s(row[9]),
            "bsa_rna": s(row[10]),
            "category": current_cat
        }
        records.append(rec)

    logger.info(f"Loaded {len(records)} records from Table S1.B")
    return records


def parse_table_s1c(xls: pd.ExcelFile, s1a_records: list) -> list:
    """Parses Table S1.C: Protein-Protein interface of the TFs."""
    logger.info("Parsing Table S1.C...")
    df = pd.read_excel(xls, sheet_name='Table S1.C', header=None)

    # Build PDB lookup from S1.A for metadata augmentation
    pdb_meta = {}
    for r in s1a_records:
        pdb_meta[r["pdb_id"]] = {
            "protein_name": r["protein_name"],
            "source_organism": r["source_organism"],
            "exp_method": r["exp_method"],
            "resolution": r["resolution"]
        }

    records = []
    current_cat = 'DNA-binding'

    for idx, row in df.iterrows():
        val0 = str(row[0]).strip() if pd.notna(row[0]) else ''
        if val0.startswith('#'):
            current_cat = val0.replace('#', '').replace('Transcription Factors', '').strip()
            continue
        if idx <= 1 or val0 == '' or val0.startswith('PDB') or val0.startswith('Table'):
            continue

        def s(v):
            if pd.isna(v) or str(v).strip() in ['', 'nan', 'None']:
                return '-'
            val_str = str(v).strip()
            try:
                flt = float(val_str)
                if flt.is_integer():
                    return str(int(flt))
                return f"{flt:.2f}"
            except ValueError:
                return val_str

        pdb_id = str(row[0]).strip()
        meta = pdb_meta.get(pdb_id, {})

        rec = {
            "pdb_id": pdb_id,
            "chain1": s(row[1]),
            "chain2": s(row[2]),
            "bsa_complex": s(row[3]),
            "bsa_protein1": s(row[4]),
            "bsa_protein2": s(row[5]),
            "category": current_cat,
            "protein_name": meta.get("protein_name", "-"),
            "source_organism": meta.get("source_organism", "-")
        }
        records.append(rec)

    logger.info(f"Loaded {len(records)} records from Table S1.C")
    return records


def get_css_content() -> str:
    """Returns the unified modern stylesheet content."""
    return """/* ==========================================================================
   TF-NRD: Transcription Factors Non-Redundant Dataset
   Modern, High-Performance Scientific Database Stylesheet
   ========================================================================== */

@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=Plus+Jakarta+Sans:wght@500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap');

:root {
  --primary: #2563eb;
  --primary-hover: #1d4ed8;
  --primary-light: #eff6ff;
  --primary-border: #bfdbfe;
  
  --secondary: #0ea5e9;
  --accent: #6366f1;
  
  --cat-dna-bg: #ecfdf5;
  --cat-dna-color: #047857;
  --cat-dna-border: #a7f3d0;
  
  --cat-rna-bg: #f5f3ff;
  --cat-rna-color: #6d28d9;
  --cat-rna-border: #ddd6fe;
  
  --cat-dnarna-bg: #fffbeb;
  --cat-dnarna-color: #b45309;
  --cat-dnarna-border: #fde68a;

  --bg-main: #f8fafc;
  --bg-card: #ffffff;
  --bg-subtle: #f1f5f9;
  
  --text-main: #0f172a;
  --text-muted: #64748b;
  --text-light: #94a3b8;
  
  --border-color: #e2e8f0;
  --border-subtle: #f1f5f9;
  
  --radius-sm: 6px;
  --radius-md: 10px;
  --radius-lg: 16px;
  --radius-full: 9999px;
  
  --shadow-sm: 0 1px 2px 0 rgba(0, 0, 0, 0.05);
  --shadow-md: 0 4px 6px -1px rgba(0, 0, 0, 0.07), 0 2px 4px -2px rgba(0, 0, 0, 0.05);
  --shadow-lg: 0 10px 15px -3px rgba(0, 0, 0, 0.08), 0 4px 6px -4px rgba(0, 0, 0, 0.04);
  --shadow-xl: 0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 8px 10px -6px rgba(0, 0, 0, 0.04);

  --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  --font-heading: 'Plus Jakarta Sans', var(--font-sans);
  --font-mono: 'JetBrains Mono', monospace;
}

*, *::before, *::after {
  box-sizing: border-box;
  margin: 0;
  padding: 0;
}

body {
  font-family: var(--font-sans);
  background-color: var(--bg-main);
  color: var(--text-main);
  line-height: 1.5;
  min-height: 100vh;
  display: flex;
  flex-direction: column;
  -webkit-font-smoothing: antialiased;
}

/* ==========================================================================
   NAVIGATION & HEADER
   ========================================================================== */
header.site-header {
  background: #0f172a;
  color: #ffffff;
  position: sticky;
  top: 0;
  z-index: 1000;
  border-bottom: 1px solid rgba(255, 255, 255, 0.1);
  box-shadow: 0 4px 20px rgba(15, 23, 42, 0.3);
  backdrop-filter: blur(8px);
}

.header-container {
  max-width: 1440px;
  margin: 0 auto;
  padding: 0 24px;
  height: 72px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 20px;
}

.brand {
  display: flex;
  align-items: center;
  gap: 12px;
  text-decoration: none;
  color: #ffffff;
}

.brand-icon {
  width: 40px;
  height: 40px;
  background: linear-gradient(135deg, #3b82f6, #06b6d4);
  border-radius: var(--radius-md);
  display: flex;
  align-items: center;
  justify-content: center;
  font-weight: 800;
  font-size: 18px;
  color: white;
  box-shadow: 0 2px 10px rgba(59, 130, 246, 0.4);
}

.brand-text h1 {
  font-family: var(--font-heading);
  font-size: 20px;
  font-weight: 700;
  letter-spacing: -0.5px;
  margin: 0;
  display: flex;
  align-items: center;
  gap: 8px;
}

.brand-text .version-tag {
  font-size: 11px;
  font-weight: 700;
  background: rgba(59, 130, 246, 0.25);
  color: #60a5fa;
  padding: 2px 8px;
  border-radius: var(--radius-full);
  border: 1px solid rgba(96, 165, 250, 0.3);
}

.brand-text p {
  font-size: 12px;
  color: #94a3b8;
  margin: 0;
}

nav.main-nav {
  display: flex;
  align-items: center;
  gap: 8px;
}

nav.main-nav a {
  color: #cbd5e1;
  text-decoration: none;
  font-size: 14px;
  font-weight: 500;
  padding: 8px 14px;
  border-radius: var(--radius-sm);
  transition: all 0.2s ease;
  display: flex;
  align-items: center;
  gap: 6px;
}

nav.main-nav a:hover {
  color: #ffffff;
  background: rgba(255, 255, 255, 0.08);
}

nav.main-nav a.active {
  color: #ffffff;
  background: #2563eb;
  font-weight: 600;
  box-shadow: 0 2px 8px rgba(37, 99, 235, 0.4);
}

.header-actions {
  display: flex;
  align-items: center;
  gap: 12px;
}

.btn-header-link {
  color: #cbd5e1;
  text-decoration: none;
  font-size: 13px;
  display: flex;
  align-items: center;
  gap: 6px;
  padding: 6px 12px;
  border: 1px solid rgba(255, 255, 255, 0.15);
  border-radius: var(--radius-sm);
  transition: all 0.2s;
}

.btn-header-link:hover {
  background: rgba(255, 255, 255, 0.1);
  color: #ffffff;
}

/* ==========================================================================
   PAGE HERO / BANNER
   ========================================================================== */
.page-hero {
  background: linear-gradient(135deg, #0f172a 0%, #1e293b 100%);
  color: #ffffff;
  padding: 40px 0 36px 0;
  border-bottom: 1px solid var(--border-color);
  position: relative;
  overflow: hidden;
}

.page-hero::after {
  content: '';
  position: absolute;
  top: -50%;
  right: -10%;
  width: 500px;
  height: 500px;
  background: radial-gradient(circle, rgba(59, 130, 246, 0.15) 0%, rgba(15, 23, 42, 0) 70%);
  pointer-events: none;
}

.container {
  max-width: 1440px;
  margin: 0 auto;
  padding: 0 24px;
}

.hero-content h2 {
  font-family: var(--font-heading);
  font-size: 30px;
  font-weight: 800;
  letter-spacing: -0.75px;
  margin-bottom: 8px;
  display: flex;
  align-items: center;
  gap: 12px;
  flex-wrap: wrap;
}

.hero-content p {
  color: #94a3b8;
  font-size: 15px;
  max-width: 850px;
  line-height: 1.6;
}

.hero-stats-row {
  display: flex;
  gap: 16px;
  margin-top: 24px;
  flex-wrap: wrap;
}

.hero-stat-pill {
  background: rgba(255, 255, 255, 0.06);
  border: 1px solid rgba(255, 255, 255, 0.1);
  padding: 8px 16px;
  border-radius: var(--radius-md);
  display: flex;
  align-items: center;
  gap: 10px;
}

.hero-stat-pill .number {
  font-family: var(--font-heading);
  font-size: 18px;
  font-weight: 700;
  color: #38bdf8;
}

.hero-stat-pill .label {
  font-size: 13px;
  color: #cbd5e1;
}

/* ==========================================================================
   DATA EXPLORER / CONTROLS BAR
   ========================================================================== */
.controls-section {
  padding: 24px 0 16px 0;
}

.controls-card {
  background: var(--bg-card);
  border: 1px solid var(--border-color);
  border-radius: var(--radius-lg);
  padding: 20px 24px;
  box-shadow: var(--shadow-sm);
  display: flex;
  flex-direction: column;
  gap: 16px;
}

.controls-top {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
  flex-wrap: wrap;
}

.search-box-wrapper {
  position: relative;
  flex: 1;
  min-width: 300px;
  max-width: 500px;
}

.search-icon {
  position: absolute;
  left: 14px;
  top: 50%;
  transform: translateY(-50%);
  color: var(--text-light);
  pointer-events: none;
}

input.search-input {
  width: 100%;
  height: 44px;
  padding: 0 40px 0 42px;
  font-family: var(--font-sans);
  font-size: 14px;
  border: 1.5px solid var(--border-color);
  border-radius: var(--radius-md);
  background: var(--bg-main);
  color: var(--text-main);
  transition: all 0.2s ease;
}

input.search-input:focus {
  outline: none;
  border-color: var(--primary);
  background: #ffffff;
  box-shadow: 0 0 0 4px rgba(37, 99, 235, 0.12);
}

.clear-search-btn {
  position: absolute;
  right: 12px;
  top: 50%;
  transform: translateY(-50%);
  background: none;
  border: none;
  color: var(--text-light);
  cursor: pointer;
  font-size: 16px;
  display: none;
  padding: 4px;
}

.clear-search-btn:hover {
  color: var(--text-main);
}

.action-buttons-group {
  display: flex;
  align-items: center;
  gap: 10px;
  flex-wrap: wrap;
}

.btn-action {
  font-family: var(--font-sans);
  font-size: 13px;
  font-weight: 600;
  height: 40px;
  padding: 0 16px;
  border-radius: var(--radius-md);
  display: inline-flex;
  align-items: center;
  gap: 8px;
  cursor: pointer;
  transition: all 0.2s ease;
  text-decoration: none;
  border: 1px solid var(--border-color);
  background: #ffffff;
  color: var(--text-main);
}

.btn-action:hover {
  background: var(--bg-subtle);
  border-color: #cbd5e1;
}

.btn-action.btn-primary {
  background: var(--primary);
  color: #ffffff;
  border-color: var(--primary);
}

.btn-action.btn-primary:hover {
  background: var(--primary-hover);
}

/* Category Filter Tabs */
.filter-tabs {
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
  align-items: center;
  border-top: 1px solid var(--border-subtle);
  padding-top: 16px;
}

.filter-tab-label {
  font-size: 13px;
  font-weight: 600;
  color: var(--text-muted);
  margin-right: 4px;
}

.filter-chip {
  font-family: var(--font-sans);
  font-size: 13px;
  font-weight: 500;
  padding: 6px 14px;
  border-radius: var(--radius-full);
  border: 1.5px solid var(--border-color);
  background: #ffffff;
  color: var(--text-muted);
  cursor: pointer;
  display: inline-flex;
  align-items: center;
  gap: 8px;
  transition: all 0.2s ease;
}

.filter-chip:hover {
  background: var(--bg-subtle);
  color: var(--text-main);
}

.filter-chip.active {
  background: #0f172a;
  color: #ffffff;
  border-color: #0f172a;
  font-weight: 600;
}

.filter-chip .chip-count {
  font-size: 11px;
  background: rgba(0, 0, 0, 0.08);
  padding: 1px 7px;
  border-radius: var(--radius-full);
  font-weight: 700;
}

.filter-chip.active .chip-count {
  background: rgba(255, 255, 255, 0.2);
  color: #ffffff;
}

/* ==========================================================================
   DATA TABLE CONTAINER & STYLES
   ========================================================================== */
.table-section {
  padding-bottom: 40px;
}

.table-card {
  background: var(--bg-card);
  border: 1px solid var(--border-color);
  border-radius: var(--radius-lg);
  box-shadow: var(--shadow-sm);
  overflow: hidden;
  display: flex;
  flex-direction: column;
}

.table-header-info {
  padding: 16px 24px;
  background: #ffffff;
  border-bottom: 1px solid var(--border-color);
  display: flex;
  justify-content: space-between;
  align-items: center;
  flex-wrap: wrap;
  gap: 12px;
}

.results-count {
  font-size: 14px;
  color: var(--text-muted);
}

.results-count strong {
  color: var(--text-main);
}

.page-size-selector {
  display: flex;
  align-items: center;
  gap: 8px;
  font-size: 13px;
  color: var(--text-muted);
}

.page-size-selector select {
  font-family: var(--font-sans);
  font-size: 13px;
  font-weight: 500;
  padding: 4px 8px;
  border-radius: var(--radius-sm);
  border: 1px solid var(--border-color);
  background: #ffffff;
  color: var(--text-main);
  cursor: pointer;
}

.responsive-table-wrapper {
  width: 100%;
  overflow-x: auto;
  -webkit-overflow-scrolling: touch;
  max-height: 75vh;
  position: relative;
}

table.data-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
  text-align: left;
}

table.data-table thead {
  position: sticky;
  top: 0;
  z-index: 20;
  background: #f8fafc;
}

table.data-table th {
  padding: 12px 14px;
  font-family: var(--font-heading);
  font-weight: 700;
  font-size: 12px;
  color: #334155;
  text-transform: uppercase;
  letter-spacing: 0.5px;
  border-bottom: 2px solid var(--border-color);
  border-right: 1px solid var(--border-subtle);
  background: #f8fafc;
  user-select: none;
  white-space: nowrap;
  cursor: pointer;
  transition: background 0.15s ease;
}

table.data-table th:hover {
  background: #edf2f7;
  color: var(--primary);
}

table.data-table th.sort-active {
  color: var(--primary);
  background: #eff6ff;
}

table.data-table th .sort-indicator {
  display: inline-block;
  margin-left: 6px;
  color: var(--text-light);
  font-size: 11px;
}

table.data-table th.sort-active .sort-indicator {
  color: var(--primary);
}

table.data-table td {
  padding: 10px 14px;
  border-bottom: 1px solid var(--border-subtle);
  border-right: 1px solid var(--border-subtle);
  vertical-align: middle;
  color: #1e293b;
}

table.data-table tbody tr {
  transition: background 0.1s ease;
  cursor: pointer;
}

table.data-table tbody tr:hover {
  background: #f1f5f9 !important;
}

table.data-table tbody tr:nth-child(even) {
  background: #fafcff;
}

/* Category Badges */
.badge-category {
  display: inline-flex;
  align-items: center;
  gap: 5px;
  padding: 3px 8px;
  border-radius: var(--radius-full);
  font-size: 11px;
  font-weight: 600;
  white-space: nowrap;
}

.badge-category.cat-dna {
  background: var(--cat-dna-bg);
  color: var(--cat-dna-color);
  border: 1px solid var(--cat-dna-border);
}

.badge-category.cat-rna {
  background: var(--cat-rna-bg);
  color: var(--cat-rna-color);
  border: 1px solid var(--cat-rna-border);
}

.badge-category.cat-dnarna {
  background: var(--cat-dnarna-bg);
  color: var(--cat-dnarna-color);
  border: 1px solid var(--cat-dnarna-border);
}

/* PDB ID Link Badge */
.pdb-badge {
  display: inline-flex;
  align-items: center;
  gap: 4px;
  font-family: var(--font-mono);
  font-size: 12px;
  font-weight: 700;
  color: var(--primary);
  background: var(--primary-light);
  border: 1px solid var(--primary-border);
  padding: 2px 7px;
  border-radius: var(--radius-sm);
  text-decoration: none;
  transition: all 0.15s ease;
}

.pdb-badge:hover {
  background: var(--primary);
  color: #ffffff;
  transform: translateY(-1px);
  box-shadow: 0 2px 6px rgba(37, 99, 235, 0.3);
}

/* Text Wrapping & Highlighting */
.cell-protein-name {
  max-width: 280px;
  white-space: normal;
  word-break: break-word;
  font-weight: 500;
}

.cell-organism {
  font-style: italic;
  color: #475569;
  white-space: nowrap;
}

.cell-mono {
  font-family: var(--font-mono);
  font-size: 12px;
}

.cell-num {
  font-family: var(--font-mono);
  text-align: right;
  font-variant-numeric: tabular-nums;
}

/* Method Tag */
.method-pill {
  font-size: 11px;
  font-weight: 600;
  padding: 2px 7px;
  border-radius: var(--radius-sm);
  background: #f1f5f9;
  color: #475569;
  border: 1px solid #e2e8f0;
}

/* ==========================================================================
   PAGINATION CONTROLS
   ========================================================================== */
.pagination-container {
  padding: 16px 24px;
  background: #ffffff;
  border-top: 1px solid var(--border-color);
  display: flex;
  align-items: center;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 16px;
}

.pagination-info {
  font-size: 13px;
  color: var(--text-muted);
}

.pagination-buttons {
  display: flex;
  align-items: center;
  gap: 4px;
}

.btn-page {
  font-family: var(--font-sans);
  font-size: 13px;
  font-weight: 600;
  min-width: 36px;
  height: 36px;
  padding: 0 8px;
  border: 1px solid var(--border-color);
  background: #ffffff;
  color: var(--text-main);
  border-radius: var(--radius-sm);
  cursor: pointer;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  transition: all 0.15s;
}

.btn-page:hover:not(:disabled) {
  background: var(--bg-subtle);
  border-color: #cbd5e1;
}

.btn-page.active {
  background: var(--primary);
  color: #ffffff;
  border-color: var(--primary);
}

.btn-page:disabled {
  opacity: 0.4;
  cursor: not-allowed;
}

/* ==========================================================================
   MODAL / DETAILS POPUP
   ========================================================================== */
.modal-overlay {
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: rgba(15, 23, 42, 0.6);
  backdrop-filter: blur(4px);
  z-index: 2000;
  display: none;
  align-items: center;
  justify-content: center;
  padding: 20px;
}

.modal-overlay.open {
  display: flex;
  animation: fadeIn 0.2s ease-out;
}

.modal-card {
  background: #ffffff;
  border-radius: var(--radius-lg);
  max-width: 700px;
  width: 100%;
  max-height: 90vh;
  overflow-y: auto;
  box-shadow: var(--shadow-xl);
  border: 1px solid var(--border-color);
  display: flex;
  flex-direction: column;
}

.modal-header {
  padding: 20px 24px;
  border-bottom: 1px solid var(--border-color);
  display: flex;
  align-items: center;
  justify-content: space-between;
  background: #f8fafc;
}

.modal-header h3 {
  font-family: var(--font-heading);
  font-size: 18px;
  font-weight: 700;
  display: flex;
  align-items: center;
  gap: 10px;
}

.modal-close-btn {
  background: none;
  border: none;
  font-size: 22px;
  color: var(--text-muted);
  cursor: pointer;
  padding: 4px;
  line-height: 1;
}

.modal-close-btn:hover {
  color: var(--text-main);
}

.modal-body {
  padding: 24px;
  display: flex;
  flex-direction: column;
  gap: 16px;
}

.modal-detail-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 12px;
}

.modal-detail-item {
  background: var(--bg-subtle);
  padding: 10px 14px;
  border-radius: var(--radius-sm);
  border: 1px solid var(--border-subtle);
}

.modal-detail-item .detail-label {
  font-size: 11px;
  font-weight: 700;
  text-transform: uppercase;
  color: var(--text-muted);
  letter-spacing: 0.5px;
  margin-bottom: 4px;
}

.modal-detail-item .detail-val {
  font-size: 14px;
  font-weight: 600;
  color: var(--text-main);
  word-break: break-word;
}

.modal-footer {
  padding: 16px 24px;
  border-top: 1px solid var(--border-color);
  display: flex;
  justify-content: flex-end;
  gap: 12px;
  background: #f8fafc;
}

/* ==========================================================================
   LANDING / HOME PAGE SPECIFIC STYLES (tfnrd.html)
   ========================================================================== */
.home-hero {
  background: linear-gradient(135deg, #0f172a 0%, #1e293b 60%, #1e1b4b 100%);
  color: #ffffff;
  padding: 70px 0 60px 0;
  text-align: center;
  position: relative;
  overflow: hidden;
}

.home-hero::before {
  content: '';
  position: absolute;
  top: -50%;
  left: 20%;
  width: 600px;
  height: 600px;
  background: radial-gradient(circle, rgba(59, 130, 246, 0.25) 0%, rgba(15, 23, 42, 0) 70%);
  pointer-events: none;
}

.home-hero h1 {
  font-family: var(--font-heading);
  font-size: 44px;
  font-weight: 800;
  letter-spacing: -1.5px;
  margin-bottom: 16px;
  line-height: 1.15;
}

.home-hero h1 span.highlight {
  background: linear-gradient(135deg, #60a5fa, #38bdf8);
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
}

.home-hero p.hero-subtitle {
  font-size: 18px;
  color: #cbd5e1;
  max-width: 800px;
  margin: 0 auto 32px auto;
  line-height: 1.6;
}

.home-cta-buttons {
  display: flex;
  justify-content: center;
  gap: 16px;
  flex-wrap: wrap;
  margin-bottom: 40px;
}

.btn-hero-primary {
  background: #2563eb;
  color: #ffffff;
  font-family: var(--font-heading);
  font-size: 15px;
  font-weight: 700;
  padding: 14px 28px;
  border-radius: var(--radius-md);
  text-decoration: none;
  display: inline-flex;
  align-items: center;
  gap: 8px;
  box-shadow: 0 4px 14px rgba(37, 99, 235, 0.5);
  transition: all 0.2s;
}

.btn-hero-primary:hover {
  background: #1d4ed8;
  transform: translateY(-2px);
  box-shadow: 0 6px 20px rgba(37, 99, 235, 0.6);
}

.btn-hero-secondary {
  background: rgba(255, 255, 255, 0.1);
  color: #ffffff;
  font-family: var(--font-heading);
  font-size: 15px;
  font-weight: 600;
  padding: 14px 28px;
  border-radius: var(--radius-md);
  text-decoration: none;
  display: inline-flex;
  align-items: center;
  gap: 8px;
  border: 1px solid rgba(255, 255, 255, 0.2);
  transition: all 0.2s;
}

.btn-hero-secondary:hover {
  background: rgba(255, 255, 255, 0.18);
  transform: translateY(-2px);
}

/* Home Statistics Grid */
.home-stats-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
  gap: 20px;
  margin-top: -30px;
  position: relative;
  z-index: 10;
}

.home-stat-card {
  background: #ffffff;
  border: 1px solid var(--border-color);
  border-radius: var(--radius-lg);
  padding: 24px;
  box-shadow: var(--shadow-md);
  text-align: center;
  transition: all 0.2s ease;
}

.home-stat-card:hover {
  transform: translateY(-4px);
  box-shadow: var(--shadow-lg);
  border-color: #cbd5e1;
}

.home-stat-card .stat-icon {
  width: 48px;
  height: 48px;
  margin: 0 auto 12px auto;
  border-radius: var(--radius-md);
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 24px;
}

.stat-icon.blue { background: #eff6ff; color: #2563eb; }
.stat-icon.green { background: #ecfdf5; color: #059669; }
.stat-icon.purple { background: #f5f3ff; color: #7c3aed; }
.stat-icon.orange { background: #fffbeb; color: #d97706; }

.home-stat-card h3 {
  font-family: var(--font-heading);
  font-size: 32px;
  font-weight: 800;
  color: var(--text-main);
  margin-bottom: 4px;
}

.home-stat-card p {
  font-size: 14px;
  color: var(--text-muted);
  font-weight: 500;
}

/* Feature Cards Section */
.features-section {
  padding: 60px 0;
}

.section-title-wrap {
  text-align: center;
  max-width: 700px;
  margin: 0 auto 40px auto;
}

.section-title-wrap h2 {
  font-family: var(--font-heading);
  font-size: 28px;
  font-weight: 800;
  letter-spacing: -0.5px;
  color: var(--text-main);
  margin-bottom: 8px;
}

.section-title-wrap p {
  color: var(--text-muted);
  font-size: 15px;
}

.dataset-cards-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
  gap: 24px;
}

.dataset-card {
  background: var(--bg-card);
  border: 1.5px solid var(--border-color);
  border-radius: var(--radius-lg);
  padding: 30px;
  display: flex;
  flex-direction: column;
  justify-content: space-between;
  box-shadow: var(--shadow-sm);
  transition: all 0.25s ease;
}

.dataset-card:hover {
  transform: translateY(-5px);
  box-shadow: var(--shadow-xl);
  border-color: var(--primary);
}

.dataset-card .card-top h3 {
  font-family: var(--font-heading);
  font-size: 20px;
  font-weight: 700;
  margin-bottom: 8px;
  color: var(--text-main);
}

.dataset-card .card-top p {
  font-size: 14px;
  color: var(--text-muted);
  line-height: 1.6;
  margin-bottom: 20px;
}

.dataset-card .meta-pills {
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
  margin-bottom: 24px;
}

.dataset-card .meta-pill {
  font-size: 12px;
  font-weight: 600;
  background: var(--bg-subtle);
  color: var(--text-muted);
  padding: 4px 10px;
  border-radius: var(--radius-sm);
}

.dataset-card a.btn-card {
  background: var(--primary);
  color: #ffffff;
  font-family: var(--font-heading);
  font-size: 14px;
  font-weight: 600;
  padding: 10px 18px;
  border-radius: var(--radius-md);
  text-decoration: none;
  text-align: center;
  transition: all 0.2s;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  gap: 8px;
}

.dataset-card a.btn-card:hover {
  background: var(--primary-hover);
}

/* ==========================================================================
   SITE FOOTER
   ========================================================================== */
footer.site-footer {
  background: #0f172a;
  color: #94a3b8;
  padding: 40px 0 24px 0;
  margin-top: auto;
  border-top: 1px solid rgba(255, 255, 255, 0.1);
  font-size: 14px;
}

.footer-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 30px;
  margin-bottom: 30px;
}

.footer-col h4 {
  font-family: var(--font-heading);
  font-size: 15px;
  font-weight: 700;
  color: #ffffff;
  margin-bottom: 12px;
}

.footer-col p {
  line-height: 1.6;
  font-size: 13px;
}

.footer-col a {
  color: #60a5fa;
  text-decoration: none;
}

.footer-col a:hover {
  text-decoration: underline;
}

.footer-bottom {
  border-top: 1px solid rgba(255, 255, 255, 0.08);
  padding-top: 20px;
  text-align: center;
  font-size: 12px;
  color: #64748b;
}

/* Animations */
@keyframes fadeIn {
  from { opacity: 0; }
  to { opacity: 1; }
}

/* Mobile responsive adjustments */
@media screen and (max-width: 900px) {
  .header-container {
    height: auto;
    padding: 14px 20px;
    flex-direction: column;
    align-items: flex-start;
  }
  nav.main-nav {
    flex-wrap: wrap;
    width: 100%;
  }
  .home-hero h1 {
    font-size: 32px;
  }
}
"""


def generate_tfnrd_home_html(s1a_rows: list, s1b_rows: list, s1c_rows: list) -> str:
    """Generates the modern landing/home page (tfnrd.html)."""
    n_complexes = len(s1a_rows)
    n_dna = sum(1 for r in s1a_rows if r['category'] == 'DNA-binding')
    n_rna = sum(1 for r in s1a_rows if r['category'] == 'RNA-binding')
    n_dnarna = sum(1 for r in s1a_rows if r['category'] == 'DNA-RNA-binding')
    n_pna_interfaces = len(s1b_rows)
    n_pp_interfaces = len(s1c_rows)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="description" content="TF-NRD: Transcription Factors Non-Redundant Dataset. Curated structural atlas of TF assemblies and macromolecular interaction interfaces across all domains of life.">
  <title>TF-NRD: Transcription Factors Non-Redundant Dataset</title>
  <link rel="stylesheet" href="css/style.css">
</head>
<body>

  <!-- SITE HEADER -->
  <header class="site-header">
    <div class="header-container">
      <a href="tfnrd.html" class="brand">
        <div class="brand-icon">TF</div>
        <div class="brand-text">
          <h1>TF-NRD <span class="version-tag">v1.0</span></h1>
          <p>Transcription Factors Non-Redundant Dataset</p>
        </div>
      </a>
      <nav class="main-nav" id="main-nav">
        <a href="tfnrd.html" class="active">Home</a>
        <a href="TFNRDv1.0.html">Dataset (Table S1.A)</a>
        <a href="TFNRDv1.0_PNA.html">TF-NA Interfaces (Table S1.B)</a>
        <a href="TFNRDv1.0_PP.html">Protein-Protein (Table S1.C)</a>
      </nav>
      <div class="header-actions">
        <a href="http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html" target="_blank" rel="noopener" class="btn-header-link" title="Official Webserver">
          <span>CSB Lab</span> ↗
        </a>
      </div>
    </div>
  </header>

  <!-- HERO BANNER -->
  <section class="home-hero">
    <div class="container">
      <h1>An Atlas of <span class="highlight">Non-Redundant</span> Transcription Factor Assemblies</h1>
      <p class="hero-subtitle">
        A curated structural benchmark of multi-entity transcription factor assemblies, macromolecular interaction interfaces (Protein–DNA, Protein–RNA, Protein–Protein), and sequence motifs across all domains of life.
      </p>
      <div class="home-cta-buttons">
        <a href="TFNRDv1.0.html" class="btn-hero-primary">
          <span>Explore Main Dataset</span> →
        </a>
        <a href="TFNRDv1.0_PNA.html" class="btn-hero-secondary">
          <span>TF-NA Interfaces</span>
        </a>
        <a href="TFNRDv1.0_PP.html" class="btn-hero-secondary">
          <span>Protein-Protein Interfaces</span>
        </a>
      </div>
    </div>
  </section>

  <!-- KEY METRICS GRID -->
  <section class="container">
    <div class="home-stats-grid">
      <div class="home-stat-card">
        <div class="stat-icon blue">🧬</div>
        <h3>{n_complexes}</h3>
        <p>Non-Redundant Structure Complexes</p>
      </div>
      <div class="home-stat-card">
        <div class="stat-icon green">📊</div>
        <h3>{n_pna_interfaces}</h3>
        <p>Unique TF-NA Interaction Interfaces</p>
      </div>
      <div class="home-stat-card">
        <div class="stat-icon purple">🤝</div>
        <h3>{n_pp_interfaces}</h3>
        <p>Protein-Protein Contact Interfaces</p>
      </div>
      <div class="home-stat-card">
        <div class="stat-icon orange">🔬</div>
        <h3>3,570</h3>
        <p>Sequence-Curated TF Entries</p>
      </div>
    </div>
  </section>

  <!-- DATASETS EXPLORER SECTION -->
  <section class="features-section container">
    <div class="section-title-wrap">
      <h2>Explore TF-NRD Benchmark Datasets</h2>
      <p>Access structured datasets generated from the manuscript supplementary data with instant search, multi-column sorting, and export capabilities.</p>
    </div>

    <div class="dataset-cards-grid">
      <!-- Card 1: Table S1.A -->
      <div class="dataset-card">
        <div class="card-top">
          <span class="badge-category cat-dna" style="margin-bottom: 12px;">Table S1.A</span>
          <h3>TF Biological Assemblies & Complexes</h3>
          <p>
            Complete non-redundant structural dataset of {n_complexes} transcription factor complexes. Includes biological assemblies, experimental methods, resolutions, oligomeric states, and complex BSA values.
          </p>
          <div class="meta-pills">
            <span class="meta-pill">{n_dna} DNA-binding</span>
            <span class="meta-pill">{n_rna} RNA-binding</span>
            <span class="meta-pill">{n_dnarna} DNA-RNA-binding</span>
          </div>
        </div>
        <a href="TFNRDv1.0.html" class="btn-card">Browse Dataset Table →</a>
      </div>

      <!-- Card 2: Table S1.B -->
      <div class="dataset-card">
        <div class="card-top">
          <span class="badge-category cat-rna" style="margin-bottom: 12px;">Table S1.B</span>
          <h3>Unique TF-NA Interfaces</h3>
          <p>
            Detailed structural interfaces for {n_pna_interfaces} TF-nucleic acid interactions. Includes interacting chain IDs, nucleic acid types (DNA, RNA), and decomposed Buried Surface Area (BSA) for complex, protein, and nucleic acids.
          </p>
          <div class="meta-pills">
            <span class="meta-pill">{n_pna_interfaces} Interfaces</span>
            <span class="meta-pill">BSA Decompositions</span>
            <span class="meta-pill">Chain Pairings</span>
          </div>
        </div>
        <a href="TFNRDv1.0_PNA.html" class="btn-card">Browse TF-NA Interfaces →</a>
      </div>

      <!-- Card 3: Table S1.C -->
      <div class="dataset-card">
        <div class="card-top">
          <span class="badge-category cat-dnarna" style="margin-bottom: 12px;">Table S1.C</span>
          <h3>Protein-Protein Interfaces of TFs</h3>
          <p>
            Atomic contact interfaces for {n_pp_interfaces} protein-protein interactions within TF assemblies. Characterizes homomeric and heteromeric co-factor interactions with Chain 1 and Chain 2 surface area contributions.
          </p>
          <div class="meta-pills">
            <span class="meta-pill">{n_pp_interfaces} PP Interfaces</span>
            <span class="meta-pill">Subunit BSA</span>
            <span class="meta-pill">Oligomer Contacts</span>
          </div>
        </div>
        <a href="TFNRDv1.0_PP.html" class="btn-card">Browse Protein-Protein Interfaces →</a>
      </div>
    </div>
  </section>

  <!-- CITATION & ABOUT -->
  <section class="container" style="padding-bottom: 60px;">
    <div style="background: #ffffff; border: 1px solid var(--border-color); border-radius: var(--radius-lg); padding: 32px; box-shadow: var(--shadow-sm);">
      <h3 style="font-family: var(--font-heading); font-size: 20px; font-weight: 700; margin-bottom: 12px;">Publication & Citation</h3>
      <p style="color: var(--text-muted); font-size: 14px; line-height: 1.7; margin-bottom: 16px;">
        Garai, S., Kant, S., & Bahadur, R. P. (2026).<br>
        <strong>An atlas of non-redundant sequences and structures of transcription factor assemblies across domains of life.</strong><br>
        <em>Computational Structural Biology Group, Department of Biotechnology, Indian Institute of Technology Kharagpur, West Bengal 721302, India.</em>
      </p>
      <div style="display: flex; gap: 12px; flex-wrap: wrap;">
        <a href="http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html" target="_blank" rel="noopener" class="btn-action btn-primary">
          <span>Visit CSB Webserver</span> ↗
        </a>
      </div>
    </div>
  </section>

  <!-- FOOTER -->
  <footer class="site-footer">
    <div class="container footer-grid">
      <div class="footer-col">
        <h4>TF-NRD Database</h4>
        <p>Non-redundant structural atlas and interface benchmark for transcription factor macromolecular assemblies across domains of life.</p>
      </div>
      <div class="footer-col">
        <h4>Quick Links</h4>
        <p><a href="TFNRDv1.0.html">Dataset Table (S1.A)</a></p>
        <p><a href="TFNRDv1.0_PNA.html">TF-NA Interfaces (S1.B)</a></p>
        <p><a href="TFNRDv1.0_PP.html">Protein-Protein (S1.C)</a></p>
      </div>
      <div class="footer-col">
        <h4>Affiliation</h4>
        <p>Computational Structural Biology Group<br>Department of Biotechnology<br>Indian Institute of Technology Kharagpur<br>Kharagpur - 721302, India</p>
      </div>
    </div>
    <div class="container footer-bottom">
      <p>&copy; 2026 Computational Structural Biology Group, IIT Kharagpur. Released under MIT License.</p>
    </div>
  </footer>

</body>
</html>
"""


def generate_table_s1a_html(records: list) -> str:
    """Generates the interactive modern webpage for Table S1.A (TFNRDv1.0.html)."""
    json_data = json.dumps(records, ensure_ascii=False)
    n_total = len(records)
    n_dna = sum(1 for r in records if r['category'] == 'DNA-binding')
    n_rna = sum(1 for r in records if r['category'] == 'RNA-binding')
    n_dnarna = sum(1 for r in records if r['category'] == 'DNA-RNA-binding')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="description" content="TF-NRD Table S1.A: Non-redundant structural dataset of transcription factor-nucleic acid complexes, biological assemblies, oligomeric states, and interface BSA.">
  <title>TF-NRD v1.0 Dataset (Table S1.A)</title>
  <link rel="stylesheet" href="css/style.css">
</head>
<body>

  <!-- SITE HEADER -->
  <header class="site-header">
    <div class="header-container">
      <a href="tfnrd.html" class="brand">
        <div class="brand-icon">TF</div>
        <div class="brand-text">
          <h1>TF-NRD <span class="version-tag">Table S1.A</span></h1>
          <p>Biological Assemblies & Complex Interfaces</p>
        </div>
      </a>
      <nav class="main-nav">
        <a href="tfnrd.html">Home</a>
        <a href="TFNRDv1.0.html" class="active">Dataset (Table S1.A)</a>
        <a href="TFNRDv1.0_PNA.html">TF-NA Interfaces (Table S1.B)</a>
        <a href="TFNRDv1.0_PP.html">Protein-Protein (Table S1.C)</a>
      </nav>
      <div class="header-actions">
        <a href="http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html" target="_blank" rel="noopener" class="btn-header-link">CSB Lab ↗</a>
      </div>
    </div>
  </header>

  <!-- PAGE HERO -->
  <section class="page-hero">
    <div class="container hero-content">
      <h2>Table S1.A: Biological Assembly of TFs & Complex Interfaces</h2>
      <p>
        Curated non-redundant structural dataset of {n_total} transcription factor complexes. Features biological assembly chain definitions, experimental methods, resolution, oligomeric states, and interface Buried Surface Area (BSA).
      </p>
      <div class="hero-stats-row">
        <div class="hero-stat-pill">
          <span class="number" id="stat-total">{n_total}</span>
          <span class="label">Total Complexes</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-dna">{n_dna}</span>
          <span class="label">DNA-Binding</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-rna">{n_rna}</span>
          <span class="label">RNA-Binding</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-dnarna">{n_dnarna}</span>
          <span class="label">DNA-RNA-Binding</span>
        </div>
      </div>
    </div>
  </section>

  <!-- CONTROLS BAR -->
  <main class="container controls-section">
    <div class="controls-card">
      <div class="controls-top">
        <div class="search-box-wrapper">
          <span class="search-icon">🔍</span>
          <input type="text" id="global-search" class="search-input" placeholder="Search PDB ID, protein name, organism, chains, oligomer..." autocomplete="off">
          <button id="clear-search" class="clear-search-btn" title="Clear search">✕</button>
        </div>
        <div class="action-buttons-group">
          <button id="btn-export-csv" class="btn-action">📥 Export CSV</button>
          <button id="btn-export-json" class="btn-action">📋 Export JSON</button>
          <button id="btn-copy-data" class="btn-action">📄 Copy</button>
        </div>
      </div>

      <!-- Category Filter Tabs -->
      <div class="filter-tabs">
        <span class="filter-tab-label">Filter Category:</span>
        <button class="filter-chip active" data-category="ALL">
          All Complexes <span class="chip-count" id="count-all">{n_total}</span>
        </button>
        <button class="filter-chip" data-category="DNA-binding">
          DNA-binding <span class="chip-count" id="count-dna">{n_dna}</span>
        </button>
        <button class="filter-chip" data-category="RNA-binding">
          RNA-binding <span class="chip-count" id="count-rna">{n_rna}</span>
        </button>
        <button class="filter-chip" data-category="DNA-RNA-binding">
          DNA-RNA-binding <span class="chip-count" id="count-dnarna">{n_dnarna}</span>
        </button>
      </div>
    </div>
  </main>

  <!-- DATA TABLE SECTION -->
  <section class="container table-section">
    <div class="table-card">
      <div class="table-header-info">
        <div class="results-count" id="results-count-text">
          Showing <strong id="showing-start">1</strong> to <strong id="showing-end">25</strong> of <strong id="showing-total">{n_total}</strong> complexes
        </div>
        <div class="page-size-selector">
          <label for="page-size">Entries per page:</label>
          <select id="page-size">
            <option value="15">15</option>
            <option value="25" selected>25</option>
            <option value="50">50</option>
            <option value="100">100</option>
            <option value="-1">All</option>
          </select>
        </div>
      </div>

      <div class="responsive-table-wrapper">
        <table class="data-table" id="dataset-table">
          <thead>
            <tr>
              <th data-key="pdb_id">PDB ID <span class="sort-indicator">↕</span></th>
              <th data-key="category">Category <span class="sort-indicator">↕</span></th>
              <th data-key="protein_name">Protein Name <span class="sort-indicator">↕</span></th>
              <th data-key="source_organism">Source Organism <span class="sort-indicator">↕</span></th>
              <th data-key="protein_length">Length <span class="sort-indicator">↕</span></th>
              <th data-key="bio_assembly">Bio Assembly <span class="sort-indicator">↕</span></th>
              <th data-key="exp_method">Method <span class="sort-indicator">↕</span></th>
              <th data-key="resolution" class="cell-num">Res (Å) <span class="sort-indicator">↕</span></th>
              <th data-key="oligomeric_state">Oligomeric State <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_complex" class="cell-num">BSA Complex (Å²) <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_protein" class="cell-num">BSA Prot <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_dna" class="cell-num">BSA DNA <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_rna" class="cell-num">BSA RNA <span class="sort-indicator">↕</span></th>
            </tr>
          </thead>
          <tbody id="table-body">
            <!-- Rendered by client JavaScript -->
          </tbody>
        </table>
      </div>

      <!-- PAGINATION -->
      <div class="pagination-container">
        <div class="pagination-info" id="pagination-summary">
          Page <strong id="current-page-num">1</strong> of <strong id="total-pages-num">1</strong>
        </div>
        <div class="pagination-buttons" id="pagination-btns">
          <!-- Page buttons rendered dynamically -->
        </div>
      </div>
    </div>
  </section>

  <!-- QUICK VIEW DETAIL MODAL -->
  <div class="modal-overlay" id="detail-modal">
    <div class="modal-card">
      <div class="modal-header">
        <h3 id="modal-title">PDB Details</h3>
        <button class="modal-close-btn" id="modal-close-btn">&times;</button>
      </div>
      <div class="modal-body" id="modal-content">
        <!-- Content filled dynamically -->
      </div>
      <div class="modal-footer">
        <a id="modal-rcsb-link" href="#" target="_blank" rel="noopener" class="btn-action btn-primary">Open in RCSB PDB ↗</a>
        <button class="btn-action" id="modal-close-action">Close</button>
      </div>
    </div>
  </div>

  <!-- SITE FOOTER -->
  <footer class="site-footer">
    <div class="container footer-bottom">
      <p>&copy; 2026 Computational Structural Biology Group, IIT Kharagpur. TF-NRD Database released under MIT License.</p>
    </div>
  </footer>

  <!-- CLIENT LOGIC & DATASET BUNDLE -->
  <script>
    const RAW_DATA = {json_data};

    let filteredData = [...RAW_DATA];
    let currentCategory = 'ALL';
    let searchQuery = '';
    let sortColumn = 'pdb_id';
    let sortDirection = 'asc';
    let currentPage = 1;
    let pageSize = 25;

    // DOM Elements
    const searchInput = document.getElementById('global-search');
    const clearSearchBtn = document.getElementById('clear-search');
    const tableBody = document.getElementById('table-body');
    const pageSizeSelect = document.getElementById('page-size');
    const paginationBtns = document.getElementById('pagination-btns');
    const showingStart = document.getElementById('showing-start');
    const showingEnd = document.getElementById('showing-end');
    const showingTotal = document.getElementById('showing-total');
    const currentPageNum = document.getElementById('current-page-num');
    const totalPagesNum = document.getElementById('total-pages-num');

    // Detail Modal Elements
    const modal = document.getElementById('detail-modal');
    const modalTitle = document.getElementById('modal-title');
    const modalContent = document.getElementById('modal-content');
    const modalRcsbLink = document.getElementById('modal-rcsb-link');
    const modalCloseBtn = document.getElementById('modal-close-btn');
    const modalCloseAction = document.getElementById('modal-close-action');

    function applyFilterAndSearch() {{
      filteredData = RAW_DATA.filter(row => {{
        const matchCategory = (currentCategory === 'ALL') || (row.category === currentCategory);
        if (!matchCategory) return false;
        if (!searchQuery) return true;
        const q = searchQuery.toLowerCase();
        return (
          row.pdb_id.toLowerCase().includes(q) ||
          row.protein_name.toLowerCase().includes(q) ||
          row.source_organism.toLowerCase().includes(q) ||
          row.bio_assembly.toLowerCase().includes(q) ||
          row.pdb_chain_id.toLowerCase().includes(q) ||
          row.oligomeric_state.toLowerCase().includes(q) ||
          row.exp_method.toLowerCase().includes(q)
        );
      }});

      // Sorting
      filteredData.sort((a, b) => {{
        let valA = a[sortColumn];
        let valB = b[sortColumn];

        let numA = parseFloat(valA);
        let numB = parseFloat(valB);

        if (!isNaN(numA) && !isNaN(numB)) {{
          return sortDirection === 'asc' ? numA - numB : numB - numA;
        }}
        valA = (valA || '').toString().toLowerCase();
        valB = (valB || '').toString().toLowerCase();
        if (valA < valB) return sortDirection === 'asc' ? -1 : 1;
        if (valA > valB) return sortDirection === 'asc' ? 1 : -1;
        return 0;
      }});

      currentPage = 1;
      renderTable();
    }}

    function renderTable() {{
      const total = filteredData.length;
      showingTotal.textContent = total;

      const size = pageSize === -1 ? total : pageSize;
      const totalPages = Math.max(1, Math.ceil(total / (size || 1)));
      currentPage = Math.min(currentPage, totalPages);

      const startIdx = (currentPage - 1) * size;
      const endIdx = pageSize === -1 ? total : Math.min(startIdx + size, total);

      showingStart.textContent = total === 0 ? 0 : startIdx + 1;
      showingEnd.textContent = endIdx;
      currentPageNum.textContent = currentPage;
      totalPagesNum.textContent = totalPages;

      const pageRows = filteredData.slice(startIdx, endIdx);

      if (pageRows.length === 0) {{
        tableBody.innerHTML = `<tr><td colspan="13" style="text-align: center; padding: 40px; color: var(--text-muted);">No complexes matching your query.</td></tr>`;
      }} else {{
        tableBody.innerHTML = pageRows.map(row => {{
          const catClass = row.category === 'DNA-binding' ? 'cat-dna' : (row.category === 'RNA-binding' ? 'cat-rna' : 'cat-dnarna');
          return `
            <tr onclick="openModal('${{row.pdb_id}}')">
              <td>
                <a href="https://www.rcsb.org/structure/${{row.pdb_id}}" target="_blank" rel="noopener" class="pdb-badge" onclick="event.stopPropagation()">
                  ${{row.pdb_id}} ↗
                </a>
              </td>
              <td><span class="badge-category ${{catClass}}">${{row.category}}</span></td>
              <td class="cell-protein-name" title="${{row.protein_name}}">${{row.protein_name}}</td>
              <td class="cell-organism">${{row.source_organism}}</td>
              <td class="cell-mono">${{row.protein_length}}</td>
              <td class="cell-mono">${{row.bio_assembly}}</td>
              <td><span class="method-pill">${{row.exp_method}}</span></td>
              <td class="cell-num">${{row.resolution}}</td>
              <td style="max-width: 180px; white-space: normal;">${{row.oligomeric_state}}</td>
              <td class="cell-num" style="font-weight: 600;">${{row.bsa_complex}}</td>
              <td class="cell-num">${{row.bsa_protein}}</td>
              <td class="cell-num">${{row.bsa_dna}}</td>
              <td class="cell-num">${{row.bsa_rna}}</td>
            </tr>
          `;
        }}).join('');
      }}

      renderPagination(totalPages);
    }}

    function renderPagination(totalPages) {{
      paginationBtns.innerHTML = '';
      if (totalPages <= 1) return;

      const addBtn = (text, page, isActive = false, isDisabled = false) => {{
        const btn = document.createElement('button');
        btn.className = `btn-page ${{isActive ? 'active' : ''}}`;
        btn.innerHTML = text;
        btn.disabled = isDisabled;
        btn.onclick = () => {{
          if (!isDisabled) {{
            currentPage = page;
            renderTable();
          }}
        }};
        paginationBtns.appendChild(btn);
      }};

      addBtn('«', 1, false, currentPage === 1);
      addBtn('‹', currentPage - 1, false, currentPage === 1);

      let start = Math.max(1, currentPage - 2);
      let end = Math.min(totalPages, currentPage + 2);

      if (start > 1) addBtn('1', 1);
      if (start > 2) {{
        const span = document.createElement('span');
        span.textContent = '...';
        span.style.padding = '0 6px';
        paginationBtns.appendChild(span);
      }}

      for (let i = start; i <= end; i++) {{
        addBtn(i, i, i === currentPage);
      }}

      if (end < totalPages - 1) {{
        const span = document.createElement('span');
        span.textContent = '...';
        span.style.padding = '0 6px';
        paginationBtns.appendChild(span);
      }}
      if (end < totalPages) addBtn(totalPages, totalPages);

      addBtn('›', currentPage + 1, false, currentPage === totalPages);
      addBtn('»', totalPages, false, currentPage === totalPages);
    }}

    function openModal(pdbId) {{
      const item = RAW_DATA.find(r => r.pdb_id === pdbId);
      if (!item) return;

      modalTitle.innerHTML = `Structure Details: <span style="color: var(--primary);">${{item.pdb_id}}</span>`;
      modalRcsbLink.href = `https://www.rcsb.org/structure/${{item.pdb_id}}`;

      modalContent.innerHTML = `
        <div class="modal-detail-grid">
          <div class="modal-detail-item">
            <div class="detail-label">Protein Name</div>
            <div class="detail-val">${{item.protein_name}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Source Organism</div>
            <div class="detail-val"><em>${{item.source_organism}}</em></div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Category</div>
            <div class="detail-val">${{item.category}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Experimental Method & Res</div>
            <div class="detail-val">${{item.exp_method}} (${{item.resolution}} Å)</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Biological Assembly (Prot:NA)</div>
            <div class="detail-val" style="font-family: var(--font-mono);">${{item.bio_assembly}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">PDB Chain IDs</div>
            <div class="detail-val" style="font-family: var(--font-mono);">${{item.pdb_chain_id}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Oligomeric State</div>
            <div class="detail-val">${{item.oligomeric_state}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Protein Length (Residues)</div>
            <div class="detail-val" style="font-family: var(--font-mono);">${{item.protein_length}}</div>
          </div>
        </div>

        <div style="margin-top: 10px; background: #eff6ff; padding: 14px; border-radius: var(--radius-md); border: 1px solid #bfdbfe;">
          <h4 style="font-size: 13px; font-weight: 700; color: #1e40af; margin-bottom: 8px;">Interface Buried Surface Area (BSA)</h4>
          <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(130px, 1fr)); gap: 10px;">
            <div><span style="font-size: 11px; color: #64748b;">Complex BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_complex}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">Protein BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_protein}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">DNA BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_dna}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">RNA BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_rna}} Å²</strong></div>
          </div>
        </div>
      `;

      modal.classList.add('open');
    }}

    function closeModal() {{
      modal.classList.remove('open');
    }}

    modalCloseBtn.onclick = closeModal;
    modalCloseAction.onclick = closeModal;
    modal.onclick = (e) => {{ if (e.target === modal) closeModal(); }};
    document.addEventListener('keydown', (e) => {{ if (e.key === 'Escape') closeModal(); }});

    // Event Listeners for Search & Filters
    searchInput.addEventListener('input', (e) => {{
      searchQuery = e.target.value.trim();
      clearSearchBtn.style.display = searchQuery ? 'block' : 'none';
      applyFilterAndSearch();
    }});

    clearSearchBtn.onclick = () => {{
      searchInput.value = '';
      searchQuery = '';
      clearSearchBtn.style.display = 'none';
      applyFilterAndSearch();
    }};

    document.querySelectorAll('.filter-chip').forEach(chip => {{
      chip.addEventListener('click', () => {{
        document.querySelectorAll('.filter-chip').forEach(c => c.classList.remove('active'));
        chip.classList.add('active');
        currentCategory = chip.getAttribute('data-category');
        applyFilterAndSearch();
      }});
    }});

    pageSizeSelect.addEventListener('change', (e) => {{
      pageSize = parseInt(e.target.value);
      currentPage = 1;
      renderTable();
    }});

    // Column Sorting
    document.querySelectorAll('table.data-table th[data-key]').forEach(th => {{
      th.addEventListener('click', () => {{
        const key = th.getAttribute('data-key');
        if (sortColumn === key) {{
          sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
        }} else {{
          sortColumn = key;
          sortDirection = 'asc';
        }}

        document.querySelectorAll('table.data-table th').forEach(h => {{
          h.classList.remove('sort-active');
          const ind = h.querySelector('.sort-indicator');
          if (ind) ind.textContent = '↕';
        }});

        th.classList.add('sort-active');
        const indicator = th.querySelector('.sort-indicator');
        if (indicator) indicator.textContent = sortDirection === 'asc' ? '▲' : '▼';

        applyFilterAndSearch();
      }});
    }});

    // Export Handlers
    document.getElementById('btn-export-csv').onclick = () => {{
      if (!filteredData.length) return;
      const headers = Object.keys(filteredData[0]);
      const csvRows = [
        headers.join(','),
        ...filteredData.map(row => headers.map(h => `"${{String(row[h] || '').replace(/"/g, '""')}}"`).join(','))
      ];
      const blob = new Blob([csvRows.join('\\n')], {{ type: 'text/csv' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `TFNRD_Table_S1A_${{currentCategory}}.csv`;
      a.click();
      URL.revokeObjectURL(url);
    }};

    document.getElementById('btn-export-json').onclick = () => {{
      if (!filteredData.length) return;
      const blob = new Blob([JSON.stringify(filteredData, null, 2)], {{ type: 'application/json' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `TFNRD_Table_S1A_${{currentCategory}}.json`;
      a.click();
      URL.revokeObjectURL(url);
    }};

    document.getElementById('btn-copy-data').onclick = () => {{
      if (!filteredData.length) return;
      const text = JSON.stringify(filteredData, null, 2);
      navigator.clipboard.writeText(text).then(() => {{
        alert('Copied ' + filteredData.length + ' records to clipboard in JSON format.');
      }});
    }};

    // Initial render
    applyFilterAndSearch();
  </script>
</body>
</html>
"""


def generate_table_s1b_html(records: list) -> str:
    """Generates the interactive modern webpage for Table S1.B (TFNRDv1.0_PNA.html)."""
    json_data = json.dumps(records, ensure_ascii=False)
    n_total = len(records)
    n_dna = sum(1 for r in records if r['category'] == 'DNA-binding')
    n_rna = sum(1 for r in records if r['category'] == 'RNA-binding')
    n_dnarna = sum(1 for r in records if r['category'] == 'DNA-RNA-binding')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="description" content="TF-NRD Table S1.B: Unique transcription factor-nucleic acid interface dataset, interacting chain pairs, and decomposed interface Buried Surface Area (BSA).">
  <title>TF-NRD TF-NA Interfaces (Table S1.B)</title>
  <link rel="stylesheet" href="css/style.css">
</head>
<body>

  <!-- SITE HEADER -->
  <header class="site-header">
    <div class="header-container">
      <a href="tfnrd.html" class="brand">
        <div class="brand-icon">TF</div>
        <div class="brand-text">
          <h1>TF-NRD <span class="version-tag">Table S1.B</span></h1>
          <p>Unique TF-NA Interface Dataset</p>
        </div>
      </a>
      <nav class="main-nav">
        <a href="tfnrd.html">Home</a>
        <a href="TFNRDv1.0.html">Dataset (Table S1.A)</a>
        <a href="TFNRDv1.0_PNA.html" class="active">TF-NA Interfaces (Table S1.B)</a>
        <a href="TFNRDv1.0_PP.html">Protein-Protein (Table S1.C)</a>
      </nav>
      <div class="header-actions">
        <a href="http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html" target="_blank" rel="noopener" class="btn-header-link">CSB Lab ↗</a>
      </div>
    </div>
  </header>

  <!-- PAGE HERO -->
  <section class="page-hero">
    <div class="container hero-content">
      <h2>Table S1.B: Unique Interface of Transcription Factors</h2>
      <p>
        Decomposed structural interface metrics for {n_total} TF-nucleic acid interfaces. Features interacting protein and nucleic acid chains, with Buried Surface Area (BSA) contributions across complex, protein, and DNA/RNA components.
      </p>
      <div class="hero-stats-row">
        <div class="hero-stat-pill">
          <span class="number" id="stat-total">{n_total}</span>
          <span class="label">Total Interfaces</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-dna">{n_dna}</span>
          <span class="label">DNA-Binding</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-rna">{n_rna}</span>
          <span class="label">RNA-Binding</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-dnarna">{n_dnarna}</span>
          <span class="label">DNA-RNA-Binding</span>
        </div>
      </div>
    </div>
  </section>

  <!-- CONTROLS BAR -->
  <main class="container controls-section">
    <div class="controls-card">
      <div class="controls-top">
        <div class="search-box-wrapper">
          <span class="search-icon">🔍</span>
          <input type="text" id="global-search" class="search-input" placeholder="Search PDB ID, protein name, organism, chains..." autocomplete="off">
          <button id="clear-search" class="clear-search-btn" title="Clear search">✕</button>
        </div>
        <div class="action-buttons-group">
          <button id="btn-export-csv" class="btn-action">📥 Export CSV</button>
          <button id="btn-export-json" class="btn-action">📋 Export JSON</button>
          <button id="btn-copy-data" class="btn-action">📄 Copy</button>
        </div>
      </div>

      <!-- Category Filter Tabs -->
      <div class="filter-tabs">
        <span class="filter-tab-label">Filter Category:</span>
        <button class="filter-chip active" data-category="ALL">
          All Interfaces <span class="chip-count" id="count-all">{n_total}</span>
        </button>
        <button class="filter-chip" data-category="DNA-binding">
          DNA-binding <span class="chip-count" id="count-dna">{n_dna}</span>
        </button>
        <button class="filter-chip" data-category="RNA-binding">
          RNA-binding <span class="chip-count" id="count-rna">{n_rna}</span>
        </button>
        <button class="filter-chip" data-category="DNA-RNA-binding">
          DNA-RNA-binding <span class="chip-count" id="count-dnarna">{n_dnarna}</span>
        </button>
      </div>
    </div>
  </main>

  <!-- DATA TABLE SECTION -->
  <section class="container table-section">
    <div class="table-card">
      <div class="table-header-info">
        <div class="results-count" id="results-count-text">
          Showing <strong id="showing-start">1</strong> to <strong id="showing-end">25</strong> of <strong id="showing-total">{n_total}</strong> interfaces
        </div>
        <div class="page-size-selector">
          <label for="page-size">Entries per page:</label>
          <select id="page-size">
            <option value="15">15</option>
            <option value="25" selected>25</option>
            <option value="50">50</option>
            <option value="100">100</option>
            <option value="-1">All</option>
          </select>
        </div>
      </div>

      <div class="responsive-table-wrapper">
        <table class="data-table" id="dataset-table">
          <thead>
            <tr>
              <th data-key="pdb_id">PDB ID <span class="sort-indicator">↕</span></th>
              <th data-key="category">Category <span class="sort-indicator">↕</span></th>
              <th data-key="protein_name">Protein Name <span class="sort-indicator">↕</span></th>
              <th data-key="source_organism">Source Organism <span class="sort-indicator">↕</span></th>
              <th data-key="protein_length">Length <span class="sort-indicator">↕</span></th>
              <th data-key="chain_protein">Chain (Prot) <span class="sort-indicator">↕</span></th>
              <th data-key="chain_dna">Chain (DNA) <span class="sort-indicator">↕</span></th>
              <th data-key="chain_rna">Chain (RNA) <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_complex" class="cell-num">BSA Complex (Å²) <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_protein" class="cell-num">BSA Prot <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_dna" class="cell-num">BSA DNA <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_rna" class="cell-num">BSA RNA <span class="sort-indicator">↕</span></th>
            </tr>
          </thead>
          <tbody id="table-body">
            <!-- Rendered by client JavaScript -->
          </tbody>
        </table>
      </div>

      <!-- PAGINATION -->
      <div class="pagination-container">
        <div class="pagination-info" id="pagination-summary">
          Page <strong id="current-page-num">1</strong> of <strong id="total-pages-num">1</strong>
        </div>
        <div class="pagination-buttons" id="pagination-btns">
          <!-- Page buttons rendered dynamically -->
        </div>
      </div>
    </div>
  </section>

  <!-- QUICK VIEW DETAIL MODAL -->
  <div class="modal-overlay" id="detail-modal">
    <div class="modal-card">
      <div class="modal-header">
        <h3 id="modal-title">Interface Details</h3>
        <button class="modal-close-btn" id="modal-close-btn">&times;</button>
      </div>
      <div class="modal-body" id="modal-content">
        <!-- Content filled dynamically -->
      </div>
      <div class="modal-footer">
        <a id="modal-rcsb-link" href="#" target="_blank" rel="noopener" class="btn-action btn-primary">Open in RCSB PDB ↗</a>
        <button class="btn-action" id="modal-close-action">Close</button>
      </div>
    </div>
  </div>

  <!-- SITE FOOTER -->
  <footer class="site-footer">
    <div class="container footer-bottom">
      <p>&copy; 2026 Computational Structural Biology Group, IIT Kharagpur. TF-NRD Database released under MIT License.</p>
    </div>
  </footer>

  <!-- CLIENT LOGIC & DATASET BUNDLE -->
  <script>
    const RAW_DATA = {json_data};

    let filteredData = [...RAW_DATA];
    let currentCategory = 'ALL';
    let searchQuery = '';
    let sortColumn = 'pdb_id';
    let sortDirection = 'asc';
    let currentPage = 1;
    let pageSize = 25;

    // DOM Elements
    const searchInput = document.getElementById('global-search');
    const clearSearchBtn = document.getElementById('clear-search');
    const tableBody = document.getElementById('table-body');
    const pageSizeSelect = document.getElementById('page-size');
    const paginationBtns = document.getElementById('pagination-btns');
    const showingStart = document.getElementById('showing-start');
    const showingEnd = document.getElementById('showing-end');
    const showingTotal = document.getElementById('showing-total');
    const currentPageNum = document.getElementById('current-page-num');
    const totalPagesNum = document.getElementById('total-pages-num');

    // Detail Modal Elements
    const modal = document.getElementById('detail-modal');
    const modalTitle = document.getElementById('modal-title');
    const modalContent = document.getElementById('modal-content');
    const modalRcsbLink = document.getElementById('modal-rcsb-link');
    const modalCloseBtn = document.getElementById('modal-close-btn');
    const modalCloseAction = document.getElementById('modal-close-action');

    function applyFilterAndSearch() {{
      filteredData = RAW_DATA.filter(row => {{
        const matchCategory = (currentCategory === 'ALL') || (row.category === currentCategory);
        if (!matchCategory) return false;
        if (!searchQuery) return true;
        const q = searchQuery.toLowerCase();
        return (
          row.pdb_id.toLowerCase().includes(q) ||
          row.protein_name.toLowerCase().includes(q) ||
          row.source_organism.toLowerCase().includes(q) ||
          row.chain_protein.toLowerCase().includes(q) ||
          row.chain_dna.toLowerCase().includes(q) ||
          row.chain_rna.toLowerCase().includes(q)
        );
      }});

      // Sorting
      filteredData.sort((a, b) => {{
        let valA = a[sortColumn];
        let valB = b[sortColumn];

        let numA = parseFloat(valA);
        let numB = parseFloat(valB);

        if (!isNaN(numA) && !isNaN(numB)) {{
          return sortDirection === 'asc' ? numA - numB : numB - numA;
        }}
        valA = (valA || '').toString().toLowerCase();
        valB = (valB || '').toString().toLowerCase();
        if (valA < valB) return sortDirection === 'asc' ? -1 : 1;
        if (valA > valB) return sortDirection === 'asc' ? 1 : -1;
        return 0;
      }});

      currentPage = 1;
      renderTable();
    }}

    function renderTable() {{
      const total = filteredData.length;
      showingTotal.textContent = total;

      const size = pageSize === -1 ? total : pageSize;
      const totalPages = Math.max(1, Math.ceil(total / (size || 1)));
      currentPage = Math.min(currentPage, totalPages);

      const startIdx = (currentPage - 1) * size;
      const endIdx = pageSize === -1 ? total : Math.min(startIdx + size, total);

      showingStart.textContent = total === 0 ? 0 : startIdx + 1;
      showingEnd.textContent = endIdx;
      currentPageNum.textContent = currentPage;
      totalPagesNum.textContent = totalPages;

      const pageRows = filteredData.slice(startIdx, endIdx);

      if (pageRows.length === 0) {{
        tableBody.innerHTML = `<tr><td colspan="12" style="text-align: center; padding: 40px; color: var(--text-muted);">No interfaces matching your query.</td></tr>`;
      }} else {{
        tableBody.innerHTML = pageRows.map(row => {{
          const catClass = row.category === 'DNA-binding' ? 'cat-dna' : (row.category === 'RNA-binding' ? 'cat-rna' : 'cat-dnarna');
          return `
            <tr onclick="openModal('${{row.pdb_id}}', '${{row.chain_protein}}', '${{row.chain_dna}}', '${{row.chain_rna}}')">
              <td>
                <a href="https://www.rcsb.org/structure/${{row.pdb_id}}" target="_blank" rel="noopener" class="pdb-badge" onclick="event.stopPropagation()">
                  ${{row.pdb_id}} ↗
                </a>
              </td>
              <td><span class="badge-category ${{catClass}}">${{row.category}}</span></td>
              <td class="cell-protein-name" title="${{row.protein_name}}">${{row.protein_name}}</td>
              <td class="cell-organism">${{row.source_organism}}</td>
              <td class="cell-mono">${{row.protein_length}}</td>
              <td class="cell-mono" style="font-weight: 700; color: #1e40af;">${{row.chain_protein}}</td>
              <td class="cell-mono" style="color: #047857;">${{row.chain_dna}}</td>
              <td class="cell-mono" style="color: #6d28d9;">${{row.chain_rna}}</td>
              <td class="cell-num" style="font-weight: 600;">${{row.bsa_complex}}</td>
              <td class="cell-num">${{row.bsa_protein}}</td>
              <td class="cell-num">${{row.bsa_dna}}</td>
              <td class="cell-num">${{row.bsa_rna}}</td>
            </tr>
          `;
        }}).join('');
      }}

      renderPagination(totalPages);
    }}

    function renderPagination(totalPages) {{
      paginationBtns.innerHTML = '';
      if (totalPages <= 1) return;

      const addBtn = (text, page, isActive = false, isDisabled = false) => {{
        const btn = document.createElement('button');
        btn.className = `btn-page ${{isActive ? 'active' : ''}}`;
        btn.innerHTML = text;
        btn.disabled = isDisabled;
        btn.onclick = () => {{
          if (!isDisabled) {{
            currentPage = page;
            renderTable();
          }}
        }};
        paginationBtns.appendChild(btn);
      }};

      addBtn('«', 1, false, currentPage === 1);
      addBtn('‹', currentPage - 1, false, currentPage === 1);

      let start = Math.max(1, currentPage - 2);
      let end = Math.min(totalPages, currentPage + 2);

      if (start > 1) addBtn('1', 1);
      if (start > 2) {{
        const span = document.createElement('span');
        span.textContent = '...';
        span.style.padding = '0 6px';
        paginationBtns.appendChild(span);
      }}

      for (let i = start; i <= end; i++) {{
        addBtn(i, i, i === currentPage);
      }}

      if (end < totalPages - 1) {{
        const span = document.createElement('span');
        span.textContent = '...';
        span.style.padding = '0 6px';
        paginationBtns.appendChild(span);
      }}
      if (end < totalPages) addBtn(totalPages, totalPages);

      addBtn('›', currentPage + 1, false, currentPage === totalPages);
      addBtn('»', totalPages, false, currentPage === totalPages);
    }}

    function openModal(pdbId, prot, dna, rna) {{
      const item = RAW_DATA.find(r => r.pdb_id === pdbId && r.chain_protein === prot && r.chain_dna === dna && r.chain_rna === rna);
      if (!item) return;

      modalTitle.innerHTML = `Interface Details: <span style="color: var(--primary);">${{item.pdb_id}}</span> (${{item.chain_protein}} : ${{item.chain_dna !== '-' ? item.chain_dna : item.chain_rna}})`;
      modalRcsbLink.href = `https://www.rcsb.org/structure/${{item.pdb_id}}`;

      modalContent.innerHTML = `
        <div class="modal-detail-grid">
          <div class="modal-detail-item">
            <div class="detail-label">Protein Name</div>
            <div class="detail-val">${{item.protein_name}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Source Organism</div>
            <div class="detail-val"><em>${{item.source_organism}}</em></div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Category</div>
            <div class="detail-val">${{item.category}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Protein Length</div>
            <div class="detail-val" style="font-family: var(--font-mono);">${{item.protein_length}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Protein Chains</div>
            <div class="detail-val" style="font-family: var(--font-mono); font-weight: 700;">${{item.chain_protein}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">DNA Chains / RNA Chains</div>
            <div class="detail-val" style="font-family: var(--font-mono);">${{item.chain_dna}} / ${{item.chain_rna}}</div>
          </div>
        </div>

        <div style="margin-top: 10px; background: #eff6ff; padding: 14px; border-radius: var(--radius-md); border: 1px solid #bfdbfe;">
          <h4 style="font-size: 13px; font-weight: 700; color: #1e40af; margin-bottom: 8px;">Decomposed Buried Surface Area (BSA)</h4>
          <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(130px, 1fr)); gap: 10px;">
            <div><span style="font-size: 11px; color: #64748b;">Complex BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_complex}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">Protein BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_protein}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">DNA BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_dna}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">RNA BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_rna}} Å²</strong></div>
          </div>
        </div>
      `;

      modal.classList.add('open');
    }}

    function closeModal() {{
      modal.classList.remove('open');
    }}

    modalCloseBtn.onclick = closeModal;
    modalCloseAction.onclick = closeModal;
    modal.onclick = (e) => {{ if (e.target === modal) closeModal(); }};
    document.addEventListener('keydown', (e) => {{ if (e.key === 'Escape') closeModal(); }});

    // Event Listeners for Search & Filters
    searchInput.addEventListener('input', (e) => {{
      searchQuery = e.target.value.trim();
      clearSearchBtn.style.display = searchQuery ? 'block' : 'none';
      applyFilterAndSearch();
    }});

    clearSearchBtn.onclick = () => {{
      searchInput.value = '';
      searchQuery = '';
      clearSearchBtn.style.display = 'none';
      applyFilterAndSearch();
    }};

    document.querySelectorAll('.filter-chip').forEach(chip => {{
      chip.addEventListener('click', () => {{
        document.querySelectorAll('.filter-chip').forEach(c => c.classList.remove('active'));
        chip.classList.add('active');
        currentCategory = chip.getAttribute('data-category');
        applyFilterAndSearch();
      }});
    }});

    pageSizeSelect.addEventListener('change', (e) => {{
      pageSize = parseInt(e.target.value);
      currentPage = 1;
      renderTable();
    }});

    // Column Sorting
    document.querySelectorAll('table.data-table th[data-key]').forEach(th => {{
      th.addEventListener('click', () => {{
        const key = th.getAttribute('data-key');
        if (sortColumn === key) {{
          sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
        }} else {{
          sortColumn = key;
          sortDirection = 'asc';
        }}

        document.querySelectorAll('table.data-table th').forEach(h => {{
          h.classList.remove('sort-active');
          const ind = h.querySelector('.sort-indicator');
          if (ind) ind.textContent = '↕';
        }});

        th.classList.add('sort-active');
        const indicator = th.querySelector('.sort-indicator');
        if (indicator) indicator.textContent = sortDirection === 'asc' ? '▲' : '▼';

        applyFilterAndSearch();
      }});
    }});

    // Export Handlers
    document.getElementById('btn-export-csv').onclick = () => {{
      if (!filteredData.length) return;
      const headers = Object.keys(filteredData[0]);
      const csvRows = [
        headers.join(','),
        ...filteredData.map(row => headers.map(h => `"${{String(row[h] || '').replace(/"/g, '""')}}"`).join(','))
      ];
      const blob = new Blob([csvRows.join('\\n')], {{ type: 'text/csv' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `TFNRD_Table_S1B_${{currentCategory}}.csv`;
      a.click();
      URL.revokeObjectURL(url);
    }};

    document.getElementById('btn-export-json').onclick = () => {{
      if (!filteredData.length) return;
      const blob = new Blob([JSON.stringify(filteredData, null, 2)], {{ type: 'application/json' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `TFNRD_Table_S1B_${{currentCategory}}.json`;
      a.click();
      URL.revokeObjectURL(url);
    }};

    document.getElementById('btn-copy-data').onclick = () => {{
      if (!filteredData.length) return;
      const text = JSON.stringify(filteredData, null, 2);
      navigator.clipboard.writeText(text).then(() => {{
        alert('Copied ' + filteredData.length + ' records to clipboard in JSON format.');
      }});
    }};

    // Initial render
    applyFilterAndSearch();
  </script>
</body>
</html>
"""


def generate_table_s1c_html(records: list) -> str:
    """Generates the interactive modern webpage for Table S1.C (TFNRDv1.0_PP.html)."""
    json_data = json.dumps(records, ensure_ascii=False)
    n_total = len(records)
    n_dna = sum(1 for r in records if r['category'] == 'DNA-binding')
    n_rna = sum(1 for r in records if r['category'] == 'RNA-binding')
    n_dnarna = sum(1 for r in records if r['category'] == 'DNA-RNA-binding')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="description" content="TF-NRD Table S1.C: Protein-Protein contact interfaces within transcription factor complexes, subunit pairings, and Buried Surface Area (BSA).">
  <title>TF-NRD Protein-Protein Interfaces (Table S1.C)</title>
  <link rel="stylesheet" href="css/style.css">
</head>
<body>

  <!-- SITE HEADER -->
  <header class="site-header">
    <div class="header-container">
      <a href="tfnrd.html" class="brand">
        <div class="brand-icon">TF</div>
        <div class="brand-text">
          <h1>TF-NRD <span class="version-tag">Table S1.C</span></h1>
          <p>Protein-Protein Interface Dataset</p>
        </div>
      </a>
      <nav class="main-nav">
        <a href="tfnrd.html">Home</a>
        <a href="TFNRDv1.0.html">Dataset (Table S1.A)</a>
        <a href="TFNRDv1.0_PNA.html">TF-NA Interfaces (Table S1.B)</a>
        <a href="TFNRDv1.0_PP.html" class="active">Protein-Protein (Table S1.C)</a>
      </nav>
      <div class="header-actions">
        <a href="http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html" target="_blank" rel="noopener" class="btn-header-link">CSB Lab ↗</a>
      </div>
    </div>
  </header>

  <!-- PAGE HERO -->
  <section class="page-hero">
    <div class="container hero-content">
      <h2>Table S1.C: Protein-Protein Interfaces of Transcription Factors</h2>
      <p>
        Detailed contact interfaces for {n_total} protein-protein interactions within TF assemblies. Features Chain 1 and Chain 2 pairings, with individual subunit and total complex Buried Surface Area (BSA).
      </p>
      <div class="hero-stats-row">
        <div class="hero-stat-pill">
          <span class="number" id="stat-total">{n_total}</span>
          <span class="label">Total PP Interfaces</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-dna">{n_dna}</span>
          <span class="label">DNA-Binding</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-rna">{n_rna}</span>
          <span class="label">RNA-Binding</span>
        </div>
        <div class="hero-stat-pill">
          <span class="number" id="stat-dnarna">{n_dnarna}</span>
          <span class="label">DNA-RNA-Binding</span>
        </div>
      </div>
    </div>
  </section>

  <!-- CONTROLS BAR -->
  <main class="container controls-section">
    <div class="controls-card">
      <div class="controls-top">
        <div class="search-box-wrapper">
          <span class="search-icon">🔍</span>
          <input type="text" id="global-search" class="search-input" placeholder="Search PDB ID, protein name, organism, chains..." autocomplete="off">
          <button id="clear-search" class="clear-search-btn" title="Clear search">✕</button>
        </div>
        <div class="action-buttons-group">
          <button id="btn-export-csv" class="btn-action">📥 Export CSV</button>
          <button id="btn-export-json" class="btn-action">📋 Export JSON</button>
          <button id="btn-copy-data" class="btn-action">📄 Copy</button>
        </div>
      </div>

      <!-- Category Filter Tabs -->
      <div class="filter-tabs">
        <span class="filter-tab-label">Filter Category:</span>
        <button class="filter-chip active" data-category="ALL">
          All Interfaces <span class="chip-count" id="count-all">{n_total}</span>
        </button>
        <button class="filter-chip" data-category="DNA-binding">
          DNA-binding <span class="chip-count" id="count-dna">{n_dna}</span>
        </button>
        <button class="filter-chip" data-category="RNA-binding">
          RNA-binding <span class="chip-count" id="count-rna">{n_rna}</span>
        </button>
        <button class="filter-chip" data-category="DNA-RNA-binding">
          DNA-RNA-binding <span class="chip-count" id="count-dnarna">{n_dnarna}</span>
        </button>
      </div>
    </div>
  </main>

  <!-- DATA TABLE SECTION -->
  <section class="container table-section">
    <div class="table-card">
      <div class="table-header-info">
        <div class="results-count" id="results-count-text">
          Showing <strong id="showing-start">1</strong> to <strong id="showing-end">25</strong> of <strong id="showing-total">{n_total}</strong> interfaces
        </div>
        <div class="page-size-selector">
          <label for="page-size">Entries per page:</label>
          <select id="page-size">
            <option value="15">15</option>
            <option value="25" selected>25</option>
            <option value="50">50</option>
            <option value="100">100</option>
            <option value="-1">All</option>
          </select>
        </div>
      </div>

      <div class="responsive-table-wrapper">
        <table class="data-table" id="dataset-table">
          <thead>
            <tr>
              <th data-key="pdb_id">PDB ID <span class="sort-indicator">↕</span></th>
              <th data-key="category">Category <span class="sort-indicator">↕</span></th>
              <th data-key="protein_name">Protein Name <span class="sort-indicator">↕</span></th>
              <th data-key="source_organism">Source Organism <span class="sort-indicator">↕</span></th>
              <th data-key="chain1">Chain 1 <span class="sort-indicator">↕</span></th>
              <th data-key="chain2">Chain 2 <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_complex" class="cell-num">BSA Complex (Å²) <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_protein1" class="cell-num">BSA Protein 1 <span class="sort-indicator">↕</span></th>
              <th data-key="bsa_protein2" class="cell-num">BSA Protein 2 <span class="sort-indicator">↕</span></th>
            </tr>
          </thead>
          <tbody id="table-body">
            <!-- Rendered by client JavaScript -->
          </tbody>
        </table>
      </div>

      <!-- PAGINATION -->
      <div class="pagination-container">
        <div class="pagination-info" id="pagination-summary">
          Page <strong id="current-page-num">1</strong> of <strong id="total-pages-num">1</strong>
        </div>
        <div class="pagination-buttons" id="pagination-btns">
          <!-- Page buttons rendered dynamically -->
        </div>
      </div>
    </div>
  </section>

  <!-- QUICK VIEW DETAIL MODAL -->
  <div class="modal-overlay" id="detail-modal">
    <div class="modal-card">
      <div class="modal-header">
        <h3 id="modal-title">Interface Details</h3>
        <button class="modal-close-btn" id="modal-close-btn">&times;</button>
      </div>
      <div class="modal-body" id="modal-content">
        <!-- Content filled dynamically -->
      </div>
      <div class="modal-footer">
        <a id="modal-rcsb-link" href="#" target="_blank" rel="noopener" class="btn-action btn-primary">Open in RCSB PDB ↗</a>
        <button class="btn-action" id="modal-close-action">Close</button>
      </div>
    </div>
  </div>

  <!-- SITE FOOTER -->
  <footer class="site-footer">
    <div class="container footer-bottom">
      <p>&copy; 2026 Computational Structural Biology Group, IIT Kharagpur. TF-NRD Database released under MIT License.</p>
    </div>
  </footer>

  <!-- CLIENT LOGIC & DATASET BUNDLE -->
  <script>
    const RAW_DATA = {json_data};

    let filteredData = [...RAW_DATA];
    let currentCategory = 'ALL';
    let searchQuery = '';
    let sortColumn = 'pdb_id';
    let sortDirection = 'asc';
    let currentPage = 1;
    let pageSize = 25;

    // DOM Elements
    const searchInput = document.getElementById('global-search');
    const clearSearchBtn = document.getElementById('clear-search');
    const tableBody = document.getElementById('table-body');
    const pageSizeSelect = document.getElementById('page-size');
    const paginationBtns = document.getElementById('pagination-btns');
    const showingStart = document.getElementById('showing-start');
    const showingEnd = document.getElementById('showing-end');
    const showingTotal = document.getElementById('showing-total');
    const currentPageNum = document.getElementById('current-page-num');
    const totalPagesNum = document.getElementById('total-pages-num');

    // Detail Modal Elements
    const modal = document.getElementById('detail-modal');
    const modalTitle = document.getElementById('modal-title');
    const modalContent = document.getElementById('modal-content');
    const modalRcsbLink = document.getElementById('modal-rcsb-link');
    const modalCloseBtn = document.getElementById('modal-close-btn');
    const modalCloseAction = document.getElementById('modal-close-action');

    function applyFilterAndSearch() {{
      filteredData = RAW_DATA.filter(row => {{
        const matchCategory = (currentCategory === 'ALL') || (row.category === currentCategory);
        if (!matchCategory) return false;
        if (!searchQuery) return true;
        const q = searchQuery.toLowerCase();
        return (
          row.pdb_id.toLowerCase().includes(q) ||
          row.protein_name.toLowerCase().includes(q) ||
          row.source_organism.toLowerCase().includes(q) ||
          row.chain1.toLowerCase().includes(q) ||
          row.chain2.toLowerCase().includes(q)
        );
      }});

      // Sorting
      filteredData.sort((a, b) => {{
        let valA = a[sortColumn];
        let valB = b[sortColumn];

        let numA = parseFloat(valA);
        let numB = parseFloat(valB);

        if (!isNaN(numA) && !isNaN(numB)) {{
          return sortDirection === 'asc' ? numA - numB : numB - numA;
        }}
        valA = (valA || '').toString().toLowerCase();
        valB = (valB || '').toString().toLowerCase();
        if (valA < valB) return sortDirection === 'asc' ? -1 : 1;
        if (valA > valB) return sortDirection === 'asc' ? 1 : -1;
        return 0;
      }});

      currentPage = 1;
      renderTable();
    }}

    function renderTable() {{
      const total = filteredData.length;
      showingTotal.textContent = total;

      const size = pageSize === -1 ? total : pageSize;
      const totalPages = Math.max(1, Math.ceil(total / (size || 1)));
      currentPage = Math.min(currentPage, totalPages);

      const startIdx = (currentPage - 1) * size;
      const endIdx = pageSize === -1 ? total : Math.min(startIdx + size, total);

      showingStart.textContent = total === 0 ? 0 : startIdx + 1;
      showingEnd.textContent = endIdx;
      currentPageNum.textContent = currentPage;
      totalPagesNum.textContent = totalPages;

      const pageRows = filteredData.slice(startIdx, endIdx);

      if (pageRows.length === 0) {{
        tableBody.innerHTML = `<tr><td colspan="9" style="text-align: center; padding: 40px; color: var(--text-muted);">No protein-protein interfaces matching your query.</td></tr>`;
      }} else {{
        tableBody.innerHTML = pageRows.map(row => {{
          const catClass = row.category === 'DNA-binding' ? 'cat-dna' : (row.category === 'RNA-binding' ? 'cat-rna' : 'cat-dnarna');
          return `
            <tr onclick="openModal('${{row.pdb_id}}', '${{row.chain1}}', '${{row.chain2}}')">
              <td>
                <a href="https://www.rcsb.org/structure/${{row.pdb_id}}" target="_blank" rel="noopener" class="pdb-badge" onclick="event.stopPropagation()">
                  ${{row.pdb_id}} ↗
                </a>
              </td>
              <td><span class="badge-category ${{catClass}}">${{row.category}}</span></td>
              <td class="cell-protein-name" title="${{row.protein_name}}">${{row.protein_name}}</td>
              <td class="cell-organism">${{row.source_organism}}</td>
              <td class="cell-mono" style="font-weight: 700; color: #1e40af;">${{row.chain1}}</td>
              <td class="cell-mono" style="font-weight: 700; color: #7c3aed;">${{row.chain2}}</td>
              <td class="cell-num" style="font-weight: 600;">${{row.bsa_complex}}</td>
              <td class="cell-num">${{row.bsa_protein1}}</td>
              <td class="cell-num">${{row.bsa_protein2}}</td>
            </tr>
          `;
        }}).join('');
      }}

      renderPagination(totalPages);
    }}

    function renderPagination(totalPages) {{
      paginationBtns.innerHTML = '';
      if (totalPages <= 1) return;

      const addBtn = (text, page, isActive = false, isDisabled = false) => {{
        const btn = document.createElement('button');
        btn.className = `btn-page ${{isActive ? 'active' : ''}}`;
        btn.innerHTML = text;
        btn.disabled = isDisabled;
        btn.onclick = () => {{
          if (!isDisabled) {{
            currentPage = page;
            renderTable();
          }}
        }};
        paginationBtns.appendChild(btn);
      }};

      addBtn('«', 1, false, currentPage === 1);
      addBtn('‹', currentPage - 1, false, currentPage === 1);

      let start = Math.max(1, currentPage - 2);
      let end = Math.min(totalPages, currentPage + 2);

      if (start > 1) addBtn('1', 1);
      if (start > 2) {{
        const span = document.createElement('span');
        span.textContent = '...';
        span.style.padding = '0 6px';
        paginationBtns.appendChild(span);
      }}

      for (let i = start; i <= end; i++) {{
        addBtn(i, i, i === currentPage);
      }}

      if (end < totalPages - 1) {{
        const span = document.createElement('span');
        span.textContent = '...';
        span.style.padding = '0 6px';
        paginationBtns.appendChild(span);
      }}
      if (end < totalPages) addBtn(totalPages, totalPages);

      addBtn('›', currentPage + 1, false, currentPage === totalPages);
      addBtn('»', totalPages, false, currentPage === totalPages);
    }}

    function openModal(pdbId, c1, c2) {{
      const item = RAW_DATA.find(r => r.pdb_id === pdbId && r.chain1 === c1 && r.chain2 === c2);
      if (!item) return;

      modalTitle.innerHTML = `Protein-Protein Interface: <span style="color: var(--primary);">${{item.pdb_id}}</span> (${{item.chain1}} &ndash; ${{item.chain2}})`;
      modalRcsbLink.href = `https://www.rcsb.org/structure/${{item.pdb_id}}`;

      modalContent.innerHTML = `
        <div class="modal-detail-grid">
          <div class="modal-detail-item">
            <div class="detail-label">Protein Name</div>
            <div class="detail-val">${{item.protein_name}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Source Organism</div>
            <div class="detail-val"><em>${{item.source_organism}}</em></div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Category</div>
            <div class="detail-val">${{item.category}}</div>
          </div>
          <div class="modal-detail-item">
            <div class="detail-label">Interacting Subunit Pair</div>
            <div class="detail-val" style="font-family: var(--font-mono); font-weight: 700;">Chain ${{item.chain1}} &harr; Chain ${{item.chain2}}</div>
          </div>
        </div>

        <div style="margin-top: 10px; background: #eff6ff; padding: 14px; border-radius: var(--radius-md); border: 1px solid #bfdbfe;">
          <h4 style="font-size: 13px; font-weight: 700; color: #1e40af; margin-bottom: 8px;">Subunit Buried Surface Area (BSA)</h4>
          <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 10px;">
            <div><span style="font-size: 11px; color: #64748b;">Total Complex BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_complex}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">Chain ${{item.chain1}} BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_protein1}} Å²</strong></div>
            <div><span style="font-size: 11px; color: #64748b;">Chain ${{item.chain2}} BSA:</span> <strong style="display: block; font-size: 16px; color: #0f172a;">${{item.bsa_protein2}} Å²</strong></div>
          </div>
        </div>
      `;

      modal.classList.add('open');
    }}

    function closeModal() {{
      modal.classList.remove('open');
    }}

    modalCloseBtn.onclick = closeModal;
    modalCloseAction.onclick = closeModal;
    modal.onclick = (e) => {{ if (e.target === modal) closeModal(); }};
    document.addEventListener('keydown', (e) => {{ if (e.key === 'Escape') closeModal(); }});

    // Event Listeners for Search & Filters
    searchInput.addEventListener('input', (e) => {{
      searchQuery = e.target.value.trim();
      clearSearchBtn.style.display = searchQuery ? 'block' : 'none';
      applyFilterAndSearch();
    }});

    clearSearchBtn.onclick = () => {{
      searchInput.value = '';
      searchQuery = '';
      clearSearchBtn.style.display = 'none';
      applyFilterAndSearch();
    }};

    document.querySelectorAll('.filter-chip').forEach(chip => {{
      chip.addEventListener('click', () => {{
        document.querySelectorAll('.filter-chip').forEach(c => c.classList.remove('active'));
        chip.classList.add('active');
        currentCategory = chip.getAttribute('data-category');
        applyFilterAndSearch();
      }});
    }});

    pageSizeSelect.addEventListener('change', (e) => {{
      pageSize = parseInt(e.target.value);
      currentPage = 1;
      renderTable();
    }});

    // Column Sorting
    document.querySelectorAll('table.data-table th[data-key]').forEach(th => {{
      th.addEventListener('click', () => {{
        const key = th.getAttribute('data-key');
        if (sortColumn === key) {{
          sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
        }} else {{
          sortColumn = key;
          sortDirection = 'asc';
        }}

        document.querySelectorAll('table.data-table th').forEach(h => {{
          h.classList.remove('sort-active');
          const ind = h.querySelector('.sort-indicator');
          if (ind) ind.textContent = '↕';
        }});

        th.classList.add('sort-active');
        const indicator = th.querySelector('.sort-indicator');
        if (indicator) indicator.textContent = sortDirection === 'asc' ? '▲' : '▼';

        applyFilterAndSearch();
      }});
    }});

    // Export Handlers
    document.getElementById('btn-export-csv').onclick = () => {{
      if (!filteredData.length) return;
      const headers = Object.keys(filteredData[0]);
      const csvRows = [
        headers.join(','),
        ...filteredData.map(row => headers.map(h => `"${{String(row[h] || '').replace(/"/g, '""')}}"`).join(','))
      ];
      const blob = new Blob([csvRows.join('\\n')], {{ type: 'text/csv' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `TFNRD_Table_S1C_${{currentCategory}}.csv`;
      a.click();
      URL.revokeObjectURL(url);
    }};

    document.getElementById('btn-export-json').onclick = () => {{
      if (!filteredData.length) return;
      const blob = new Blob([JSON.stringify(filteredData, null, 2)], {{ type: 'application/json' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `TFNRD_Table_S1C_${{currentCategory}}.json`;
      a.click();
      URL.revokeObjectURL(url);
    }};

    document.getElementById('btn-copy-data').onclick = () => {{
      if (!filteredData.length) return;
      const text = JSON.stringify(filteredData, null, 2);
      navigator.clipboard.writeText(text).then(() => {{
        alert('Copied ' + filteredData.length + ' records to clipboard in JSON format.');
      }});
    }};

    // Initial render
    applyFilterAndSearch();
  </script>
</body>
</html>
"""


def build_webpages():
    """Builds all HTML pages and writes them to target directories."""
    if not EXCEL_PATH.exists():
        logger.error(f"Supplementary Excel file not found: {EXCEL_PATH}")
        sys.exit(1)

    logger.info(f"Loading master supplementary file: {EXCEL_PATH}")
    xls = pd.ExcelFile(EXCEL_PATH)

    s1a_records = parse_table_s1a(xls)
    s1b_records = parse_table_s1b(xls)
    s1c_records = parse_table_s1c(xls, s1a_records)

    logger.info("Generating modern HTML and CSS content...")
    css_text = get_css_content()
    tfnrd_html = generate_tfnrd_home_html(s1a_records, s1b_records, s1c_records)
    s1a_html = generate_table_s1a_html(s1a_records)
    s1b_html = generate_table_s1b_html(s1b_records)
    s1c_html = generate_table_s1c_html(s1c_records)

    pages = {
        "tfnrd.html": tfnrd_html,
        "TFNRDv1.0.html": s1a_html,
        "TFNRDv1.0_PNA.html": s1b_html,
        "TFNRDv1.0_PP.html": s1c_html,
        "css/style.css": css_text
    }

    for target_dir in TARGET_DIRS:
        target_dir.mkdir(parents=True, exist_ok=True)
        (target_dir / "css").mkdir(parents=True, exist_ok=True)

        for filename, content in pages.items():
            out_file = target_dir / filename
            out_file.write_text(content, encoding='utf-8')
            logger.info(f"Saved: {out_file} ({len(content):,} bytes)")

    logger.info("Successfully generated all TF-NRD webpages in TFNRDv1.0 and root directory!")


if __name__ == "__main__":
    build_webpages()
