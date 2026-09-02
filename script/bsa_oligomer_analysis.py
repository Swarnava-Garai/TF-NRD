#!/usr/bin/env python3
"""
bsa_oligomer_analysis.py
========================
Comprehensive statistical analysis and visualization suite for Buried Surface Area (BSA)
parameters categorized across 3 oligomeric classes:
  1. Monomer
  2. Homodimer
  3. Heterodimer

Generates Figure 5 for the TF-NRD manuscript:
  - Focused BSA Complex violin plots with Mean ± SD & N sample badges
  - Combined publication-ready 2-panel Figure 5 (Panel A: Protein-NA, Panel B: Protein-Protein)

Supports both:
  - Protein-Nucleic Acid complexes ('TF_nucleic_acid_whole_with_oligostate.xlsx')
  - Protein-Protein complexes ('TF_protein_protein_with_oligostate.xlsx')

Results and figures are saved into dedicated separate subdirectories under results/:
  - results/Interface/protein_nucleic_acid/ (and Figures/)
  - results/Interface/protein_protein/ (and Figures/)
  - results/Figures/
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

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

DEFAULT_INPUT_NA = PROJECT_ROOT / "input_data" / "bsa" / "TF_nucleic_acid_whole_with_oligostate.xlsx"
DEFAULT_INPUT_PP = PROJECT_ROOT / "input_data" / "bsa" / "TF_protein_protein_with_oligostate.xlsx"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "Interface"


def resolve_path(path_input: Any) -> Path:
    """
    Intelligently resolves relative/absolute file paths whether invoked from repo root,
    script dir, or workspace root.
    """
    if isinstance(path_input, Path):
        p = path_input
    else:
        p = Path(str(path_input).strip())

    if p.is_absolute() and p.exists():
        return p

    candidates = [
        p.resolve(),
        PROJECT_ROOT / p,
        PROJECT_ROOT / str(p).replace("TF-NRD/", "", 1),
        Path.cwd() / p,
        Path.cwd() / str(p).replace("TF-NRD/", "", 1),
    ]

    for cand in candidates:
        if cand.exists():
            return cand.resolve()

    return p.resolve()


class BSAOligomerAnalyzer:
    """
    Comprehensive Statistical & Visual Analyzer for BSA distributions across oligomeric classes.
    """

    CLASS_ORDER = ['Monomer', 'Homodimer', 'Heterodimer']

    CLASS_PALETTE = {
        'Monomer': '#2b5c8f',             # Deep Steel Blue
        'Homodimer': '#d95f02',           # Burnt Orange
        'Heterodimer': '#7570b3'          # Royal Purple
    }

    def __init__(self, excel_path: Path, dataset_name: str = "protein_nucleic_acid", base_output_dir: Path = DEFAULT_OUTPUT_DIR):
        self.excel_path = resolve_path(excel_path)
        self.dataset_name = dataset_name
        self.output_dir = Path(base_output_dir).resolve() / dataset_name
        self.figures_dir = self.output_dir / "Figures"

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)

        self.df = None
        self.monomer_df = None
        self.homodimer_df = None
        self.heterodimer_df = None
        self.class_dfs = {}
        self.bsa_columns = []

    def load_and_preprocess(self) -> pd.DataFrame:
        """
        Loads the dataset, standardizes oligomeric state column,
        categorizes into Monomer, Homodimer, Heterodimer, and filters higher order oligomers.
        """
        if not self.excel_path.exists():
            raise FileNotFoundError(f"Input file not found at: {self.excel_path}")

        logger.info(f"Loading {self.dataset_name} dataset from: {self.excel_path}")
        ext = self.excel_path.suffix.lower()
        if ext in ('.xlsx', '.xls', '.xlsm'):
            self.df = pd.read_excel(self.excel_path)
        elif ext == '.csv':
            self.df = pd.read_csv(self.excel_path)
        elif ext == '.tsv':
            self.df = pd.read_csv(self.excel_path, sep='\t')
        elif ext == '.json':
            self.df = pd.read_json(self.excel_path)
        else:
            raise ValueError(f"Unsupported format: {ext}")

        logger.info(f"Initial raw record count: {len(self.df)}")

        # Find oligomeric state column
        oligo_col = None
        for col in self.df.columns:
            if col.strip().lower() in ['oligomeric state', 'oligomer_state', 'oligomeric_state', 'oligostate']:
                oligo_col = col
                break

        if not oligo_col:
            raise ValueError(f"Could not find Oligomeric State column in {self.excel_path}. Columns: {list(self.df.columns)}")

        # Clean string values
        self.df['Oligomeric_State_Clean'] = self.df[oligo_col].astype(str).str.strip().str.capitalize()

        # Map to 3 classes
        def classify_oligo(val: str) -> Optional[str]:
            val_lower = val.lower()
            if 'monomer' in val_lower:
                return 'Monomer'
            elif 'homodimer' in val_lower:
                return 'Homodimer'
            elif 'heterodimer' in val_lower:
                return 'Heterodimer'
            else:
                return None  # Higher order oligomers excluded

        self.df['Oligomer_Class'] = self.df['Oligomeric_State_Clean'].apply(classify_oligo)

        # Log higher order exclusion count
        excluded_count = self.df['Oligomer_Class'].isna().sum()
        logger.info(f"Excluding {excluded_count} higher order oligomer records (analyzing Monomer, Homodimer, Heterodimer only).")

        self.df = self.df[self.df['Oligomer_Class'].notna()].copy()
        logger.info(f"Filtered dataset record count: {len(self.df)}")

        # Identify BSA numeric columns
        self.bsa_columns = []
        for col in self.df.columns:
            col_lower = col.lower()
            if 'bsa' in col_lower or 'buried' in col_lower or 'area' in col_lower:
                if np.issubdtype(self.df[col].dtype, np.number) or pd.to_numeric(self.df[col], errors='coerce').notna().any():
                    self.df[col] = pd.to_numeric(self.df[col], errors='coerce')
                    self.bsa_columns.append(col)

        # Deduplicate and sort BSA columns
        self.bsa_columns = list(dict.fromkeys(self.bsa_columns))
        logger.info(f"Identified BSA numerical columns for analysis: {self.bsa_columns}")

        # Partition by class
        for c in self.CLASS_ORDER:
            sub = self.df[self.df['Oligomer_Class'] == c].copy()
            self.class_dfs[c] = sub
            logger.info(f"  - {c}: {len(sub)} entries")

        self.monomer_df = self.class_dfs['Monomer']
        self.homodimer_df = self.class_dfs['Homodimer']
        self.heterodimer_df = self.class_dfs['Heterodimer']

        return self.df

    @staticmethod
    def analyze_column_and_outliers(df: pd.DataFrame, column: str, iqr_factor: float = 1.5, verbose: bool = False) -> Tuple[Dict[str, Any], pd.DataFrame]:
        """
        Computes descriptive statistics (N, Mean, SD, CV, Median, Min, Max, IQR) and flags IQR-based outliers.
        """
        if column not in df.columns or df.empty:
            return {
                'total_count': 0, 'mean': np.nan, 'std': np.nan, 'cv_pct': np.nan,
                'median': np.nan, 'min': np.nan, 'max': np.nan, 'q1': np.nan, 'q3': np.nan,
                'iqr': np.nan, 'lower_bound': np.nan, 'upper_bound': np.nan,
                'outlier_count': 0, 'outlier_pct': 0.0
            }, pd.DataFrame()

        series = df[column].dropna()
        n = len(series)

        if n == 0:
            return {
                'total_count': 0, 'mean': np.nan, 'std': np.nan, 'cv_pct': np.nan,
                'median': np.nan, 'min': np.nan, 'max': np.nan, 'q1': np.nan, 'q3': np.nan,
                'iqr': np.nan, 'lower_bound': np.nan, 'upper_bound': np.nan,
                'outlier_count': 0, 'outlier_pct': 0.0
            }, pd.DataFrame()

        mean_val = series.mean()
        std_val = series.std()
        cv_val = (std_val / mean_val * 100) if mean_val != 0 else np.nan
        median_val = series.median()
        min_val = series.min()
        max_val = series.max()

        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1

        lower_bound = q1 - iqr_factor * iqr
        upper_bound = q3 + iqr_factor * iqr

        outlier_mask = (df[column] < lower_bound) | (df[column] > upper_bound)
        outliers = df[outlier_mask].copy()

        stats_dict = {
            'total_count': n,
            'mean': mean_val,
            'std': std_val,
            'cv_pct': cv_val,
            'median': median_val,
            'min': min_val,
            'max': max_val,
            'q1': q1,
            'q3': q3,
            'iqr': iqr,
            'lower_bound': lower_bound,
            'upper_bound': upper_bound,
            'outlier_count': len(outliers),
            'outlier_pct': (len(outliers) / n * 100) if n > 0 else 0.0
        }

        if verbose:
            logger.info(f"--- Analysis for '{column}' ---")
            logger.info(f"Count: {n} | Mean: {mean_val:.2f} | Std: {std_val:.2f} | CV (%): {cv_val:.2f}%")
            logger.info(f"Median: {median_val:.2f} | Q1: {q1:.2f} | Q3: {q3:.2f} | IQR: {iqr:.2f}")
            logger.info(f"Outlier bounds: < {lower_bound:.2f} or > {upper_bound:.2f}")
            logger.info(f"Outliers found: {len(outliers)} ({stats_dict['outlier_pct']:.2f}%)")

        return stats_dict, outliers

    def generate_summary_tables(self) -> Dict[str, Any]:
        """
        Generates comparative summary statistics, crosstabs, and hypothesis test tables.
        Exports to CSV and JSON formats.
        """
        results = {}
        available_cols = self.bsa_columns

        for col in available_cols:
            summary_rows = []
            for class_name in self.CLASS_ORDER:
                sub_df = self.class_dfs[class_name]
                stat_res, _ = self.analyze_column_and_outliers(sub_df, col, verbose=False)
                stat_res['Oligomer_Class'] = class_name
                stat_res['Parameter'] = col
                summary_rows.append(stat_res)

            summary_df = pd.DataFrame(summary_rows)
            ordered_cols = ['Parameter', 'Oligomer_Class', 'total_count', 'mean', 'std', 'cv_pct',
                            'median', 'min', 'max', 'q1', 'q3', 'iqr', 'lower_bound', 'upper_bound',
                            'outlier_count', 'outlier_pct']
            summary_df = summary_df[ordered_cols]
            results[f"summary_{col}"] = summary_df

            csv_path = self.output_dir / f"summary_stats_{col}.csv"
            json_path = self.output_dir / f"summary_stats_{col}.json"
            summary_df.to_csv(csv_path, index=False)
            summary_df.to_json(json_path, orient="records", indent=2)
            logger.info(f"Saved summary table for {col} -> {csv_path}")

        # Flag outliers across the combined dataframe per Oligomer_Class
        for col in available_cols:
            def flag_outliers(s):
                clean_s = s.dropna()
                if len(clean_s) == 0:
                    return pd.Series(False, index=s.index)
                q1 = clean_s.quantile(0.25)
                q3 = clean_s.quantile(0.75)
                iqr = q3 - q1
                return (s < (q1 - 1.5 * iqr)) | (s > (q3 + 1.5 * iqr))

            self.df[f'{col}_Is_Outlier'] = self.df.groupby('Oligomer_Class')[col].transform(flag_outliers)

            # Crosstab counts
            ct_counts = pd.crosstab(
                self.df['Oligomer_Class'],
                self.df[f'{col}_Is_Outlier'].rename({False: 'Normal', True: 'Outlier'}),
                margins=True, margins_name="Total"
            )

            # Crosstab percentages
            ct_pct = pd.crosstab(
                self.df['Oligomer_Class'],
                self.df[f'{col}_Is_Outlier'].rename({False: 'Normal (%)', True: 'Outlier (%)'}),
                normalize='index'
            ) * 100

            results[f"crosstab_counts_{col}"] = ct_counts
            results[f"crosstab_pct_{col}"] = ct_pct

            ct_counts.to_csv(self.output_dir / f"crosstab_counts_{col}.csv")
            ct_pct.round(2).to_csv(self.output_dir / f"crosstab_pct_{col}.csv")

        # Statistical hypothesis testing (Kruskal-Wallis & ANOVA)
        test_rows = []
        for col in available_cols:
            groups = [sub_df[col].dropna().values for sub_df in self.class_dfs.values() if len(sub_df[col].dropna()) > 0]
            if len(groups) > 1:
                kw_stat, kw_p = stats.kruskal(*groups)
                anova_stat, anova_p = stats.f_oneway(*groups)
                test_rows.append({
                    'Parameter': col,
                    'Kruskal_Wallis_H': float(kw_stat),
                    'Kruskal_Wallis_p': float(kw_p),
                    'ANOVA_F': float(anova_stat),
                    'ANOVA_p': float(anova_p),
                    'Significant_Difference': bool(kw_p < 0.05)
                })

        test_df = pd.DataFrame(test_rows)
        results['statistical_tests'] = test_df
        test_csv = self.output_dir / "statistical_tests_overall.csv"
        test_json = self.output_dir / "statistical_tests_overall.json"
        test_df.to_csv(test_csv, index=False)
        test_df.to_json(test_json, orient="records", indent=2)
        logger.info(f"Saved statistical tests -> {test_csv} and {test_json}")

        return results

    @staticmethod
    def format_col_label(col_name: str) -> str:
        """Formats column name for clean display on plot axes (e.g. BSA_Complex -> BSA Complex)."""
        label = col_name.replace('_', ' ')
        label = label.replace('bsa', 'BSA').replace('Bsa', 'BSA')
        label = label.replace('dna', 'DNA').replace('Dna', 'DNA')
        label = label.replace('rna', 'RNA').replace('Rna', 'RNA')
        label = label.replace('protein', 'Protein')
        label = label.replace('Protein 1', 'Protein 1').replace('protein 2', 'Protein 2')
        return label

    def plot_bsa_distributions(self):
        """
        Generates publication-quality focused BSA Complex violin plot at 600 DPI.
        """
        available_cols = self.bsa_columns
        if not available_cols:
            logger.warning("No BSA columns available for plotting.")
            return

        active_classes = [c for c in self.CLASS_ORDER if c in self.df['Oligomer_Class'].values and len(self.df[self.df['Oligomer_Class'] == c]) > 0]
        palette = [self.CLASS_PALETTE[c] for c in active_classes]

        sns.set_theme(style="ticks", font_scale=1.15)
        matplotlib.rcParams['font.sans-serif'] = "DejaVu Sans"
        matplotlib.rcParams['font.family'] = "sans-serif"

        # Focused plot for primary BSA_Complex / Total BSA
        primary_col = next((c for c in available_cols if 'complex' in c.lower() or 'total' in c.lower()), available_cols[0])
        clean_df = self.df[['Oligomer_Class', primary_col]].dropna()

        fig, ax = plt.subplots(figsize=(8.0, 6.5))
        sns.violinplot(
            data=clean_df,
            x='Oligomer_Class',
            y=primary_col,
            order=active_classes,
            palette=palette,
            inner=None,
            cut=0,
            alpha=0.50,
            linewidth=1.6,
            ax=ax
        )
        sns.boxplot(
            data=clean_df,
            x='Oligomer_Class',
            y=primary_col,
            order=active_classes,
            width=0.20,
            boxprops=dict(facecolor='white', edgecolor='black', linewidth=1.4, alpha=0.92),
            medianprops=dict(color='black', linewidth=2.2),
            whiskerprops=dict(color='black', linewidth=1.3),
            capprops=dict(color='black', linewidth=1.3),
            showfliers=False,
            ax=ax
        )
        sns.stripplot(
            data=clean_df,
            x='Oligomer_Class',
            y=primary_col,
            order=active_classes,
            color='#111111',
            size=5.2,
            jitter=0.18,
            alpha=0.80,
            linewidth=0.35,
            edgecolor='black',
            ax=ax
        )

        # Annotate Mean ± SD and sample size (N)
        y_max = clean_df[primary_col].max()
        y_min = clean_df[primary_col].min()
        y_range = y_max - y_min
        ax.set_ylim(bottom=max(0, y_min - 0.05 * y_range), top=y_max + 0.20 * y_range)

        for i, c in enumerate(active_classes):
            sub_vals = clean_df[clean_df['Oligomer_Class'] == c][primary_col].dropna()
            m_val = sub_vals.mean()
            sd_val = sub_vals.std()
            n_val = len(sub_vals)

            badge_text = f"Mean: {round(m_val):,} ± {round(sd_val):,} Å²\n(N = {n_val})"
            ax.text(
                i,
                y_max + 0.04 * y_range,
                badge_text,
                ha='center',
                va='bottom',
                fontsize=11.0,
                fontweight='bold',
                bbox=dict(
                    boxstyle='round,pad=0.45',
                    facecolor='white',
                    edgecolor=self.CLASS_PALETTE.get(c, '#333333'),
                    alpha=0.95,
                    linewidth=1.6
                )
            )

        clean_y_label = self.format_col_label(primary_col)
        ax.set_xlabel("Oligomeric Class", fontsize=13, fontweight='bold')
        ax.set_ylabel(f"{clean_y_label} (Å²)", fontsize=13, fontweight='bold')
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=11)
        ax.grid(True, linestyle='--', alpha=0.4)
        sns.despine(top=True, right=True)
        plt.tight_layout()

        # Save focused violin plots (both specific name and Figure 5 filename)
        out_png = self.figures_dir / f"bsa_complex_violin_{self.dataset_name}.png"
        out_pdf = self.figures_dir / f"bsa_complex_violin_{self.dataset_name}.pdf"
        plt.savefig(out_png, dpi=600, bbox_inches='tight')
        plt.savefig(out_pdf, dpi=600, bbox_inches='tight')

        fig5_single_png = self.figures_dir / f"Figure 5_{self.dataset_name}.png"
        fig5_single_pdf = self.figures_dir / f"Figure 5_{self.dataset_name}.pdf"
        plt.savefig(fig5_single_png, dpi=600, bbox_inches='tight')
        plt.savefig(fig5_single_pdf, dpi=600, bbox_inches='tight')
        plt.close()

        logger.info(f"Saved focused violin plot -> {out_png}")
        logger.info(f"Saved Figure 5 dataset plot -> {fig5_single_png}")


def plot_combined_figure_5(na_analyzer: BSAOligomerAnalyzer, pp_analyzer: BSAOligomerAnalyzer, output_dirs: List[Path]):
    """
    Generates a unified publication Figure 5 combining:
      - (A) Protein-Nucleic Acid Complexes BSA_Complex
      - (B) Protein-Protein Complexes BSA_Complex
    Side-by-side in a crisp 2-panel layout at 600 DPI.
    """
    sns.set_theme(style="ticks", font_scale=1.15)
    matplotlib.rcParams['font.sans-serif'] = "DejaVu Sans"
    matplotlib.rcParams['font.family'] = "sans-serif"

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.8))

    analyzers = [
        (na_analyzer, axes[0], "A", "Protein–Nucleic Acid Complexes"),
        (pp_analyzer, axes[1], "B", "Protein–Protein Complexes")
    ]

    for analyzer, ax, panel_tag, panel_title in analyzers:
        if analyzer.df is None or analyzer.df.empty:
            continue

        primary_col = next((c for c in analyzer.bsa_columns if 'complex' in c.lower() or 'total' in c.lower()), analyzer.bsa_columns[0])
        clean_df = analyzer.df[['Oligomer_Class', primary_col]].dropna()

        active_classes = [c for c in analyzer.CLASS_ORDER if c in clean_df['Oligomer_Class'].values and len(clean_df[clean_df['Oligomer_Class'] == c]) > 0]
        palette = [analyzer.CLASS_PALETTE[c] for c in active_classes]

        # 1. Violin Plot
        sns.violinplot(
            data=clean_df,
            x='Oligomer_Class',
            y=primary_col,
            order=active_classes,
            palette=palette,
            inner=None,
            cut=0,
            alpha=0.50,
            linewidth=1.6,
            ax=ax
        )

        # 2. Box Plot
        sns.boxplot(
            data=clean_df,
            x='Oligomer_Class',
            y=primary_col,
            order=active_classes,
            width=0.20,
            boxprops=dict(facecolor='white', edgecolor='black', linewidth=1.4, alpha=0.92),
            medianprops=dict(color='black', linewidth=2.2),
            whiskerprops=dict(color='black', linewidth=1.3),
            capprops=dict(color='black', linewidth=1.3),
            showfliers=False,
            ax=ax
        )

        # 3. Strip Plot
        sns.stripplot(
            data=clean_df,
            x='Oligomer_Class',
            y=primary_col,
            order=active_classes,
            color='#111111',
            size=5.2,
            jitter=0.18,
            alpha=0.80,
            linewidth=0.35,
            edgecolor='black',
            ax=ax
        )

        # Annotate Mean ± SD and sample size (N)
        y_max = clean_df[primary_col].max()
        y_min = clean_df[primary_col].min()
        y_range = y_max - y_min
        ax.set_ylim(bottom=max(0, y_min - 0.05 * y_range), top=y_max + 0.22 * y_range)

        for i, c in enumerate(active_classes):
            sub_vals = clean_df[clean_df['Oligomer_Class'] == c][primary_col].dropna()
            m_val = sub_vals.mean()
            sd_val = sub_vals.std()
            n_val = len(sub_vals)

            badge_text = f"Mean: {round(m_val):,} ± {round(sd_val):,} Å²\n(N = {n_val})"
            ax.text(
                i,
                y_max + 0.04 * y_range,
                badge_text,
                ha='center',
                va='bottom',
                fontsize=11.0,
                fontweight='bold',
                bbox=dict(
                    boxstyle='round,pad=0.45',
                    facecolor='white',
                    edgecolor=analyzer.CLASS_PALETTE.get(c, '#333333'),
                    alpha=0.95,
                    linewidth=1.6
                )
            )

        ax.set_title(f"({panel_tag}) {panel_title}", fontsize=14, fontweight='bold', pad=12)
        ax.set_xlabel("Oligomeric Class", fontsize=13, fontweight='bold')
        ax.set_ylabel("Total Complex BSA (Å²)", fontsize=13, fontweight='bold')
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=11)
        ax.grid(True, linestyle='--', alpha=0.4)
        sns.despine(top=True, right=True, ax=ax)

    plt.tight_layout()

    for out_d in output_dirs:
        out_d.mkdir(parents=True, exist_ok=True)
        fig5_png = out_d / "Figure_5.png"
        fig5_pdf = out_d / "Figure_5.pdf"
        plt.savefig(fig5_png, dpi=600, bbox_inches='tight')
        plt.savefig(fig5_pdf, dpi=600, bbox_inches='tight')
        logger.info(f"Saved Unified 2-Panel Figure 5 -> {fig5_png}")

    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description="BSA Distribution & Outlier Analysis across Oligomeric States (Figure 5 Generator).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=['protein_nucleic_acid', 'protein_protein', 'both'],
        default='both',
        help="Dataset mode: 'protein_nucleic_acid', 'protein_protein', or 'both'."
    )
    parser.add_argument(
        "--input_na",
        type=Path,
        default=DEFAULT_INPUT_NA,
        help="Path to Protein-Nucleic Acid Excel dataset."
    )
    parser.add_argument(
        "--input_pp",
        type=Path,
        default=DEFAULT_INPUT_PP,
        help="Path to Protein-Protein Excel dataset."
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Base directory where output folders will be saved."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    analyzers_dict = {}

    if args.mode in ['protein_nucleic_acid', 'both']:
        na_analyzer = BSAOligomerAnalyzer(excel_path=args.input_na, dataset_name='protein_nucleic_acid', base_output_dir=args.output_dir)
        na_analyzer.load_and_preprocess()
        na_analyzer.generate_summary_tables()
        na_analyzer.plot_bsa_distributions()
        analyzers_dict['protein_nucleic_acid'] = na_analyzer

    if args.mode in ['protein_protein', 'both']:
        pp_analyzer = BSAOligomerAnalyzer(excel_path=args.input_pp, dataset_name='protein_protein', base_output_dir=args.output_dir)
        pp_analyzer.load_and_preprocess()
        pp_analyzer.generate_summary_tables()
        pp_analyzer.plot_bsa_distributions()
        analyzers_dict['protein_protein'] = pp_analyzer

    # If both are run, generate the unified 2-panel Figure 5
    if 'protein_nucleic_acid' in analyzers_dict and 'protein_protein' in analyzers_dict:
        out_dirs = [
            Path(args.output_dir) / "Figures",
            PROJECT_ROOT / "results" / "Figures"
        ]
        plot_combined_figure_5(analyzers_dict['protein_nucleic_acid'], analyzers_dict['protein_protein'], out_dirs)

    logger.info("All Figure 5 visualizations and statistical analyses completed successfully!")


if __name__ == '__main__':
    main()
