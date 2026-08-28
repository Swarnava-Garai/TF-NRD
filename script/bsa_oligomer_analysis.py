#!/usr/bin/env python3
"""
bsa_oligomer_analysis.py
========================
Comprehensive statistical analysis and visualization suite for Buried Surface Area (BSA)
parameters categorized across 3 oligomeric classes:
  1. Monomer
  2. Homodimer
  3. Heterodimer

Note: Higher order oligomers are excluded from analysis because BSA for higher order
oligomers was calculated only for selected chains which are Transcription Factors.

Supports both:
  - Protein-Nucleic Acid complexes ('TF_nucleic_acid_whole_with_oligostate.xlsx')
  - Protein-Protein complexes ('TF_protein_protein_with_oligostate.xlsx')

Results and figures are saved into dedicated separate subdirectories under results/:
  - results/protein_nucleic_acid/ (and Figures/)
  - results/protein_protein/ (and Figures/)

Author: Antigravity AI / PhD Projects Team
Location: TF-NRD/script/bsa_oligomer_analysis.py
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)


def resolve_path(path_str: str) -> str:
    """
    Intelligently resolves relative/absolute file paths whether invoked from repo root, 
    script dir, or workspace root.
    """
    if not path_str:
        return path_str

    if os.path.isabs(path_str) and os.path.exists(path_str):
        return path_str

    clean_path = path_str.strip()

    candidates = [
        os.path.abspath(clean_path),
        os.path.join(REPO_ROOT, clean_path),
        os.path.join(REPO_ROOT, clean_path.replace("TF-NRD/", "", 1)) if clean_path.startswith("TF-NRD/") else clean_path,
        os.path.join(os.getcwd(), clean_path),
        os.path.join(os.getcwd(), clean_path.replace("TF-NRD/", "", 1)) if clean_path.startswith("TF-NRD/") else clean_path,
    ]

    for cand in candidates:
        if os.path.exists(cand):
            return os.path.abspath(cand)

    return os.path.abspath(clean_path)


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

    def __init__(self, excel_path: str, dataset_name: str = "protein_nucleic_acid", base_output_dir: str = "results/Interface"):
        self.excel_path = resolve_path(excel_path)
        self.dataset_name = dataset_name
        
        base_out = resolve_path(base_output_dir)
        if not os.path.exists(os.path.dirname(base_out)):
            base_out = os.path.abspath(os.path.join(REPO_ROOT, base_output_dir.replace("TF-NRD/", "", 1) if base_output_dir.startswith("TF-NRD/") else base_output_dir))
            
        self.output_dir = os.path.join(base_out, dataset_name)
        self.figures_dir = os.path.join(self.output_dir, "Figures")
        
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.figures_dir, exist_ok=True)
        
        self.df = None
        self.monomer_df = None
        self.homodimer_df = None
        self.heterodimer_df = None
        self.class_dfs = {}
        self.bsa_columns = []

    def load_and_preprocess(self) -> pd.DataFrame:
        """
        Loads dataset, detects BSA columns, and categorizes records into 3 distinct oligomeric classes.
        Higher order oligomers are excluded.
        """
        if not os.path.exists(self.excel_path):
            raise FileNotFoundError(f"Input file not found at: {self.excel_path}")

        print(f"\n=======================================================")
        print(f" Processing Dataset: {self.dataset_name.upper()} ")
        print(f" File: {self.excel_path}")
        print(f" Output directory: {self.output_dir}")
        print(f"=======================================================")

        self.df = pd.read_excel(self.excel_path)
        
        # Identify state column ('States' or 'Oligomeric_State')
        state_col = None
        for col_name in ['States', 'Oligomeric_State', 'Oligostate', 'State']:
            if col_name in self.df.columns:
                state_col = col_name
                break

        if state_col is None:
            raise KeyError(f"Could not find an oligomeric state column in: {self.df.columns.tolist()}")

        print(f"[INFO] Using state column '{state_col}'.")
        print(f"[INFO] Raw unique states: {self.df[state_col].dropna().unique().tolist()}")

        # Dynamically detect all BSA columns (e.g. BSA_Complex, BSA_Protein, BSA_Protein_1, BSA_protein_2, BSA_DNA, BSA_RNA)
        self.bsa_columns = [c for c in self.df.columns if c.lower().startswith('bsa_')]
        print(f"[INFO] Detected BSA parameters: {self.bsa_columns}")

        # Coerce BSA columns to numeric (handling any strings like 'None' or '-')
        for bcol in self.bsa_columns:
            self.df[bcol] = pd.to_numeric(self.df[bcol], errors='coerce')

        # Assign 3 target classes (Higher order oligomers are excluded)
        def classify_state(state):
            if pd.isna(state):
                return np.nan
            state_str = str(state).strip()
            if state_str == 'Monomer':
                return 'Monomer'
            elif state_str == 'Homodimer':
                return 'Homodimer'
            elif state_str == 'Heterodimer':
                return 'Heterodimer'
            else:
                return np.nan

        self.df['Oligomer_Class'] = self.df[state_col].apply(classify_state)

        # Exclude higher order oligomers / unclassified states from analysis
        before_count = len(self.df)
        self.df = self.df[self.df['Oligomer_Class'].isin(self.CLASS_ORDER)].copy()
        print(f"[INFO] Filtered out {before_count - len(self.df)} higher order oligomer / unclassified records. Remaining: {len(self.df)}")

        # Slice dataframes for the 3 target classes
        self.monomer_df = self.df[self.df['Oligomer_Class'] == 'Monomer'].copy()
        self.homodimer_df = self.df[self.df['Oligomer_Class'] == 'Homodimer'].copy()
        self.heterodimer_df = self.df[self.df['Oligomer_Class'] == 'Heterodimer'].copy()

        self.class_dfs = {
            'Monomer': self.monomer_df,
            'Homodimer': self.homodimer_df,
            'Heterodimer': self.heterodimer_df
        }

        print("\n[INFO] Class record counts:")
        for name in self.CLASS_ORDER:
            sub_df = self.class_dfs[name]
            print(f"  - {name:22s}: {len(sub_df)} records")

        return self.df

    @staticmethod
    def analyze_column_and_outliers(df: pd.DataFrame, column: str, iqr_factor: float = 1.5, verbose: bool = False):
        """
        Calculates summary statistics and IQR-based outliers for a specific DataFrame column.

        Parameters:
        -----------
        df : pd.DataFrame
            Input DataFrame.
        column : str
            Column name to analyze.
        iqr_factor : float
            IQR multiplier (default 1.5).
        verbose : bool
            Whether to print human-readable output.

        Returns:
        --------
        stats : dict
            Summary metrics dictionary.
        outliers : pd.DataFrame
            Outlier rows subset.
        """
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
            print(f"--- Analysis for '{column}' ---")
            print(f"Count: {n} | Mean: {mean_val:.2f} | Std: {std_val:.2f} | CV (%): {cv_val:.2f}%")
            print(f"Median: {median_val:.2f} | Q1: {q1:.2f} | Q3: {q3:.2f} | IQR: {iqr:.2f}")
            print(f"Outlier bounds: < {lower_bound:.2f} or > {upper_bound:.2f}")
            print(f"Outliers found: {len(outliers)} ({stats_dict['outlier_pct']:.2f}%)\n")

        return stats_dict, outliers

    def generate_summary_tables(self) -> dict:
        """
        Generates comparative summary statistics, crosstabs, and hypothesis test tables.
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

            csv_path = os.path.join(self.output_dir, f"summary_stats_{col}.csv")
            summary_df.to_csv(csv_path, index=False)
            print(f"[SAVED] Summary table for {col} -> {csv_path}")

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

            ct_counts.to_csv(os.path.join(self.output_dir, f"crosstab_counts_{col}.csv"))
            ct_pct.round(2).to_csv(os.path.join(self.output_dir, f"crosstab_pct_{col}.csv"))

        # Statistical hypothesis testing (Kruskal-Wallis & ANOVA)
        test_rows = []
        for col in available_cols:
            groups = [sub_df[col].dropna().values for sub_df in self.class_dfs.values() if len(sub_df[col].dropna()) > 0]
            if len(groups) > 1:
                kw_stat, kw_p = stats.kruskal(*groups)
                anova_stat, anova_p = stats.f_oneway(*groups)
                test_rows.append({
                    'Parameter': col,
                    'Kruskal_Wallis_H': kw_stat,
                    'Kruskal_Wallis_p': kw_p,
                    'ANOVA_F': anova_stat,
                    'ANOVA_p': anova_p,
                    'Significant_Difference': (kw_p < 0.05)
                })

        test_df = pd.DataFrame(test_rows)
        results['statistical_tests'] = test_df
        test_df.to_csv(os.path.join(self.output_dir, "statistical_tests_overall.csv"), index=False)
        print(f"[SAVED] Overall statistical tests -> {os.path.join(self.output_dir, 'statistical_tests_overall.csv')}")

        return results

    @staticmethod
    def format_col_label(col_name: str) -> str:
        """Formats column name for clean display on plot axes (e.g. BSA_Complex -> BSA Complex)."""
        label = col_name.replace('_', ' ')
        label = label.replace('bsa', 'BSA').replace('Bsa', 'BSA')
        label = label.replace('dna', 'DNA').replace('Dna', 'DNA')
        label = label.replace('rna', 'RNA').replace('Rna', 'RNA')
        label = label.replace('protein', 'Protein')
        label = label.replace('Protein 1', 'Protein 1').replace('protein 2', 'Protein 2').replace('Protein 2', 'Protein 2')
        return label

    def plot_bsa_distributions(self):
        """
        Generates publication-quality Violin plots overlaid with Boxplots and jittered points.
        """
        available_cols = self.bsa_columns
        num_cols = len(available_cols)

        if num_cols == 0:
            print("[WARN] No BSA columns found to plot.")
            return

        sns.set_theme(style="whitegrid", font_scale=1.1)

        # Determine subplot grid layout
        if num_cols == 1:
            nrows, ncols = 1, 1
            fig_size = (8, 6)
        elif num_cols == 2:
            nrows, ncols = 1, 2
            fig_size = (14, 6)
        elif num_cols == 3:
            nrows, ncols = 1, 3
            fig_size = (18, 5.5)
        elif num_cols == 4:
            nrows, ncols = 2, 2
            fig_size = (16, 12)
        else:
            ncols = 3
            nrows = (num_cols + 2) // 3
            fig_size = (18, 5.5 * nrows)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=fig_size, squeeze=False)
        axes = axes.flatten()

        for idx, col in enumerate(available_cols):
            ax = axes[idx]
            clean_df = self.df.dropna(subset=[col, 'Oligomer_Class'])

            # Dynamically determine active classes with data (N >= 2)
            active_classes = [c for c in self.CLASS_ORDER if c in clean_df['Oligomer_Class'].values and len(clean_df[clean_df['Oligomer_Class'] == c]) >= 2]

            # 1. Violin Plot
            sns.violinplot(
                data=clean_df,
                x='Oligomer_Class',
                y=col,
                order=active_classes,
                hue='Oligomer_Class',
                legend=False,
                palette=self.CLASS_PALETTE,
                inner=None,
                cut=0,
                density_norm='width',
                alpha=0.55,
                ax=ax
            )

            # 2. Overlay Boxplot inside Violin
            sns.boxplot(
                data=clean_df,
                x='Oligomer_Class',
                y=col,
                order=active_classes,
                width=0.18,
                showcaps=True,
                boxprops=dict(facecolor='white', edgecolor='black', alpha=0.9, linewidth=1.6),
                whiskerprops=dict(color='black', linewidth=1.6),
                capprops=dict(color='black', linewidth=1.6),
                medianprops=dict(color='#d90429', linewidth=2.8),
                showfliers=False,
                ax=ax
            )

            # 3. Jittered Strip plot for individual points (Solid, crisp publication-ready dots)
            sns.stripplot(
                data=clean_df,
                x='Oligomer_Class',
                y=col,
                order=active_classes,
                color='#111111',
                size=4.0,
                jitter=0.18,
                alpha=0.75,
                linewidth=0.3,
                edgecolor='black',
                ax=ax
            )

            clean_y_label = self.format_col_label(col)
            ax.set_xlabel("Oligomeric Class", fontsize=12, fontweight='bold')
            ax.set_ylabel(f"{clean_y_label} (Å²)", fontsize=12, fontweight='bold')
            ax.tick_params(axis='x', rotation=15 if len(active_classes) > 2 else 0)

        # Hide extra subplots if grid has empty cells
        for idx in range(num_cols, len(axes)):
            fig.delaxes(axes[idx])

        plt.tight_layout()

        multi_png = os.path.join(self.figures_dir, f"bsa_distributions_{self.dataset_name}.png")
        multi_pdf = os.path.join(self.figures_dir, f"bsa_distributions_{self.dataset_name}.pdf")
        plt.savefig(multi_png, dpi=300, bbox_inches='tight')
        plt.savefig(multi_pdf, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[SAVED] Multi-panel Violin Plot -> {multi_png}")

        # Dedicated Focused Plot for BSA_Complex
        if 'BSA_Complex' in self.df.columns:
            plt.figure(figsize=(9, 6.5))
            col = 'BSA_Complex'
            clean_df = self.df.dropna(subset=[col, 'Oligomer_Class'])

            active_classes = [c for c in self.CLASS_ORDER if c in clean_df['Oligomer_Class'].values and len(clean_df[clean_df['Oligomer_Class'] == c]) >= 2]

            ax = sns.violinplot(
                data=clean_df,
                x='Oligomer_Class',
                y=col,
                order=active_classes,
                hue='Oligomer_Class',
                legend=False,
                palette=self.CLASS_PALETTE,
                inner=None,
                cut=0,
                density_norm='width',
                alpha=0.55
            )

            sns.boxplot(
                data=clean_df,
                x='Oligomer_Class',
                y=col,
                order=active_classes,
                width=0.15,
                showcaps=True,
                boxprops=dict(facecolor='white', edgecolor='black', alpha=0.9, linewidth=1.8),
                whiskerprops=dict(color='black', linewidth=1.8),
                capprops=dict(color='black', linewidth=1.8),
                medianprops=dict(color='#d90429', linewidth=3.0),
                showfliers=False,
                ax=ax
            )

            sns.stripplot(
                data=clean_df,
                x='Oligomer_Class',
                y=col,
                order=active_classes,
                color='#111111',
                size=4.8,
                jitter=0.18,
                alpha=0.8,
                linewidth=0.35,
                edgecolor='black',
                ax=ax
            )

            clean_y_label = self.format_col_label(col)
            plt.xlabel("Oligomeric Class", fontsize=13, fontweight='bold')
            plt.ylabel(f"{clean_y_label} (Å²)", fontsize=13, fontweight='bold')
            plt.xticks(fontsize=11, fontweight='bold')
            plt.yticks(fontsize=11)
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()

            single_png = os.path.join(self.figures_dir, f"bsa_complex_violin_{self.dataset_name}.png")
            single_pdf = os.path.join(self.figures_dir, f"bsa_complex_violin_{self.dataset_name}.pdf")
            plt.savefig(single_png, dpi=300, bbox_inches='tight')
            plt.savefig(single_pdf, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[SAVED] Focused BSA_Complex Violin Plot -> {single_png}")


def main():
    parser = argparse.ArgumentParser(description="BSA Distribution & Outlier Analysis across Oligomeric States.")
    parser.add_argument(
        "--mode", 
        type=str,
        choices=['protein_nucleic_acid', 'protein_protein', 'both'],
        default='both',
        help="Dataset mode: 'protein_nucleic_acid', 'protein_protein', or 'both'."
    )
    parser.add_argument(
        "--input_na", 
        type=str, 
        default="input_data/bsa/TF_nucleic_acid_whole_with_oligostate.xlsx",
        help="Path to Protein-Nucleic Acid Excel dataset."
    )
    parser.add_argument(
        "--input_pp", 
        type=str, 
        default="input_data/bsa/TF_protein_protein_with_oligostate.xlsx",
        help="Path to Protein-Protein Excel dataset."
    )
    parser.add_argument(
        "--output_dir", 
        type=str, 
        default="results/Interface",
        help="Base directory where output folders (protein_nucleic_acid, protein_protein) will be saved."
    )

    args = parser.parse_args()

    tasks = []
    if args.mode in ['protein_nucleic_acid', 'both']:
        tasks.append(('protein_nucleic_acid', args.input_na))
    if args.mode in ['protein_protein', 'both']:
        tasks.append(('protein_protein', args.input_pp))

    for dataset_name, input_file in tasks:
        analyzer = BSAOligomerAnalyzer(excel_path=input_file, dataset_name=dataset_name, base_output_dir=args.output_dir)
        analyzer.load_and_preprocess()
        analyzer.generate_summary_tables()
        analyzer.plot_bsa_distributions()

    print("\n[SUCCESS] All dataset analyses and visualizations completed successfully!")


if __name__ == '__main__':
    main()
