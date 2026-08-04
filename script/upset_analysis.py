#!/usr/bin/env python3
"""
TF-NRD UpSet Plot Analysis Suite
=================================
Processes domain and motif distributions across subcellular localization categories
('ON', 'OC', 'NC', 'NO', 'CO', 'OO') and generates publication-quality UpSet plots.

Outputs high-resolution 600 DPI PNG and vector SVG plots into:
  - results/Figures/TFNRD_Domain_UpSet.png (.svg)
  - results/Figures/TFNRD_Motif_distribution.png (.svg)

Author: Antigravity Team / PhD Project Suite
"""

import argparse
import logging
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import upsetplot.plotting as upp
from upsetplot import UpSet, from_contents

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("UpSetAnalyzer")

# --- Patch upsetplot for Pandas 3.0 / Matplotlib 3.11 compatibility ---
_orig_plot_matrix = upp.UpSet.plot_matrix

def _safe_plot_matrix(self, ax):
    """
    Custom plot_matrix wrapper to fix Pandas 3.0 Copy-on-Write inplace fillna bug
    and Matplotlib 3.11 text rendering array scalar conversion error.
    """
    ax = self._reorient(ax)
    data = self.intersections
    n_cats = data.index.nlevels
    inclusion = data.index.to_frame().values

    other_dots = str(self._other_dots_color) if isinstance(self._other_dots_color, (int, float)) else self._other_dots_color

    styles = [
        [
            self.subset_styles[i]
            if inclusion[i, j]
            else {'facecolor': other_dots, 'linewidth': 0}
            for j in range(n_cats)
        ]
        for i in range(len(data))
    ]
    styles = sum(styles, [])
    style_columns = {
        'facecolor': 'facecolors',
        'edgecolor': 'edgecolors',
        'linewidth': 'linewidths',
        'linestyle': 'linestyles',
    }
    df_styles = pd.DataFrame(styles)
    for col in style_columns.keys():
        if col not in df_styles.columns:
            df_styles[col] = np.nan

    df_styles['linewidth'] = df_styles['linewidth'].fillna(1)
    df_styles['facecolor'] = df_styles['facecolor'].fillna(self._facecolor)
    df_styles['edgecolor'] = df_styles['edgecolor'].fillna(df_styles['facecolor'])
    df_styles['linestyle'] = df_styles['linestyle'].fillna('solid')

    x = np.repeat(np.arange(len(data)), n_cats)
    y = np.tile(np.arange(n_cats), len(data))

    if self._element_size is not None:
        s = (self._element_size * 0.35) ** 2
    else:
        s = 200

    kw = {style_columns[col]: df_styles[col].values for col in style_columns}
    ax.scatter(
        *self._swapaxes(x, y),
        s=s,
        zorder=10,
        **kw,
    )

    for i in range(len(data)):
        active = np.where(inclusion[i])[0]
        if len(active) > 1:
            ax.plot(
                *self._swapaxes([i, i], [active[0], active[-1]]),
                color=self._facecolor,
                linewidth=2,
                zorder=9,
            )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

upp.UpSet.plot_matrix = _safe_plot_matrix


def _fix_text_positions(plot_res):
    """Fix array text positions for Matplotlib 3.11 compatibility."""
    for ax_name, ax in plot_res.items():
        for t in ax.texts:
            x, y = t.get_position()
            t.set_position((float(np.ravel(x)[0]), float(np.ravel(y)[0])))


class UpSetAnalyzer:
    """
    Feature analysis suite for domain and motif UpSet plot generation across
    subcellular localization categories.
    """

    def __init__(self, input_dir: str = "input_data", output_dir: str = "results/Figures"):
        self.base_dir = Path(__file__).resolve().parent.parent
        self.input_dir = self.base_dir / input_dir
        self.output_dir = self.base_dir / output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot_domain_upset(self) -> Path:
        """
        Maps structure subcellular localization categories ('ON', 'OC', 'NC', 'NO', 'CO', 'OO')
        with PFAM domain accessions from standard_start_end_domain.xlsx and plots UpSet plot.
        """
        str_subcell_file = self.input_dir / 'subcellular_location' / 'Structure_subcellular_location_detailed.xlsx'
        dom_file = self.input_dir / 'domain' / 'standard_start_end_domain.xlsx'

        if not str_subcell_file.exists() or not dom_file.exists():
            logger.warning(f"Required files for Domain UpSet plot not found: {str_subcell_file} or {dom_file}")
            return None

        str_subcell = pd.read_excel(str_subcell_file)
        dom_df = pd.read_excel(dom_file)

        on_pdbs = set(str_subcell['Only_Nucleus'].dropna())
        oc_pdbs = set(str_subcell['Only_Cytoplasm'].dropna())
        nc_pdbs = set(str_subcell['Nucleus_Cytoplasm_Both'].dropna())
        no_pdbs = set(str_subcell['Nucleus_with_Other_not_Cytoplasm'].dropna())
        co_pdbs = set(str_subcell['Cytoplasm_with_Other_not_Nucleus'].dropna())

        def get_str_loc(pdb):
            if pdb in on_pdbs: return 'ON'
            if pdb in oc_pdbs: return 'OC'
            if pdb in nc_pdbs: return 'NC'
            if pdb in no_pdbs: return 'NO'
            if pdb in co_pdbs: return 'CO'
            return 'OO'

        dom_df['Subcellular_Location'] = dom_df['PDB'].apply(get_str_loc)

        domain_sets = {
            'ON': set(dom_df.loc[dom_df['Subcellular_Location'] == 'ON', 'PFAM_ACCESSION'].dropna()),
            'OC': set(dom_df.loc[dom_df['Subcellular_Location'] == 'OC', 'PFAM_ACCESSION'].dropna()),
            'NC': set(dom_df.loc[dom_df['Subcellular_Location'] == 'NC', 'PFAM_ACCESSION'].dropna()),
            'NO': set(dom_df.loc[dom_df['Subcellular_Location'] == 'NO', 'PFAM_ACCESSION'].dropna()),
            'CO': set(dom_df.loc[dom_df['Subcellular_Location'] == 'CO', 'PFAM_ACCESSION'].dropna()),
            'OO': set(dom_df.loc[dom_df['Subcellular_Location'] == 'OO', 'PFAM_ACCESSION'].dropna())
        }

        domain_sets = {k: v for k, v in domain_sets.items() if len(v) > 0}
        if not domain_sets:
            logger.warning("No domain data available for UpSet plot")
            return None

        upset_data = from_contents(domain_sets)

        plt.rcParams.update({
            "font.size": 14,
            "font.weight": "regular",
            "axes.labelsize": 16,
            "axes.titlesize": 18,
            "figure.dpi": 600
        })

        fig = plt.figure(figsize=(12, 6))
        upset = UpSet(
            upset_data,
            orientation='horizontal',
            show_counts='%d',
            sort_by='cardinality',
            min_subset_size=1,
            facecolor='darkblue',
            shading_color='lightgray',
            element_size=None,
            other_dots_color='0.4',
            intersection_plot_elements=10,
            totals_plot_elements=2
        )

        upset.style_subsets(min_degree=1, facecolor="firebrick", edgecolor="black", linewidth=0.5)
        upset.style_subsets(min_degree=2, facecolor="darkgreen", edgecolor="black", linewidth=0.05)
        upset.style_subsets(min_degree=3, facecolor="gold", edgecolor="black", linewidth=0.5)

        plot_res = upset.plot(fig)
        _fix_text_positions(plot_res)

        plot_res['intersections'].set_ylabel('Number of PFAM Domains', fontweight='bold', fontsize=14)
        plot_res['totals'].set_xlabel('Total Domains \nper Category', fontweight='bold', fontsize=12)

        p_png = self.output_dir / 'TFNRD_Domain_UpSet.png'
        p_svg = self.output_dir / 'TFNRD_Domain_UpSet.svg'

        plt.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.15)
        plt.savefig(p_png, format='png', dpi=600, bbox_inches='tight', facecolor='white')
        plt.savefig(p_svg, format='svg', dpi=600, bbox_inches='tight', facecolor='white')
        plt.close()

        logger.info(f"Successfully generated Domain UpSet Plot: {p_png}")
        return p_png

    def plot_motif_upset(self) -> Path:
        """
        Maps sequence subcellular location categories ('ON', 'OC', 'NC', 'NO', 'CO', 'OO')
        with InterPro motif accessions from nr_sequence_dataset_motif_details.xlsx and plots UpSet plot.
        """
        motif_file = self.input_dir / 'motif' / 'nr_sequence_dataset_motif_details.xlsx'
        seq_subcell_file = self.input_dir / 'subcellular_location' / 'Sequence_subcellular_location_detailed.xlsx'

        if not motif_file.exists() or not seq_subcell_file.exists():
            logger.warning(f"Required files for Motif UpSet plot not found: {motif_file} or {seq_subcell_file}")
            return None

        seq_subcell = pd.read_excel(seq_subcell_file)
        df_m = pd.read_excel(motif_file)
        df_m['Entry_clean'] = df_m['Entry'].astype(str).str.strip()

        seq_on = set(seq_subcell['Only_Nucleus'].dropna().astype(str).str.strip())
        seq_oc = set(seq_subcell['Only_Cytoplasm'].dropna().astype(str).str.strip())
        seq_nc = set(seq_subcell['Nucleus_Cytoplasm_Both'].dropna().astype(str).str.strip())
        seq_no = set(seq_subcell['Nucleus_with_Other_not_Cytoplasm'].dropna().astype(str).str.strip())
        seq_co = set(seq_subcell['Cytoplasm_with_Other_not_Nucleus'].dropna().astype(str).str.strip())

        def get_seq_loc(entry):
            if entry in seq_on: return 'ON'
            if entry in seq_oc: return 'OC'
            if entry in seq_nc: return 'NC'
            if entry in seq_no: return 'NO'
            if entry in seq_co: return 'CO'
            return 'OO'

        df_m['Subcellular_Location'] = df_m['Entry_clean'].apply(get_seq_loc)
        df_m['Inter_Pro'] = df_m['InterPro'].astype(str).str.split(';').str[0].str.strip()
        df_m = df_m[df_m['Inter_Pro'].notna() & (df_m['Inter_Pro'] != 'nan')]

        motif_sets = {
            'ON': set(df_m.loc[df_m['Subcellular_Location'] == 'ON', 'Inter_Pro']),
            'OC': set(df_m.loc[df_m['Subcellular_Location'] == 'OC', 'Inter_Pro']),
            'NC': set(df_m.loc[df_m['Subcellular_Location'] == 'NC', 'Inter_Pro']),
            'NO': set(df_m.loc[df_m['Subcellular_Location'] == 'NO', 'Inter_Pro']),
            'CO': set(df_m.loc[df_m['Subcellular_Location'] == 'CO', 'Inter_Pro']),
            'OO': set(df_m.loc[df_m['Subcellular_Location'] == 'OO', 'Inter_Pro'])
        }

        motif_sets = {k: v for k, v in motif_sets.items() if len(v) > 0}
        if not motif_sets:
            logger.warning("No motif data available for UpSet plot")
            return None

        upset_data_m = from_contents(motif_sets)

        plt.rcParams.update({
            "font.size": 14,
            "font.weight": "regular",
            "axes.labelsize": 16,
            "axes.titlesize": 18,
            "figure.dpi": 600
        })

        fig = plt.figure(figsize=(12, 6))
        upset_m = UpSet(
            upset_data_m,
            orientation='horizontal',
            show_counts='%d',
            sort_by='cardinality',
            min_subset_size=1,
            facecolor='darkblue',
            shading_color='lightgray',
            element_size=None,
            other_dots_color='0.4',
            intersection_plot_elements=10,
            totals_plot_elements=2
        )

        upset_m.style_subsets(min_degree=1, facecolor="firebrick", edgecolor="black", linewidth=0.5)
        upset_m.style_subsets(min_degree=2, facecolor="darkgreen", edgecolor="black", linewidth=0.05)
        upset_m.style_subsets(min_degree=3, facecolor="gold", edgecolor="black", linewidth=0.5)

        plot_res_m = upset_m.plot(fig)
        _fix_text_positions(plot_res_m)

        plot_res_m['intersections'].set_ylabel('Number of unique Motifs', fontweight='bold', fontsize=14)
        plot_res_m['totals'].set_xlabel('Total Motifs \nper Category', fontweight='bold', fontsize=12)

        p_png = self.output_dir / 'TFNRD_Motif_distribution.png'
        p_svg = self.output_dir / 'TFNRD_Motif_distribution.svg'

        plt.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.15)
        plt.savefig(p_png, format='png', dpi=600, bbox_inches='tight', facecolor='white')
        plt.savefig(p_svg, format='svg', dpi=600, bbox_inches='tight', facecolor='white')
        plt.close()

        logger.info(f"Successfully generated Motif UpSet Plot: {p_png}")
        return p_png

    def run_all(self):
        """Executes UpSet domain and motif plot generation."""
        logger.info("--- Running UpSet Domain & Motif Analysis ---")
        self.plot_domain_upset()
        self.plot_motif_upset()


def main():
    parser = argparse.ArgumentParser(description="TF-NRD UpSet Plot Analysis Suite")
    parser.add_argument("--all", action="store_true", help="Run UpSet analysis for domain and motif datasets")
    parser.add_argument("--domain-upset", action="store_true", help="Generate UpSet plot for PFAM domains")
    parser.add_argument("--motif-upset", action="store_true", help="Generate UpSet plot for InterPro motifs")
    parser.add_argument("--input-dir", default="input_data", help="Directory containing input data")
    parser.add_argument("--output-dir", default="results/Figures", help="Directory to save figures")

    args = parser.parse_args()
    analyzer = UpSetAnalyzer(input_dir=args.input_dir, output_dir=args.output_dir)

    if args.all:
        analyzer.run_all()
    else:
        if args.domain_upset:
            analyzer.plot_domain_upset()
        if args.motif_upset:
            analyzer.plot_motif_upset()
        if not (args.domain_upset or args.motif_upset):
            analyzer.run_all()


if __name__ == "__main__":
    main()
