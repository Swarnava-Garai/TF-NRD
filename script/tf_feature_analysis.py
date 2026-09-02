"""
tf_feature_analysis.py
----------------------
Transcription Factor Feature Analysis Suite for TF-NRD.
Parses, analyzes, and visualizes Subcellular Localization, PFAM Structure Domains,
Sequence Motifs, Logomaker Sequence Logos, KEGG Disease Pathways, and Organism Data.
Generates publication-quality figures (600 DPI) into structured dna/ and rna/ directories.
"""

from pathlib import Path
import logging
import argparse
import sys
import re
import textwrap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
PROJECT_DIR = SCRIPT_DIR.parent
INPUT_DATA = PROJECT_DIR / "input_data"
RESULTS_DIR = PROJECT_DIR / "results"
PLOT_DIR = RESULTS_DIR / "Figures"
LOG_FILE = SCRIPT_PATH.with_suffix(".log")

# Logging setup
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

# Global plotting style configuration
PLT_STYLE = {
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "figure.dpi": 600
}

# Subcellular Location Mapping
SUBCELLULAR_LOCATION_MAP = {
    "ON": "Only Nucleus",
    "OC": "Only Cytoplasm",
    "NC": "Nucleus and Cytoplasm",
    "NO": "Nucleus and Other (Excluding Cytoplasm)",
    "CO": "Cytoplasm and Other (Excluding Nucleus)",
    "OO": "Other Miscellaneous Locations"
}


class TFFeatureAnalyzer:
    """
    Modular feature analysis suite for Transcription Factor sequence and structure datasets.
    """

    def __init__(self, input_dir: Path = INPUT_DATA, plot_dir: Path = PLOT_DIR):
        self.input_dir = Path(input_dir)
        self.plot_dir = Path(plot_dir)
        self.plot_dir.mkdir(parents=True, exist_ok=True)
        plt.rcParams.update(PLT_STYLE)

    def _get_target_dir(self, category: str, subfolder: str = "") -> Path:
        """
        Helper method to get or create target directory (e.g. plot_dir/dna, plot_dir/rna, plot_dir/all).
        """
        target = self.plot_dir / str(category).lower()
        if subfolder:
            target = target / subfolder
        target.mkdir(parents=True, exist_ok=True)
        return target

    def plot_subcellular_location(self, category: str = "all") -> Path:
        """
        Generates comparison bar plot for Sequence TFs vs Structure TFs across subcellular locations
        with a broken Y-axis to cleanly display high (2158) vs lower (7-552) values matching Figure 6.
        """
        out_dir = self._get_target_dir(category)
        data = {
            'Subcellular_Location': ["ON", "OC", "NC", "NO", "CO", "OO"],
            'Count_Seq': [2158, 118, 552, 158, 64, 52],
            'Count_Str': [107, 37, 52, 15, 7, 7]
        }
        df = pd.DataFrame(data)

        # Create broken y-axis subplot grid
        fig, (ax1, ax2) = plt.subplots(
            2, 1,
            sharex=True,
            figsize=(8.5, 6),
            dpi=600,
            gridspec_kw={'height_ratios': [1, 3.2], 'hspace': 0.08}
        )

        x = np.arange(len(df))
        bar_width = 0.36

        # Colors matching Figure_6.png: Dark Navy Blue for Sequence TFs, Terracotta Red for Structure TFs
        color_seq = "#1F386B"
        color_str = "#C0432E"

        # Plot bars on both axes
        for ax in (ax1, ax2):
            ax.bar(x - bar_width/2, df['Count_Seq'], width=bar_width, label='Sequence TFs', color=color_seq)
            ax.bar(x + bar_width/2, df['Count_Str'], width=bar_width, label='Structure TFs', color=color_str)

        # Set Y-axis limits for broken axis
        ax1.set_ylim(2000, 2200)  # Top panel for 2,158
        ax2.set_ylim(0, 650)      # Bottom panel for 0-600

        # Y-ticks
        ax1.set_yticks([2000, 2100, 2200])
        ax2.set_yticks([0, 100, 200, 300, 400, 500, 600])
        ax1.tick_params(labelsize=13)
        ax2.tick_params(labelsize=13)

        # Hide spines between ax1 and ax2
        ax1.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        ax1.xaxis.tick_top()
        ax1.tick_params(top=False, labeltop=False)  # don't put tick labels at the top
        ax2.xaxis.tick_bottom()

        # Grid lines
        ax1.yaxis.grid(True, linestyle='--', linewidth=0.5, alpha=0.5, color='#CCCCCC')
        ax2.yaxis.grid(True, linestyle='--', linewidth=0.5, alpha=0.5, color='#CCCCCC')

        # Add slant break marks on Y-axis
        d = .015  # diagonal break mark length
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False, linewidth=1.2)
        ax1.plot((-d, +d), (-d * 3.2, +d * 3.2), **kwargs)

        kwargs.update(transform=ax2.transAxes)
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)

        # Add data labels on top of bars
        for i in range(len(df)):
            seq_val = df.loc[i, 'Count_Seq']
            str_val = df.loc[i, 'Count_Str']

            # Label for Sequence TF
            if seq_val > 2000:
                ax1.text(
                    x[i] - bar_width/2, seq_val + 15,
                    f"{seq_val:,}", ha='center', va='bottom', fontsize=13, color='black'
                )
            else:
                ax2.text(
                    x[i] - bar_width/2, seq_val + 10,
                    f"{seq_val:,}", ha='center', va='bottom', fontsize=13, color='black'
                )

            # Label for Structure TF
            ax2.text(
                x[i] + bar_width/2, str_val + 10,
                f"{str_val:,}", ha='center', va='bottom', fontsize=13, color='black'
            )

        # Set X-ticks
        ax2.set_xticks(x)
        ax2.set_xticklabels(df['Subcellular_Location'], fontsize=14)
        ax2.set_xlabel('Subcellular Location', fontweight='bold', fontsize=15, labelpad=8)

        # Y-axis label centered vertically across both subplots
        fig.text(0.02, 0.5, 'Number of Transcription Factors', va='center', rotation='vertical', fontweight='bold', fontsize=15)

        # Add legend in ax1 (top right)
        ax1.legend(frameon=False, fontsize=13.5, loc='upper right')

        # Add Subcellular Location Map legend text as a unified block starting in ax1
        legend_lines = [f"{k}: {v}" for k, v in SUBCELLULAR_LOCATION_MAP.items()]
        legend_text = "\n".join(legend_lines)

        ax1.text(
            0.65, 2185, legend_text,
            ha='left', va='top',
            fontsize=11.5,
            linespacing=1.28,
            clip_on=False
        )

        plt.subplots_adjust(left=0.12, right=0.96, top=0.96, bottom=0.12)

        # Figure 6: Subcellular distribution of TFs
        output_fig = out_dir / "Figure_6.png"
        output_pdf = out_dir / "Figure_6.pdf"
        plt.savefig(output_fig, bbox_inches="tight", dpi=600)
        plt.savefig(output_pdf, bbox_inches="tight", dpi=600)
        plt.savefig(out_dir / "Subcellular_location_bar_plot.png", bbox_inches="tight", dpi=600)

        main_fig = self.plot_dir / "Figure_6.png"
        main_pdf = self.plot_dir / "Figure_6.pdf"
        plt.savefig(main_fig, bbox_inches="tight", dpi=600)
        plt.savefig(main_pdf, bbox_inches="tight", dpi=600)
        plt.close()

        logger.info(f"Successfully generated Subcellular Location Bar Plot Figure 6 ({category.upper()}): {main_fig}")
        return main_fig

    def plot_top_pfam_domains(self, category: str = "all", top_n: int = 10) -> Path:
        """
        Generates horizontal bar plot of top N PFAM structure domains.
        """
        out_dir = self._get_target_dir(category)
        domain_file = self.input_dir / "domains" / "domain_stat_structure_dataset_total.xlsx"
        if not domain_file.exists():
            domain_file = PROJECT_DIR / "results" / "domains" / "domain_stat_structure_dataset_total.xlsx"
        if not domain_file.exists():
            domain_file = PROJECT_DIR / "results" / "domains" / "domain_stat_structure_dataset_TF_NRD_final_377.xlsx"
        if not domain_file.exists():
            domain_file = self.input_dir / "domain" / "domain_stat_structure_dataset_total.xlsx"

        if not domain_file.exists():
            logger.error(f"PFAM Domain file not found in input_data or results: {domain_file}")
            return None

        df = pd.read_excel(domain_file)
        name_col = "PFAM_NAME" if "PFAM_NAME" in df.columns else (df.columns[1] if len(df.columns) > 1 else df.columns[0])
        val_col = "Number_of_PDBs" if "Number_of_PDBs" in df.columns else (df.columns[2] if len(df.columns) > 2 else df.columns[-1])

        df_filtered = df.sort_values(by=val_col, ascending=False).head(top_n)

        labels = [str(label) for label in df_filtered[name_col]]
        values = [int(v) for v in df_filtered[val_col]]

        bar_color = "#1F3B73" if category.lower() == "dna" else ("#C44536" if category.lower() == "rna" else "#4C72B0")

        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=600)
        ax.barh(labels, values, color=bar_color, edgecolor="black", linewidth=0.6)

        offset = max(values) * 0.01
        for i, v in enumerate(values):
            ax.text(v + offset, i, str(v), va='center', ha='left', fontsize=10)

        ax.set_xlim(0, max(values) * 1.15)
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.invert_yaxis()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel("Number of PDB Structures", fontweight='bold', fontsize=14)
        ax.set_ylabel("PFAM Domains", fontweight='bold', fontsize=14)

        plt.tight_layout()
        output_fig = out_dir / "Top10_PFAM_bar_plot.png"
        output_pdf = out_dir / "Top10_PFAM_bar_plot.pdf"
        plt.savefig(output_fig, bbox_inches="tight", dpi=600)
        plt.savefig(output_pdf, bbox_inches="tight", dpi=600)
        plt.close()

        logger.info(f"Successfully generated Top PFAM Bar Plot ({category.upper()}): {output_fig}")
        return output_fig

    def plot_top_motifs(self, category: str = "all", top_n: int = 10) -> Path:
        """
        Generates horizontal bar plot of top N sequence motifs for specific category ('dna', 'rna', or 'all').
        """
        out_dir = self._get_target_dir(category)

        if category.lower() == "dna":
            motif_file = self.input_dir / "motif" / "Motifs_dna_519_dataset.csv"
        elif category.lower() == "rna":
            motif_file = self.input_dir / "motif" / "Motifs_rna_359_dataset.csv"
        else:
            motif_file = PROJECT_DIR / "results" / "motif" / "nr_sequence_dataset_motif_stat.xlsx"
            if not motif_file.exists():
                motif_file = self.input_dir / "motif" / "nr_sequence_dataset_motif_stat.xlsx"

        if not motif_file.exists():
            logger.warning(f"Motif dataset file missing for {category}: {motif_file}")
            motif_file = self.input_dir / "motif" / "nr_sequence_dataset_motif_stat.xlsx"

        if motif_file.suffix == ".csv":
            df_raw = pd.read_csv(motif_file)
            stat = df_raw.groupby('Motif_Name')['UniProt ID'].nunique().reset_index()
            stat.columns = ['Motif_Name', 'Number_of_UniProt_IDs']
            df_filtered = stat.sort_values(by="Number_of_UniProt_IDs", ascending=False).head(top_n)
        else:
            df = pd.read_excel(motif_file)
            df = df.rename(columns={df.columns[0]: "Motif_Name", df.columns[1]: "Number_of_UniProt_IDs"})
            df_filtered = df.sort_values(by="Number_of_UniProt_IDs", ascending=False).head(top_n)

        labels = [textwrap.fill(str(label), 25) for label in df_filtered["Motif_Name"]]
        values = df_filtered["Number_of_UniProt_IDs"].tolist()

        bar_color = "#1F3B73" if category.lower() == "dna" else ("#C44536" if category.lower() == "rna" else "#4C72B0")

        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=600)
        ax.barh(labels, values, color=bar_color, edgecolor="black", linewidth=0.6)

        offset = max(values) * 0.01
        for i, v in enumerate(values):
            ax.text(v + offset, i, str(v), va='center', ha='left', fontsize=10)

        ax.set_xlim(0, max(values) * 1.15)
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.invert_yaxis()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel("Number of UniProt IDs", fontweight='bold', fontsize=14)
        ax.set_ylabel(f"{category.upper()} Motif Name", fontweight='bold', fontsize=14)

        plt.tight_layout()
        output_fig = out_dir / "Top10_Motif_bar_plot.png"
        plt.savefig(output_fig, bbox_inches="tight", dpi=600)
        plt.close()

        logger.info(f"Successfully generated Top Motifs Bar Plot ({category.upper()}): {output_fig}")
        return output_fig

    def plot_motif_sequence_logos(self, category: str = "all", min_instances: int = 3) -> list:
        """
        Extracts motif signatures and generates Logomaker sequence logos into category subdirectories ('dna', 'rna', 'all').
        """
        out_logo_dir = self._get_target_dir(category, subfolder="motifs")

        if category.lower() == "dna":
            motif_csv = self.input_dir / "motif" / "Motifs_dna_519_dataset.csv"
            df_motifs = pd.read_csv(motif_csv) if motif_csv.exists() else pd.DataFrame()
        elif category.lower() == "rna":
            motif_csv = self.input_dir / "motif" / "Motifs_rna_359_dataset.csv"
            df_motifs = pd.read_csv(motif_csv) if motif_csv.exists() else pd.DataFrame()
        else:
            dna_csv = self.input_dir / "motif" / "Motifs_dna_519_dataset.csv"
            rna_csv = self.input_dir / "motif" / "Motifs_rna_359_dataset.csv"
            dfs = []
            if dna_csv.exists():
                dfs.append(pd.read_csv(dna_csv))
            if rna_csv.exists():
                dfs.append(pd.read_csv(rna_csv))
            df_motifs = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

        if df_motifs.empty:
            # Fallback to extracting from fasta & motif_details.xlsx
            fasta_file = self.input_dir / "sequences" / "TF_NRD_Sequence_3570.fasta"
            if not fasta_file.exists():
                fasta_file = PROJECT_DIR / "results" / "Sequences" / "TF_Sequence_dataset_3570.fasta"
            motif_details_file = self.input_dir / "motif" / "nr_sequence_dataset_motif_details.xlsx"
            if not fasta_file.exists() or not motif_details_file.exists():
                logger.error(f"Missing FASTA or motif details file for {category} logos.")
                return []

            seq_dict = {}
            with open(fasta_file, "r", encoding="utf-8", errors="replace") as f:
                cur_header = None
                cur_seq = []
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if cur_header:
                            acc = cur_header.split('|')[1] if '|' in cur_header else cur_header.split()[0]
                            seq_dict[acc] = "".join(cur_seq)
                        cur_header = line
                        cur_seq = []
                    else:
                        cur_seq.append(line)
                if cur_header:
                    acc = cur_header.split('|')[1] if '|' in cur_header else cur_header.split()[0]
                    seq_dict[acc] = "".join(cur_seq)

            df_det = pd.read_excel(motif_details_file)

            records = []
            for idx, row in df_det.iterrows():
                acc = str(row.get('Entry', '')).strip()
                motif_str = str(row.get('Motif', ''))
                if acc not in seq_dict or pd.isna(row.get('Motif')):
                    continue

                full_seq = seq_dict[acc]
                matches = re.findall(r'MOTIF\s+(\d+)\.\.(\d+);\s*/note=\"([^\"]+)\"', motif_str)
                for start, end, note in matches:
                    s_idx = int(start) - 1
                    e_idx = int(end)
                    sub_seq = full_seq[s_idx:e_idx]
                    records.append({
                        'Motif_Name': note,
                        'Motif_signature': sub_seq
                    })
            df_motifs = pd.DataFrame(records)

        if df_motifs.empty:
            logger.warning(f"No motif signatures found for {category}.")
            return []

        generated_logos = []
        for name, group in df_motifs.groupby('Motif_Name'):
            lengths = group['Motif_signature'].dropna().str.len()
            if lengths.empty:
                continue

            common_len = lengths.mode()[0]
            seqs = group[group['Motif_signature'].str.len() == common_len]['Motif_signature'].dropna().tolist()
            seqs = [s.upper() for s in seqs if isinstance(s, str)]
            seqs = [re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', s) for s in seqs]

            if len(seqs) < min_instances:
                continue

            try:
                counts = logomaker.alignment_to_matrix(seqs, to_type='information')
                fig, ax = plt.subplots(figsize=(8, 3), dpi=600)
                logo = logomaker.Logo(counts, ax=ax, shade_below=0.5, fade_below=0.5, color_scheme='chemistry')
                logo.style_spines(visible=False)
                logo.style_spines(spines=['left', 'bottom'], visible=True)

                ax.set_ylabel('Information (bits)', fontsize=12, fontweight='bold')
                ax.set_xlabel('Position', fontsize=12, fontweight='bold')
                ax.set_title(f'{category.upper()} Motif Logo: {name}', fontsize=13, fontweight='bold', pad=10)

                n_pos = counts.shape[0]
                ax.set_xticks(range(n_pos))
                ax.set_xticklabels(range(1, n_pos + 1), fontsize=10, fontweight='bold')
                ax.set_xlim(-0.5, n_pos - 0.5)

                safe_name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
                p_svg = out_logo_dir / f'{safe_name}_logo.svg'
                p_png = out_logo_dir / f'{safe_name}_logo.png'

                plt.tight_layout()
                plt.savefig(p_svg, format='svg', bbox_inches='tight', dpi=600)
                plt.savefig(p_png, format='png', bbox_inches='tight', dpi=600)
                plt.close()

                generated_logos.append(p_png)
                logger.info(f"Generated {category.upper()} motif sequence logo for '{name}': {p_png}")
            except Exception as e:
                logger.warning(f"Error rendering logo for '{name}' ({category}): {e}")

        return generated_logos

    def plot_kegg_pathways(self, category: str = "all", top_n: int = 10) -> Path:
        """
        Generates horizontal bar plot of top 10 KEGG disease pathways (Figure 9).
        Viral carcinogenesis (28) and Neutrophil extracellular trap formation (28) top the distribution.
        """
        out_dir = self._get_target_dir(category)
        kegg_file = self.input_dir / "kegg" / "TF_NRD_KEGG_pathways_combined_final.xlsx"
        if not kegg_file.exists():
            kegg_file = PROJECT_DIR / "results" / "kegg" / "TF_NRD_KEGG_pathways_classified.xlsx"

        if not kegg_file.exists():
            logger.error(f"KEGG pathway file not found: {kegg_file}")
            return None

        df = pd.read_excel(kegg_file)
        df = df.rename(columns={df.columns[0]: "KEGG_ID", df.columns[1]: "Pathway_Name", df.columns[2]: "Number_of_Hits"})
        df_filtered = df.sort_values(by="Number_of_Hits", ascending=False).head(top_n)

        labels = [textwrap.fill(str(label), 32) for label in df_filtered["Pathway_Name"]]
        values = df_filtered["Number_of_Hits"].tolist()

        bar_color = "#1F3B73" if category.lower() == "dna" else ("#C44536" if category.lower() == "rna" else "#4C72B0")

        fig, ax = plt.subplots(figsize=(8.0, 5.5), dpi=600)
        ax.barh(labels, values, color=bar_color, edgecolor="black", linewidth=0.6)

        offset = max(values) * 0.015
        for i, v in enumerate(values):
            ax.text(v + offset, i, str(v), va='center', ha='left', fontsize=11, fontweight='bold')

        ax.set_xlim(0, max(values) * 1.15)
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y', labelsize=11)
        ax.invert_yaxis()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel("Number of Transcription Factors", fontweight='bold', fontsize=13)
        ax.set_ylabel("Disease Pathway", fontweight='bold', fontsize=13)

        plt.subplots_adjust(left=0.45)
        plt.tight_layout()

        # Figure 9: Functional enrichment of TFs in human disease pathways
        output_fig = out_dir / "Figure_9.png"
        output_pdf = out_dir / "Figure_9.pdf"
        plt.savefig(output_fig, bbox_inches="tight", dpi=600)
        plt.savefig(output_pdf, bbox_inches="tight", dpi=600)
        plt.savefig(out_dir / "KEGG_pathway_bar_plot.png", bbox_inches="tight", dpi=600)

        main_fig = self.plot_dir / "Figure_9.png"
        main_pdf = self.plot_dir / "Figure_9.pdf"
        plt.savefig(main_fig, bbox_inches="tight", dpi=600)
        plt.savefig(main_pdf, bbox_inches="tight", dpi=600)
        plt.close()

        logger.info(f"Successfully generated KEGG Pathway Bar Plot Figure 9 ({category.upper()}): {main_fig}")
        return main_fig

    def plot_domain_upset(self, category: str = "all") -> Path:
        """Generates UpSet plot for PFAM domain distribution across subcellular locations."""
        try:
            from script.upset_analysis import UpSetAnalyzer
        except ImportError:
            from upset_analysis import UpSetAnalyzer
        upset_analyzer = UpSetAnalyzer(input_dir=self.input_dir, output_dir=self.plot_dir)
        return upset_analyzer.plot_domain_upset(category=category)

    def plot_motif_upset(self, category: str = "all") -> Path:
        """Generates UpSet plot for motif distribution across subcellular locations."""
        try:
            from script.upset_analysis import UpSetAnalyzer
        except ImportError:
            from upset_analysis import UpSetAnalyzer
        upset_analyzer = UpSetAnalyzer(input_dir=self.input_dir, output_dir=self.plot_dir)
        return upset_analyzer.plot_motif_upset(category=category)

    def run_interface_analysis(self) -> dict:
        """
        Executes structural interface property analysis (BSA, FNP, FBU, LD)
        and exports multi-sheet Excel workbook for PNA and PP datasets.
        """
        try:
            from script.tf_interface_analysis import TFInterfaceAnalyzer
        except ImportError:
            from tf_interface_analysis import TFInterfaceAnalyzer
        out_dir = self.plot_dir.parent / "Interface"
        interface_analyzer = TFInterfaceAnalyzer(output_dir=out_dir)
        return interface_analyzer.run_all()

    def run_all(self, include_motifs: bool = False):
        """
        Runs feature and interface analyses and generates figures into dna/, rna/, and all/ directories.
        Sequence motif figures are excluded by default and from --all unless include_motifs=True.
        """
        logger.info("Executing TF-NRD Feature Analysis Pipeline for DNA, RNA, and Combined datasets...")

        for category in ["dna", "rna", "all"]:
            logger.info(f"Processing Category: {category.upper()}")
            self.plot_subcellular_location(category=category)
            self.plot_top_pfam_domains(category=category)
            self.plot_kegg_pathways(category=category)
            self.plot_domain_upset(category=category)
            if include_motifs:
                self.plot_top_motifs(category=category)
                self.plot_motif_sequence_logos(category=category)
                self.plot_motif_upset(category=category)

        logger.info("Running Interface Analysis Pipeline...")
        self.run_interface_analysis()

        logger.info("TF-NRD Feature Analysis & Interface Pipeline Completed Successfully for DNA, RNA, and Combined Datasets!")


def main():
    parser = argparse.ArgumentParser(description="TF-NRD Feature Analysis & Visualization Suite.")
    parser.add_argument("--all", action="store_true", help="Run feature analyses (subcellular, PFAM domains, KEGG, domain UpSet, interface) excluding sequence motifs into dna, rna, and all dirs")
    parser.add_argument("--subcellular", action="store_true", help="Generate Subcellular Localization plot")
    parser.add_argument("--domains", action="store_true", help="Generate Top PFAM Domains plot")
    parser.add_argument("--motifs", action="store_true", help="Generate Top Motifs plot (sequence motif, separate argument)")
    parser.add_argument("--motif-logos", action="store_true", help="Generate Logomaker Motif Sequence Logos (sequence motif, separate argument)")
    parser.add_argument("--kegg", action="store_true", help="Generate KEGG Disease Pathways plot")
    parser.add_argument("--domain-upset", action="store_true", help="Generate PFAM Domain UpSet plot")
    parser.add_argument("--motif-upset", action="store_true", help="Generate Motif UpSet plot (sequence motif, separate argument)")
    parser.add_argument("--interface", action="store_true", help="Run Interface Analysis and export multi-sheet Excel for PNA & PP")
    parser.add_argument("--category", choices=["dna", "rna", "all"], default="all", help="Dataset category focus ('dna', 'rna', or 'all')")
    parser.add_argument("-i", "--input-dir", default=None, help="Input data directory path")
    parser.add_argument("-o", "--output-dir", default=None, help="Output plot directory path")

    args = parser.parse_args()

    input_d = Path(args.input_dir) if args.input_dir else INPUT_DATA
    plot_d = Path(args.output_dir) if args.output_dir else PLOT_DIR

    analyzer = TFFeatureAnalyzer(input_dir=input_d, plot_dir=plot_d)

    if args.all:
        analyzer.run_all(include_motifs=False)
        cat = args.category
        if args.motifs:
            for c in (["dna", "rna", "all"] if cat == "all" else [cat]):
                analyzer.plot_top_motifs(category=c)
        if args.motif_logos:
            for c in (["dna", "rna", "all"] if cat == "all" else [cat]):
                analyzer.plot_motif_sequence_logos(category=c)
        if args.motif_upset:
            for c in (["dna", "rna", "all"] if cat == "all" else [cat]):
                analyzer.plot_motif_upset(category=c)
    else:
        cat = args.category
        ran_any = False
        if args.subcellular:
            analyzer.plot_subcellular_location(category=cat)
            ran_any = True
        if args.domains:
            analyzer.plot_top_pfam_domains(category=cat)
            ran_any = True
        if args.kegg:
            analyzer.plot_kegg_pathways(category=cat)
            ran_any = True
        if args.domain_upset:
            analyzer.plot_domain_upset(category=cat)
            ran_any = True
        if args.interface:
            analyzer.run_interface_analysis()
            ran_any = True

        # Sequence motif options are separate arguments
        if args.motifs:
            analyzer.plot_top_motifs(category=cat)
            ran_any = True
        if args.motif_logos:
            analyzer.plot_motif_sequence_logos(category=cat)
            ran_any = True
        if args.motif_upset:
            analyzer.plot_motif_upset(category=cat)
            ran_any = True

        if not ran_any:
            analyzer.run_all(include_motifs=False)


if __name__ == "__main__":
    main()



