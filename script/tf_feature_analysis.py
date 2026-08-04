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

# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('tf_feature_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("TFFeatureAnalyzer")

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent.resolve()
INPUT_DATA = PROJECT_DIR / "input_data"
RESULTS_DIR = PROJECT_DIR / "results"
PLOT_DIR = RESULTS_DIR / "Figures"

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
        Generates comparison bar plot for Sequence TFs vs Structure TFs across subcellular locations.
        """
        out_dir = self._get_target_dir(category)
        data = {
            'Subcellular_Location': ["ON", "OC", "NC", "NO", "CO", "OO"],
            'Count_Seq': [2158, 118, 552, 158, 64, 52],
            'Count_Str': [107, 37, 52, 15, 7, 7]
        }
        df = pd.DataFrame(data)

        fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
        x = np.arange(len(df))
        bar_width = 0.38

        bar_color_seq = "#1F3B73" if category.lower() != "rna" else "#C44536"
        bars_seq = ax.bar(x - bar_width/2, df['Count_Seq'], width=bar_width, label='Sequence TFs', color=bar_color_seq)
        bars_str = ax.bar(x + bar_width/2, df['Count_Str'], width=bar_width, label='Structure TFs', color="#2B4C7E")

        for bars in [bars_seq, bars_str]:
            for bar in bars:
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width()/2,
                    height + max(df['Count_Seq']) * 0.01,
                    f"{int(height):,}",
                    ha='center',
                    va='bottom',
                    fontsize=12
                )

        ax.set_xticks(x)
        ax.set_xticklabels(df['Subcellular_Location'], fontsize=14)
        ax.set_ylabel('Number of Transcription Factors', fontweight='bold', fontsize=14)
        ax.set_xlabel('Subcellular Location', fontweight='bold', fontsize=14)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.grid(True, linestyle='--', linewidth=0.6, alpha=0.6)
        ax.legend(frameon=False, fontsize=14)

        legend_text = "\n".join([f"{k}: {v}" for k, v in SUBCELLULAR_LOCATION_MAP.items()])
        plt.figtext(0.52, 0.75, legend_text, ha='left', va='top', fontsize=10, linespacing=1.4)

        plt.tight_layout()
        output_fig = out_dir / "Subcellular_location_bar_plot.png"
        plt.savefig(output_fig, bbox_inches="tight", dpi=600)
        plt.close()

        logger.info(f"Successfully generated Subcellular Location Bar Plot ({category.upper()}): {output_fig}")
        return output_fig

    def plot_top_pfam_domains(self, category: str = "all", top_n: int = 10) -> Path:
        """
        Generates horizontal bar plot of top N PFAM structure domains.
        """
        out_dir = self._get_target_dir(category)
        domain_file = self.input_dir / "domain" / "domain_stat_structure_dataset_total.xlsx"
        if not domain_file.exists():
            logger.error(f"PFAM Domain file not found: {domain_file}")
            return None

        df = pd.read_excel(domain_file)
        df = df.rename(columns={df.columns[0]: "PFAM_NAME", df.columns[2]: "Number_of_PDBs"})
        df_filtered = df.sort_values(by="Number_of_PDBs", ascending=True).tail(top_n)

        labels = df_filtered["PFAM_NAME"]
        values = df_filtered["Number_of_PDBs"].tolist()

        bar_color = "#1F3B73" if category.lower() == "dna" else ("#C44536" if category.lower() == "rna" else "#4C72B0")

        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=600)
        ax.barh(labels, values, color=bar_color, edgecolor="black", linewidth=0.6)

        offset = max(values) * 0.01
        for i, v in enumerate(values):
            ax.text(v + offset, i, str(v), va='center', ha='left', fontsize=10)

        ax.set_xlim(0, max(values) * 1.15)
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel("Number of PDB Structures", fontweight='bold', fontsize=14)
        ax.set_ylabel("PFAM Domains", fontweight='bold', fontsize=14)

        plt.tight_layout()
        output_fig = out_dir / "Top10_PFAM_bar_plot.png"
        plt.savefig(output_fig, bbox_inches="tight", dpi=600)
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
        elif category.lower() == "rna":
            motif_csv = self.input_dir / "motif" / "Motifs_rna_359_dataset.csv"
        else:
            motif_csv = None

        if motif_csv and motif_csv.exists():
            df_motifs = pd.read_csv(motif_csv)
        else:
            # Fallback to extracting from fasta & motif_details.xlsx
            fasta_file = self.input_dir / "sequences" / "nr_final_sequence_dataset_3570.fasta"
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

    def plot_kegg_pathways(self, category: str = "all", min_hits: int = 10) -> Path:
        """
        Generates horizontal bar plot of top KEGG disease pathways.
        """
        out_dir = self._get_target_dir(category)
        kegg_file = self.input_dir / "kegg" / "TF_NRD_KEGG_pathways_combined_final.xlsx"
        if not kegg_file.exists():
            logger.error(f"KEGG pathway file not found: {kegg_file}")
            return None

        df = pd.read_excel(kegg_file)
        df = df.rename(columns={df.columns[0]: "KEGG_ID", df.columns[1]: "Pathway_Name", df.columns[2]: "Number_of_Hits"})
        df_filtered = df[df["Number_of_Hits"] > min_hits].sort_values(by="Number_of_Hits", ascending=False)

        labels = [textwrap.fill(str(label), 25) for label in df_filtered["Pathway_Name"]]
        values = df_filtered["Number_of_Hits"].tolist()

        bar_color = "#1F3B73" if category.lower() == "dna" else ("#C44536" if category.lower() == "rna" else "#4C72B0")

        fig, ax = plt.subplots(figsize=(7, 5), dpi=600)
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
        ax.set_xlabel("Number of Transcription Factors", fontweight='bold', fontsize=14)
        ax.set_ylabel("Disease Pathway", fontweight='bold', fontsize=14)

        plt.subplots_adjust(left=0.45)
        plt.tight_layout()
        output_fig = out_dir / "KEGG_pathway_bar_plot.png"
        plt.savefig(output_fig, bbox_inches="tight", dpi=600)
        plt.close()

        logger.info(f"Successfully generated KEGG Pathway Bar Plot ({category.upper()}): {output_fig}")
        return output_fig

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

    def run_all(self):
        """
        Runs all feature and interface analyses and generates figures into dna/, rna/, and all/ directories.
        """
        logger.info("Executing Complete TF-NRD Feature Analysis Pipeline for DNA, RNA, and Combined datasets...")

        for category in ["dna", "rna", "all"]:
            logger.info(f"Processing Category: {category.upper()}")
            self.plot_subcellular_location(category=category)
            self.plot_top_pfam_domains(category=category)
            self.plot_top_motifs(category=category)
            self.plot_motif_sequence_logos(category=category)
            self.plot_kegg_pathways(category=category)
            self.plot_domain_upset(category=category)
            self.plot_motif_upset(category=category)

        logger.info("Running Interface Analysis Pipeline...")
        self.run_interface_analysis()

        logger.info("TF-NRD Feature Analysis & Interface Pipeline Completed Successfully for DNA, RNA, and Combined Datasets!")


def main():
    parser = argparse.ArgumentParser(description="TF-NRD Feature Analysis & Visualization Suite.")
    parser.add_argument("--all", action="store_true", help="Run all feature analyses and generate figures into dna, rna, and all dirs")
    parser.add_argument("--subcellular", action="store_true", help="Generate Subcellular Localization plot")
    parser.add_argument("--domains", action="store_true", help="Generate Top PFAM Domains plot")
    parser.add_argument("--motifs", action="store_true", help="Generate Top Motifs plot")
    parser.add_argument("--motif-logos", action="store_true", help="Generate Logomaker Motif Sequence Logos")
    parser.add_argument("--kegg", action="store_true", help="Generate KEGG Disease Pathways plot")
    parser.add_argument("--domain-upset", action="store_true", help="Generate PFAM Domain UpSet plot")
    parser.add_argument("--motif-upset", action="store_true", help="Generate Motif UpSet plot")
    parser.add_argument("--interface", action="store_true", help="Run Interface Analysis and export multi-sheet Excel for PNA & PP")
    parser.add_argument("--category", choices=["dna", "rna", "all"], default="all", help="Dataset category focus ('dna', 'rna', or 'all')")
    parser.add_argument("-i", "--input-dir", default=None, help="Input data directory path")
    parser.add_argument("-o", "--output-dir", default=None, help="Output plot directory path")

    args = parser.parse_args()

    input_d = Path(args.input_dir) if args.input_dir else INPUT_DATA
    plot_d = Path(args.output_dir) if args.output_dir else PLOT_DIR

    analyzer = TFFeatureAnalyzer(input_dir=input_d, plot_dir=plot_d)

    if args.all:
        analyzer.run_all()
    else:
        cat = args.category
        if args.subcellular:
            analyzer.plot_subcellular_location(category=cat)
        if args.domains:
            analyzer.plot_top_pfam_domains(category=cat)
        if args.motifs:
            analyzer.plot_top_motifs(category=cat)
        if args.motif_logos:
            analyzer.plot_motif_sequence_logos(category=cat)
        if args.kegg:
            analyzer.plot_kegg_pathways(category=cat)
        if args.domain_upset:
            analyzer.plot_domain_upset(category=cat)
        if args.motif_upset:
            analyzer.plot_motif_upset(category=cat)
        if args.interface:
            analyzer.run_interface_analysis()
        if not any([args.subcellular, args.domains, args.motifs, args.motif_logos, args.kegg, args.domain_upset, args.motif_upset, args.interface]):
            analyzer.run_all()


if __name__ == "__main__":
    main()


