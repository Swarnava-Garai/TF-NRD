# TF-NRD Scripts Documentation

This directory contains the core modular Python scripts for the **TF-NRD (Transcription Factor Non-Redundant Dataset)** pipeline, including sequence curation, length quality filtering, multi-entity metadata processing, feature visualization, and Logomaker sequence logo generation.

---

## 📁 Directory Structure

```text
TF-NRD/script/
├── README.md                    # Documentation for script directory
├── sequence_curator.py          # Sequence curation & quality filtering suite (SequenceCurator class)
├── tf_feature_analysis.py       # Feature analysis & visualization suite (TFFeatureAnalyzer class)
└── upset_analysis.py            # UpSet plot visualization suite (UpSetAnalyzer class)
```

---

## 1. `sequence_curator.py`

[`sequence_curator.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/sequence_curator.py) provides the **`SequenceCurator`** static method class, designed for bioinformatics dataset curation, sequence quality control, multi-entity RCSB PDB report handling, chain ID aggregation, and complex table summarization.

### Class Architecture (`SequenceCurator`)

```python
from script.sequence_curator import SequenceCurator
```

#### Key `@staticmethod` Functions

- **`SequenceCurator.detect_molecule_type(sequence: str) -> str`**
  Detects molecule type (`'DNA'`, `'RNA'`, or `'PROTEIN'`) based on nucleotide and amino acid character sets.

- **`SequenceCurator.parse_fasta(file_path: Path) -> list`**
  Parses a FASTA file into a list of `[(header, sequence), ...]` tuples.

- **`SequenceCurator.extract_chain_id(header: str) -> str`**
  Extracts chain identifiers from FASTA header strings.

- **`SequenceCurator.filter_fasta_by_length(input_fasta: Path, min_protein_len=30, min_na_len=5) -> dict`**
  Quality checks FASTA records and flags sequence length violations (protein $\le 30$ aa, nucleic acid $\le 5$ nt).

- **`SequenceCurator.read_rcsb_csv(csv_path: Path) -> pd.DataFrame`**
  Intelligently reads RCSB PDB custom report CSV files, detecting 1-header or 2-header formats.

- **`SequenceCurator.normalize_columns(df: pd.DataFrame) -> pd.DataFrame`**
  Normalizes varying custom column headers (`Auth Asym ID` $\rightarrow$ `auth_asym_id`, `Accession Code(s)` $\rightarrow$ `uniprot_id`, etc.).

- **`SequenceCurator.preprocess_df(df: pd.DataFrame) -> pd.DataFrame`**
  Cleans whitespace and forward-fills entry-level metadata across multi-entity PDB rows.

- **`SequenceCurator.merge_metadata_report(seq_df: pd.DataFrame, metadata_csv: Path) -> pd.DataFrame`**
  Merges structure metadata report columns (Resolution, Source Organism, Title) on `entry_id`.

- **`SequenceCurator.filter_complete_complexes(df: pd.DataFrame, min_protein_len=30, min_na_len=5) -> pd.DataFrame`**
  Filters Entry IDs to keep complexes containing **both** a qualifying protein (length $\ge 30$) AND nucleic acid (length $\ge 5$).

- **`SequenceCurator.combine_asym_ids(df: pd.DataFrame, condition_dict=None, chain_sep='') -> pd.DataFrame`**
  Merges unique chain IDs for identical entity sequences within a PDB entry (e.g. `A`, `B` $\rightarrow$ `AB`).

- **`SequenceCurator.build_complex_table(df: pd.DataFrame) -> pd.DataFrame`**
  Generates a 1-row-per-Entry-ID summary table with columns: `['Entry ID', 'Refinement Resolution (Å)', 'Source Organism', 'protein_name', 'protein_length', 'protein_chain', 'DNA_name', 'DNA_length', 'DNA_chain', 'RNA_name', 'RNA_length', 'RNA_chain', 'Chains']`.

- **`SequenceCurator.export_fasta(df, output_file, mol_type='polypeptide(L)', line_width=0)`**
  Exports sequence datasets to FASTA format.

- **`SequenceCurator.export_csv(df, output_file, mol_type='polypeptide(L)')`**
  Exports cleaned datasets to CSV format.

- **`SequenceCurator.export_complex_table(df, output_file)`**
  Exports complex summary tables to CSV format.

### Command-Line Usage

```bash
# 1. Curation & Complex Table Generation
python script/sequence_curator.py --filter-complexes --complex-table-csv output_data/TFNRDv1.0_complex_table.csv

# 2. FASTA Quality Check (Filter length <= 30 for protein, <= 5 for NA)
python script/sequence_curator.py --check-fasta output_data/TFNRDv1.0_Sequence_367_oneline.fasta
```

---

## 2. `tf_feature_analysis.py`

[`tf_feature_analysis.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/tf_feature_analysis.py) provides the **`TFFeatureAnalyzer`** suite for feature analysis, publication-quality plotting (600 DPI), and Logomaker sequence logo generation.

### Class Architecture (`TFFeatureAnalyzer`)

```python
from script.tf_feature_analysis import TFFeatureAnalyzer

analyzer = TFFeatureAnalyzer(input_dir="input_data", plot_dir="results/Figures")
```

#### Key Methods

- **`analyzer.plot_subcellular_location(category='all')`**
  Generates comparison bar plots comparing Sequence TFs vs Structure TFs across 6 subcellular localization categories (`ON`: Only Nucleus, `OC`: Only Cytoplasm, `NC`: Nucleus and Cytoplasm, `NO`: Nucleus and Other, `CO`: Cytoplasm and Other, `OO`: Other Miscellaneous Locations).

- **`analyzer.plot_top_pfam_domains(category='all', top_n=10)`**
  Generates horizontal bar plots of top 10 PFAM structure domains.

- **`analyzer.plot_top_motifs(category='all', top_n=10)`**
  Generates horizontal bar plots of top 10 sequence motifs for DNA, RNA, or combined datasets.

- **`analyzer.plot_motif_sequence_logos(category='all', min_instances=3)`**
  Extracts motif sequence signatures and renders position-frequency information sequence logos using `logomaker` (vector SVG and 600 DPI PNG format).

- **`analyzer.plot_kegg_pathways(category='all', min_hits=10)`**
  Generates horizontal bar plots of top KEGG disease pathways.

- **`analyzer.run_all(include_motifs=False)`**
  Executes feature analysis suite (subcellular, PFAM domains, KEGG pathways, domain UpSet, interface) across `dna/`, `rna/`, and `all/` categories. Sequence motif figures are excluded by default and from `--all` unless `include_motifs=True` or explicit motif flags are passed.

### Target Directory Output Structure

Figures are automatically saved into category subdirectories:
- **DNA Figures**: [`results/Figures/dna/`](file:///home/labuser/Projects/PhD_projects/TF-NRD/results/Figures/dna/) & [`results/Figures/dna/motifs/`](file:///home/labuser/Projects/PhD_projects/TF-NRD/results/Figures/dna/motifs/)
- **RNA Figures**: [`results/Figures/rna/`](file:///home/labuser/Projects/PhD_projects/TF-NRD/results/Figures/rna/) & [`results/Figures/rna/motifs/`](file:///home/labuser/Projects/PhD_projects/TF-NRD/results/Figures/rna/motifs/)
- **Combined Figures**: [`results/Figures/all/`](file:///home/labuser/Projects/PhD_projects/TF-NRD/results/Figures/all/)

### Command-Line Usage

```bash
# Run feature analysis suite (subcellular, PFAM domains, KEGG, domain UpSet, interface) excluding sequence motifs
python script/tf_feature_analysis.py --all

# Sequence motif figures are kept as totally separate arguments:
python script/tf_feature_analysis.py --motifs --motif-logos --motif-upset --category dna
```

---

## 3. `upset_analysis.py`

[`upset_analysis.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/upset_analysis.py) provides the **`UpSetAnalyzer`** suite for visualizing PFAM domain and InterPro motif distributions across subcellular localization categories (`ON`, `OC`, `NC`, `NO`, `CO`, `OO`).

### Class Architecture (`UpSetAnalyzer`)

```python
from script.upset_analysis import UpSetAnalyzer

analyzer = UpSetAnalyzer(input_dir="input_data", output_dir="results/Figures")
```

#### Key Methods

- **`analyzer.plot_domain_upset(category='all')`**
  Maps structure subcellular localization categories with PFAM domain accessions from `standard_start_end_domain.xlsx` and generates UpSet plots (`TFNRD_Domain_UpSet.png` & `TFNRD_Domain_UpSet.svg`).

- **`analyzer.plot_motif_upset(category='all')`**
  Maps sequence subcellular location categories with InterPro motif accessions from `nr_sequence_dataset_motif_details.xlsx` and generates UpSet plots (`TFNRD_Motif_distribution.png` & `TFNRD_Motif_distribution.svg`).

- **`analyzer.run_all()`**
  Executes UpSet domain and motif plotting across `dna/`, `rna/`, and `all/` categories.

### Command-Line Usage

```bash
# Run UpSet domain & motif analysis across all categories
python script/upset_analysis.py --all

# Run specific domain UpSet plot for RNA category
python script/upset_analysis.py --domain-upset --category rna
```

---

## 4. Metadata Curation Suite (`utils.py`)

[`utils.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/utils/utils.py) provides Section 12 RCSB PDB custom report curation, chain combining, entity counting, oligomeric state extraction, BSA calculation, and complex table construction.

### Key Functions in `utils.py`

- **`read_rcsb_custom_report(csv_path: str | Path) -> pd.DataFrame`**
  Intelligently reads RCSB PDB custom report CSV files, automatically detecting single-header vs two-header formats.

- **`combine_asym_ids(df: pd.DataFrame, condition_dict=None, chain_sep='') -> pd.DataFrame`**
  Combines asymmetric chain IDs for identical macromolecule entities within a PDB entry.

- **`filter_complete_complexes(df: pd.DataFrame, min_protein_len=30, min_na_len=5) -> pd.DataFrame`**
  Filters PDB entries to retain complete complexes containing qualifying protein ($\ge 30$ AA) and nucleic acid ($\ge 5$ nt).

- **`abbreviate_organism(name: str) -> str`**
  Abbreviates binomial species names (e.g., `Homo sapiens` $\rightarrow$ `H. sapiens`, `Escherichia coli` $\rightarrow$ `E. coli`).

- **`calculate_bsa_from_int(filepath: str | Path) -> float`**
  Calculates total Buried Surface Area (BSA in $\text{\AA}^2$) from PRince `.int` interface output files.

- **`build_complex_table(df: pd.DataFrame, prince_results_dir=None) -> pd.DataFrame`**
  Constructs a 1-row-per-Entry-ID complex summary table containing:
  - Entity counts: **`Protein Entities`**, **`RNA Entities`**, **`DNA Entities`**
  - Assembly metadata: **`Oligomeric State`**, `Refinement Resolution (\AA)`, `Source Organism`
  - Interface metrics: **`BSA_complex`**, **`BSA_protein`**, **`BSA_DNA`**, **`BSA_RNA`**

- **`clean_and_merge_custom_reports(struct_csv_path, seq_csv_path, prince_results_dir=None) -> (pd.DataFrame, pd.DataFrame)`**
  End-to-end master pipeline reading custom RCSB reports, forward-filling multi-entity fields, merging sequence metadata, combining chain IDs, filtering complexes, and exporting cleaned metadata and complex summary tables.

### Python Usage Example

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path("../utils").resolve()))
from utils import clean_and_merge_custom_reports

struct_csv = "TFNRD_EM_structure_custom_report.csv"
seq_csv = "TFNRD_EM_sequence.csv"

# Execute master cleaning & merging pipeline
merged_df, complex_df = clean_and_merge_custom_reports(struct_csv, seq_csv)

# Export cleaned outputs
merged_df.to_csv("TFNRDv1.0_EM_cleaned_merged_metadata.csv", index=False)
complex_df.to_csv("TFNRDv1.0_EM_complex_table.csv", index=False)
```

---

## 🚀 Environment Requirements

- Python 3.8+
- `pandas`
- `numpy`
- `matplotlib`
- `upsetplot`
- `logomaker`
- `openpyxl`

