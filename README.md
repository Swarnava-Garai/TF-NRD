# TF-NRD: Transcription Factors Non-Redundant Dataset

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**TF-NRD (Transcription Factors Non-Redundant Dataset)** is a curated, non-redundant atlas of transcription factor (TF) sequences and multi-entity biological structures designed to support structural, functional, and evolutionary analysis of TF assemblies across all domains of life.

---

## 🔬 Overview

Transcription factors (TFs) are key regulators of gene expression. They form intricate biological assemblies interacting with diverse macromolecular partners, including **DNA**, **RNA**, **DNA–RNA hybrids**, and **protein co-regulators**.

The **TF-NRD pipeline** provides standardized datasets, structural interface properties, domain/motif distributions, and automated visualization suites to facilitate:

- **Macromolecular Interface Characterization**: Quantifying Buried Surface Area (BSA), Fraction of Nonpolar Atoms (FNP), Fraction of Buried Uncharged Residues (FBU), and Linker Density (LD).
- **Subcellular & Taxonomic Profiling**: Categorizing TFs across 6 subcellular locations (`ON`, `OC`, `NC`, `NO`, `CO`, `OO`) and taxonomic domains of life (Eukaryota, Bacteria, Archaea, Viruses).
- **Functional Enrichment & Motif Signature Analysis**: Mapping PFAM domains, InterPro sequence motifs, Logomaker sequence logos, and KEGG human disease pathways.
- **Computational Modeling & Machine Learning**: Providing benchmark datasets for predicting protein–nucleic acid and protein–protein interaction interfaces.

---

## ✨ Key Highlights

- **Comprehensive Multi-Entity Dataset**: Encompasses sequence-curated TFs ($n = 3,570$) and non-redundant structure complexes ($n = 377$).
- **Interface Property Atlas**: Detailed metrics for Protein–Nucleic Acid (PNA) and Protein–Protein (PP) interfaces.
- **Logomaker Sequence Logos**: Automated generation of publication-grade sequence logos for sequence motif signatures.
- **UpSet Intersection Visualizations**: High-resolution UpSet plots displaying domain and motif co-occurrence across subcellular locations.
- **Manuscript Supplementary Data**: Complete Excel workbook and structured JSON exports for all 17 manuscript tables (`Table S1.A` to `Table S15`).

---

## 📁 Directory Structure

```text
TF-NRD/
├── README.md                          # Main project documentation
├── Supplementary/                     # Manuscript supplementary data (Excel & JSON)
│   ├── Supplementary.xlsx             # Master supplementary Excel workbook (17 tables)
│   └── *.json                         # JSON exports for all supplementary tables
├── script/                            # Core Python processing pipelines & visualizers
│   ├── README.md                      # Documentation for script directory
│   ├── blast2nr.py                    # Sequence similarity search & clustering module
│   ├── bsa_oligomer_analysis.py       # BSA & oligomeric state distribution analyzer
│   ├── cross_validate_interface.py    # Interface metric cross-validation module
│   ├── disease_annotation.py          # Disease association processing module
│   ├── final_compare_difference_similarity.py # Interface atom commonality & uniqueness
│   ├── generate_interface_file_path.py# Path generator for interface calculations
│   ├── generate_supplementary_json.py# Exports Supplementary XLSX tables to JSON
│   ├── Interface_calculations.py     # Structural BSA & interaction atomic parser
│   ├── Interface_features.py         # Interface property metric extraction suite
│   ├── kegg_classify.py               # KEGG disease pathway classification script
│   ├── kingdom_domain_of_life_classification.py # Taxonomy domain of life classifier
│   ├── mapper_TF_pdb_Pfam_Rfam.py     # TF PDB to Pfam/Rfam mapping script
│   ├── sequence_curator.py            # Sequence curation suite (SequenceCurator class)
│   ├── sequence_dataset_motif_stat.py # Sequence motif statistics calculator
│   ├── structural_dataset_domain_stat.py # Structural domain statistics calculator
│   ├── subcellular_localization.py    # Subcellular localization detailed classifier
│   ├── tf_feature_analysis.py         # Master feature analysis & plotting suite
│   ├── tf_interface_analysis.py       # Structural interface property analysis suite
│   ├── tf_interface_common_uniqueness.py # Common & unique interface residue analyzer
│   └── upset_analysis.py              # UpSet plot visualization suite (UpSetAnalyzer)
├── utils/                             # Shared metadata curation & PDB utilities
│   ├── mmcif_clean_reader.py          # mmCIF reader & custom report cleaner
│   └── utils.py                       # Master metadata cleaning & complex table builder
├── input_data/                        # Raw input datasets & metadata (Pending verification push)
│   ├── RCSB_PDB/                      # RCSB PDB custom metadata reports
│   ├── UniProtKB/                     # UniProtKB sequence annotations
│   ├── blast/                         # BLAST similarity search databases
│   ├── bsa/                           # PRince interface BSA calculation files (.int)
│   ├── disease/                       # Disease association raw datasets
│   ├── domain_of_life/                # Organism domain of life taxonomy datasets
│   ├── domains/                       # PFAM domain mapping files
│   ├── kegg/                          # KEGG disease pathway annotations
│   ├── motif/                         # PROSITE and InterPro motif datasets
│   ├── sequences/                     # FASTA sequence datasets
│   └── subcellular_location/          # Subcellular location detailed datasets
├── results/                           # Output plots, tables, and statistics (Pending verification push)
│   ├── Figures/                       # Publication figures (600 DPI PNG & SVG in dna/, rna/, all/)
│   ├── Interface/                     # PNA & PP interface property tables and statistical outputs
│   ├── Sequences/                     # Cleaned FASTA files and sequence dataset statistics
│   ├── Structures/                    # Non-redundant structure summaries and complex tables
│   ├── blast/                         # Sequence clustering & blast outputs
│   ├── disease/                       # Disease association summary tables
│   ├── domain_of_life/                # Taxonomy classification outputs
│   ├── domains/                       # PFAM domain distribution tables
│   ├── kegg/                          # KEGG pathway enrichment tables
│   ├── motif/                         # Motif statistics and sequence logos
│   └── subcellular_location/          # Subcellular localization outputs
└── notebook/                          # Exploration & metadata processing Jupyter notebooks
    ├── Parse_supp.ipynb               # Supplementary table parsing notebook
    └── clean_rcsb_metadata.ipynb      # RCSB metadata cleaning workflow notebook
```

---

## 📂 Detailed Folder Descriptions

### 1. `script/`

Contains modular, executable Python scripts powering the TF-NRD pipeline:

- **`sequence_curator.py`**: Static method suite for quality filtering FASTA records (protein length $\ge 30$, nucleic acid $\ge 5$), chain ID normalization, and multi-entity complex table construction.
- **`tf_feature_analysis.py`**: Master visualization suite generating broken y-axis subcellular localization plots (Figure 6), horizontal top-10 PFAM domain bar plots, Logomaker sequence logos, and KEGG pathway bar plots (Figure 9).
- **`tf_interface_analysis.py`**: Analyzes structural interface property metrics (BSA, FNP, FBU, LD) across Protein–Nucleic Acid (PNA) and Protein–Protein (PP) interfaces.
- **`upset_analysis.py`**: Builds publication-quality UpSet plots visualizing domain and motif overlaps across subcellular locations.
- **`generate_supplementary_json.py`**: Automates JSON export for all 17 tables in `Supplementary.xlsx`.

### 2. `utils/`

Shared core modules imported across analysis scripts:

- **`utils.py`**: Custom RCSB PDB report parser, chain ID combiner, oligomeric state builder, and master cleaning pipeline (`clean_and_merge_custom_reports`).
- **`mmcif_clean_reader.py`**: Robust parser for multi-entity mmCIF/PDB custom metadata reports.

### 3. `Supplementary/`

Contains the complete manuscript supplementary dataset:

- **`Supplementary.xlsx`**: Master Excel workbook containing 17 sheets (`Table S1.A` to `Table S15`).
- **`*.json`**: Machine-readable JSON exports corresponding to each sheet in the workbook (`Table S1.A.json` – `Table S15.json`).

### 4. `input_data/` *(To be pushed after verification)*

Contains raw input datasets used for feature extraction and structural analysis:

- **`RCSB_PDB/`**: Raw custom reports from RCSB PDB.
- **`UniProtKB/`**: UniProt annotations and cross-references.
- **`bsa/`**: Raw PRince interface atomic surface area calculation (`.int`) files.
- **`domains/`**, **`motif/`**, **`kegg/`**, **`disease/`**, **`subcellular_location/`**, **`domain_of_life/`**: Domain-specific input tables.

### 5. `results/` *(To be pushed after verification)*

Contains all generated outputs:

- **`Figures/`**: 600 DPI publication figures (PNG and SVG format) organized by category (`dna/`, `rna/`, `all/`).
- **`Interface/`**: Exported multi-sheet Excel files, statistical tests, and summary metrics for PNA and PP interfaces.
- **`Structures/`**, **`Sequences/`**, **`domains/`**, **`kegg/`**, **`motif/`**: Processed summary tables.

### 6. `notebook/`

Jupyter notebooks for interactive exploration, metadata validation, and supplementary table parsing (`Parse_supp.ipynb`, `clean_rcsb_metadata.ipynb`).

---

## 🚀 Quick Start

### Environment Requirements

- Python 3.8+
- `pandas`, `numpy`, `matplotlib`, `logomaker`, `upsetplot`, `openpyxl`

### Execution Examples

```bash
# 1. Run main feature analysis suite (subcellular, PFAM domains, KEGG, domain UpSet, interface)
python script/tf_feature_analysis.py --all

# 2. Run sequence motif analysis explicitly
python script/tf_feature_analysis.py --motifs --motif-logos --motif-upset --category dna

# 3. Execute UpSet domain & motif analysis
python script/upset_analysis.py --all

# 4. Generate JSON exports from Supplementary.xlsx
python script/generate_supplementary_json.py
```

---

## 🌐 Website & Data Availability

The TF-NRD dataset and webserver interface are publicly available at:

🔗 **[http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html](http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html)**

---

## 📖 Reference & Citation

Garai, S., Kant, S., & Bahadur, R. P. (2026).  
**An atlas of non-redundant sequences and structures of transcription factor assemblies across domains of life.**

If you use TF-NRD datasets or scripts in your research, please cite the reference above.

---

## 📜 License

This project is licensed under the [MIT License](LICENSE).
