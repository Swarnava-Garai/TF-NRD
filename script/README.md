# TF-NRD Scripts Documentation

This directory contains the core modular Python scripts for the **TF-NRD (Transcription Factor Non-Redundant Dataset)** pipeline, including sequence curation, length quality filtering, multi-entity metadata processing, structural interface property analysis, feature visualization, Logomaker sequence logo generation, and UpSet intersection plotting.

---

## 📁 Directory Structure (Alphabetical Order)

```text
TF-NRD/script/
├── README.md                              # Documentation for script directory
├── blast2nr.py                            # Sequence similarity search & BLAST NR clustering module
├── bsa_oligomer_analysis.py               # BSA & oligomeric state distribution analyzer
├── cross_validate_interface.py            # Interface metric cross-validation module
├── disease_annotation.py                  # Disease association processing module
├── final_compare_difference_similarity.py # Interface commonality & uniqueness comparison script
├── generate_interface_file_path.py        # Path generator for interface calculation files
├── generate_supplementary_json.py         # Converts Supplementary XLSX tables to JSON
├── Interface_calculations.py             # Atomic BSA calculation & interaction parser
├── Interface_features.py                 # Interface property metric extraction suite
├── kegg_classify.py                       # KEGG disease pathway classification script
├── kingdom_domain_of_life_classification.py # Taxonomy domain of life classifier
├── mapper_TF_pdb_Pfam_Rfam.py             # TF PDB to Pfam/Rfam mapping module
├── pdb_generate_chain_details.py          # PDB chain & interaction distance KDTree analyzer
├── sequence_curator.py                    # Sequence curation suite (SequenceCurator class)
├── sequence_dataset_motif_stat.py         # Sequence motif statistics calculator
├── structural_dataset_domain_stat.py       # Structural domain statistics calculator
├── subcellular_localization.py            # Subcellular localization detailed classifier
├── tf_feature_analysis.py                 # Feature analysis & visualization suite (TFFeatureAnalyzer)
├── tf_interface_analysis.py               # Structural interface property suite (TFInterfaceAnalyzer)
├── tf_interface_common_uniqueness.py       # Common & unique interface atom analyzer
└── upset_analysis.py                      # UpSet plot visualization suite (UpSetAnalyzer)
```

---

## 📖 Alphabetical Script Reference

---

### 1. `blast2nr.py`

[`blast2nr.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/blast2nr.py) performs sequence similarity searches and non-redundant clustering against NCBI reference databases.

- **Key Functions**:
  - `run_blastp()`: Executes local BLAST searching with sequence identity and $e$-value thresholds.
  - `parse_blast_results()`: Parses BLAST XML/TSV output files to extract identity percentages and coverage metrics.
  - `cluster_sequences()`: Groups sequence hits into non-redundant similarity clusters.
- **CLI Usage**:
  ```bash
  python script/blast2nr.py -i input_data/sequences/ -o results/blast/
  ```

---

### 2. `bsa_oligomer_analysis.py`

[`bsa_oligomer_analysis.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/bsa_oligomer_analysis.py) provides the **`BSAOligomerAnalyzer`** class for analyzing Buried Surface Area (BSA) distributions and classifying oligomeric states across TF complexes.

- **Key Functions & Classes**:
  - `BSAOligomerAnalyzer`: Core analyzer for loading BSA tables and performing statistical calculations.
  - `classify_oligo(val: str)`: Categorizes complex oligomeric states (monomer, homodimer, heterodimer, tetramer, etc.).
  - `analyze_column_and_outliers()`: Identifies BSA distribution outliers using IQR bounds ($1.5 \times \text{IQR}$).
- **CLI Usage**:
  ```bash
  python script/bsa_oligomer_analysis.py --dataset protein_nucleic_acid
  ```

---

### 3. `cross_validate_interface.py`

[`cross_validate_interface.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/cross_validate_interface.py) cross-validates PRince interface calculation (`.int`) files across dataset iterations.

- **Key Functions**:
  - `discover_interface_groups()`: Groups `.int` calculation files by PDB ID and interface type.
  - `check_insertion_codes()`: Validates PDB insertion codes in interface residue definitions.
  - `check_residue_numbering()`: Checks residue alignment and computes residue index offsets across run groups.
- **CLI Usage**:
  ```bash
  python script/cross_validate_interface.py -i input_data/bsa/
  ```

---

### 4. `disease_annotation.py`

[`disease_annotation.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/disease_annotation.py) processes human disease annotations associated with transcription factors.

- **Key Functions**:
  - `extract_clean_disease_names()`: Parses free-text disease entries to extract standardized disease nomenclature.
  - `process_disease_annotations()`: Cross-references UniProt disease annotations with OMIM and KEGG identifiers.
- **CLI Usage**:
  ```bash
  python script/disease_annotation.py -i input_data/disease/ -o results/disease/
  ```

---

### 5. `final_compare_difference_similarity.py`

[`final_compare_difference_similarity.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/final_compare_difference_similarity.py) performs atomic-level pairwise comparison between two interface calculation files.

- **Key Functions**:
  - `parse_atoms()`: Extracts atomic coordinate records from `.int` interface files.
  - `compare_int_files()`: Calculates total common atoms, atoms unique to File 1, and atoms unique to File 2.
- **CLI Usage**:
  ```bash
  python script/final_compare_difference_similarity.py file1.int file2.int
  ```

---

### 6. `generate_interface_file_path.py`

[`generate_interface_file_path.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/generate_interface_file_path.py) generates structured path manifest text files mapping PDB IDs to their corresponding `.int` interface calculation files.

- **Key Functions**:
  - `generate_int_paths()`: Scans target directories and outputs a plain-text path reference file.
- **CLI Usage**:
  ```bash
  python script/generate_interface_file_path.py
  ```

---

### 7. `generate_supplementary_json.py`

[`generate_supplementary_json.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/generate_supplementary_json.py) converts all 17 tables in the manuscript Excel workbook ([`Supplementary.xlsx`](file:///home/labuser/Projects/PhD_projects/TF-NRD/Supplementary/Supplementary.xlsx)) into machine-readable JSON files in [`Supplementary/`](file:///home/labuser/Projects/PhD_projects/TF-NRD/Supplementary/).

- **Key Functions**:
  - `convert_excel_sheets_to_json()`: Handles multi-level headers, title rows, and column deduplication for sheets `Table S1.A` to `Table S15`.
- **CLI Usage**:
  ```bash
  python script/generate_supplementary_json.py
  ```

---

### 8. `Interface_calculations.py`

[`Interface_calculations.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/Interface_calculations.py) calculates solvent-accessible surface area (ASA) differences and atomic interaction contact lines.

- **Key Functions**:
  - `fetch_atomline()`: Extracts PDB atomic coordinates for specified protein and nucleic acid chain groups.
  - `generate_interface_atomfile()`: Computes $\Delta\text{ASA}$ surface area loss upon complex formation using a distance cutoff ($6.0\text{ \AA}$).
  - `calc_interface_area()`: Integrates surface area metrics to determine total Buried Surface Area (BSA).
- **CLI Usage**:
  ```bash
  python script/Interface_calculations.py -i input_data/bsa/
  ```

---

### 9. `Interface_features.py`

[`Interface_features.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/Interface_features.py) serves as a multi-threaded batch runner for PRince interface feature extraction across structural datasets.

- **Key Functions**:
  - `run_prince_protein_nucleic_acid()`: Executes PRince calculation for Protein–Nucleic Acid (PNA) complexes.
  - `run_prince_protein_protein()`: Executes PRince calculation for Protein–Protein (PP) complexes.
  - `execute_batch_mode()`: Parallelizes batch interface processing using `ProcessPoolExecutor`.
- **CLI Usage**:
  ```bash
  python script/Interface_features.py -i input_data/bsa/ -o results/Interface/
  ```

---

### 10. `kegg_classify.py`

[`kegg_classify.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/kegg_classify.py) categorizes TF disease associations into human disease pathway hierarchies using the KEGG REST API and BRITE database.

- **Key Functions**:
  - `fetch_and_parse_brite_hierarchy()`: Downloads and parses KEGG BRITE hierarchy mappings.
  - `get_disease_category_online()`: Queries online KEGG records for specific pathway IDs (`hsa05203`, etc.).
  - `process_kegg_classification()`: Generates consolidated KEGG disease pathway summary tables.
- **CLI Usage**:
  ```bash
  python script/kegg_classify.py -i input_data/kegg/ -o results/kegg/
  ```

---

### 11. `kingdom_domain_of_life_classification.py`

[`kingdom_domain_of_life_classification.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/kingdom_domain_of_life_classification.py) classifies TF source organisms into biological domains of life.

- **Key Functions**:
  - `extract_domain_from_lineage()`: Maps organism lineage strings to 4 primary domain categories (**Eukaryota**, **Bacteria**, **Archaea**, **Viruses**).
  - `classify_domains_of_life()`: Processes UniProt taxonomy metadata to export organism domain breakdown tables.
- **CLI Usage**:
  ```bash
  python script/kingdom_domain_of_life_classification.py -i input_data/domain_of_life/ -o results/domain_of_life/
  ```

---

### 12. `mapper_TF_pdb_Pfam_Rfam.py`

[`mapper_TF_pdb_Pfam_Rfam.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/mapper_TF_pdb_Pfam_Rfam.py) maps PDB structure chain identifiers to Pfam domain accessions (`PF00170`, etc.) and Rfam RNA family IDs.

- **Key Functions**:
  - `map_pdb_to_pfam()`: Maps PDB protein chains to Pfam domain names and residue start/end boundaries.
  - `map_pdb_to_rfam()`: Maps nucleic acid chains to Rfam family annotations.
- **CLI Usage**:
  ```bash
  python script/mapper_TF_pdb_Pfam_Rfam.py -i input_data/domains/ -o results/domains/
  ```

---

### 13. `pdb_generate_chain_details.py`

[`pdb_generate_chain_details.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/pdb_generate_chain_details.py) parses mmCIF structure files to determine physical chain interactions and experimental methods.

- **Key Functions**:
  - `extract_exp_method()`: Parses experimental resolution and method (X-ray, Cryo-EM, NMR) from mmCIF header tags.
  - `analyze_cif_chains_kdtree()`: Uses `scipy.spatial.KDTree` spatial indexing ($5.0\text{ \AA}$ cutoff) to identify contacting protein, DNA, and RNA chain pairs.
- **CLI Usage**:
  ```bash
  python script/pdb_generate_chain_details.py -i input_data/RCSB_PDB/
  ```

---

### 14. `sequence_curator.py`

[`sequence_curator.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/sequence_curator.py) provides the **`SequenceCurator`** static method class for dataset sequence quality control and metadata normalization.

- **Key `@staticmethod` Functions**:
  - `SequenceCurator.detect_molecule_type()`: Classifies sequence records (`DNA`, `RNA`, `PROTEIN`) based on character set compositions.
  - `SequenceCurator.filter_fasta_by_length()`: Enforces length thresholds (protein $\ge 30$ aa, nucleic acid $\ge 5$ nt).
  - `SequenceCurator.parse_fasta()`: Reads FASTA records into `(header, sequence)` tuples.
  - `SequenceCurator.read_rcsb_csv()`: Normalizes RCSB multi-entity metadata reports.

---

### 15. `sequence_dataset_motif_stat.py`

[`sequence_dataset_motif_stat.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/sequence_dataset_motif_stat.py) calculates sequence motif frequency statistics across sequence-curated TF datasets.

- **Key Functions**:
  - `extract_all_motifs()`: Extracts InterPro and PROSITE motif signatures from UniProt annotations.
  - `process_motif_statistics()`: Computes total motif occurrence counts, sequence coverage, and output statistics.
- **CLI Usage**:
  ```bash
  python script/sequence_dataset_motif_stat.py -i input_data/motif/ -o results/motif/
  ```

---

### 16. `structural_dataset_domain_stat.py`

[`structural_dataset_domain_stat.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/structural_dataset_domain_stat.py) calculates PFAM domain distribution statistics across non-redundant structural TF datasets.

- **Key Functions**:
  - `process_structural_domain_stats()`: Aggregates domain frequencies, unique PDB counts per domain, and exports structural domain summary tables.
- **CLI Usage**:
  ```bash
  python script/structural_dataset_domain_stat.py -i input_data/domains/ -o results/domains/
  ```

---

### 17. `subcellular_localization.py`

[`subcellular_localization.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/subcellular_localization.py) categorizes TF subcellular localization annotations into 6 detailed categories.

- **Key Functions**:
  - `clean_subcellular_text()`: Standardizes free-text UniProt subcellular location string entries.
  - `classify_subcellular_category()`: Assigns TFs to subcellular categories:
    - **`ON`**: Only Nucleus
    - **`OC`**: Only Cytoplasm
    - **`NC`**: Nucleus and Cytoplasm
    - **`NO`**: Nucleus and Other
    - **`CO`**: Cytoplasm and Other
    - **`OO`**: Only Other locations
- **CLI Usage**:
  ```bash
  python script/subcellular_localization.py -i input_data/subcellular_location/ -o results/subcellular_location/
  ```

---

### 18. `tf_feature_analysis.py`

[`tf_feature_analysis.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/tf_feature_analysis.py) provides the **`TFFeatureAnalyzer`** master class for generating feature plots and publication figures.

- **Key Methods**:
  - `plot_subcellular_location()`: Generates broken y-axis bar plots for subcellular locations (Figure 6).
  - `plot_top_pfam_domains()`: Renders top-10 PFAM domain horizontal bar plots.
  - `plot_top_motifs()`: Renders top-10 sequence motif horizontal bar plots (*separate argument*).
  - `plot_motif_sequence_logos()`: Generates Logomaker sequence logos (*separate argument*).
  - `plot_kegg_pathways()`: Generates top-10 KEGG disease pathway bar plots (Figure 9).
  - `run_all(include_motifs=False)`: Orchestrates complete feature analysis pipeline excluding sequence motifs by default and from `--all`.
- **CLI Usage**:
  ```bash
  # Run feature analyses (excludes sequence motifs)
  python script/tf_feature_analysis.py --all

  # Run sequence motifs explicitly
  python script/tf_feature_analysis.py --motifs --motif-logos --motif-upset --category dna
  ```

---

### 19. `tf_interface_analysis.py`

[`tf_interface_analysis.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/tf_interface_analysis.py) provides the **`TFInterfaceAnalyzer`** class for extracting structural interface property metrics.

- **Key Functions**:
  - `calculate_interface_metrics()`: Calculates interface parameters: Buried Surface Area (BSA), Fraction of Nonpolar Atoms (FNP), Fraction of Buried Uncharged residues (FBU), and Linker Density (LD).
  - `determine_interface_formed()`: Identifies interface macromolecular types (Protein–DNA, Protein–RNA, Protein–Protein).
  - `TFInterfaceAnalyzer.run_all()`: Exports multi-sheet Excel workbooks and statistical summaries to `results/Interface/`.
- **CLI Usage**:
  ```bash
  python script/tf_interface_analysis.py
  ```

---

### 20. `tf_interface_common_uniqueness.py`

[`tf_interface_common_uniqueness.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/tf_interface_common_uniqueness.py) identifies common and unique interacting interface residue atoms across structural assemblies.

- **Key Functions**:
  - `auto_discover_int_pairs()`: Automatically finds pairs of `.int` calculation files for the same PDB structure.
  - `compare_int_files()`: Identifies atomic coordinate overlaps and unique atom sets.
  - `save_common_atoms()`: Writes common interface residue output tables.
- **CLI Usage**:
  ```bash
  python script/tf_interface_common_uniqueness.py -i input_data/bsa/
  ```

---

### 21. `upset_analysis.py`

[`upset_analysis.py`](file:///home/labuser/Projects/PhD_projects/TF-NRD/script/upset_analysis.py) provides the **`UpSetAnalyzer`** class for generating publication-quality UpSet plots.

- **Key Methods**:
  - `UpSetAnalyzer.plot_domain_upset()`: Renders UpSet plot of PFAM domain distributions across subcellular categories (`ON`, `OC`, `NC`, `NO`, `CO`, `OO`).
  - `UpSetAnalyzer.plot_motif_upset()`: Renders UpSet plot of sequence motif distributions across subcellular categories.
  - `_safe_plot_matrix()` & `_safe_label_sizes()`: Internal monkey-patches resolving Pandas 2.2+ Copy-on-Write and Matplotlib scalar label formatting compatibility.
- **CLI Usage**:
  ```bash
  python script/upset_analysis.py --all
  ```
