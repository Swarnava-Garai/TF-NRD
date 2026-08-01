"""
sequence_curator.py
-------------------
Bioinformatics Sequence Curation & Filtering Suite for Protein and Nucleic Acid Datasets.
Provides the SequenceCurator class with static methods for reading RCSB/NCBI custom reports,
molecule type detection, sequence length quality filtering, chain identifier aggregation,
FASTA parsing/export, and complex table summarization.
"""

from pathlib import Path
import logging
import argparse
import sys
import pandas as pd
import numpy as np

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('sequence_curator.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("SequenceCurator")

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent.resolve()
INPUT_DIR = PROJECT_DIR / "input_data"
OUTPUT_DIR = PROJECT_DIR / "output_data"

# Standard column alias mappings for RCSB PDB custom reports
COLUMN_ALIASES = {
    'entry_id': ['Entry ID', 'pdb_id', 'PDB_ID', 'Structure ID', 'Identifier'],
    'entity_id': ['Entity ID', 'entity_id', 'Entity_ID'],
    'asym_id': ['Asym ID', 'asym_id'],
    'auth_asym_id': ['Auth Asym ID', 'auth_asym_id', 'Chain ID', 'Chain', 'chain_id'],
    'uniprot_id': ['Accession Code(s)', 'UniProt ID', 'uniprot_id', 'Accession Code', 'Database Accession'],
    'sequence': ['Sequence', 'Polymer Entity Sequence', 'sequence'],
    'sequence_length': ['Polymer Entity Sequence Length', 'Sequence Length', 'sequence_length', 'Length'],
    'macromolecule_type': ['Entity Macromolecule Type', 'Polymer Type', 'Entity Polymer Type', 'macromolecule_type', 'Type'],
    'molecular_weight': ['Molecular Weight (Entity)', 'Molecular Weight', 'molecular_weight'],
    'protein_name': ['Protein Name', 'Macromolecule Name', 'Entity Description', 'Structure Title', 'protein_name'],
    'resolution': ['Refinement Resolution (Å)', 'Resolution (Å)', 'Refinement Resolution', 'resolution'],
    'organism': ['Source Organism', 'Organism', 'source_organism'],
    'expression_host': ['Expression Host', 'expression_host']
}


class SequenceCurator:
    """
    Utility suite for sequence curation, quality checking, chain identifier merging,
    macromolecule classification, and report generation. Implemented as static methods.
    """

    @staticmethod
    def detect_molecule_type(sequence: str) -> str:
        """
        Detects molecule type from raw sequence characters ('DNA', 'RNA', or 'PROTEIN').
        """
        seq = str(sequence).upper().replace(" ", "").replace("\n", "").replace("\r", "")
        if not seq or seq in ['NAN', 'NONE']:
            return "UNKNOWN"

        dna_set = set("ATGC")
        rna_set = set("AUGC")
        seq_chars = set(seq)

        if seq_chars.issubset(dna_set):
            return "DNA"
        elif seq_chars.issubset(rna_set):
            return "RNA"
        else:
            return "PROTEIN"

    @staticmethod
    def extract_chain_id(header: str) -> str:
        """
        Extracts chain identifier from a FASTA header string.
        """
        parts = str(header).split("|")
        for part in parts:
            if "CHAIN" in part.upper():
                return part.strip()
        
        # Fallback: check if header contains : separator e.g. 1A02:FJ
        if ":" in header.split()[0]:
            return header.split()[0].split(":")[1]
        return "Unknown"

    @staticmethod
    def parse_fasta(file_path: Path) -> list:
        """
        Parses a FASTA file into a list of (header, sequence) tuples.
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")

        records = []
        with open(file_path, "r", encoding="utf-8", errors="replace") as f:
            header = None
            seq_lines = []

            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        records.append((header, "".join(seq_lines)))
                    header = line[1:]
                    seq_lines = []
                else:
                    seq_lines.append(line)

            if header:
                records.append((header, "".join(seq_lines)))

        return records

    @staticmethod
    def filter_fasta_by_length(
        input_fasta: Path,
        output_valid_fasta: Path = None,
        output_invalid_fasta: Path = None,
        min_protein_len: int = 30,
        min_na_len: int = 5
    ) -> dict:
        """
        Quality checks a FASTA file and filters out entries that do not satisfy length criteria.
        Protein <= min_protein_len or Nucleic Acid <= min_na_len are flagged.
        """
        records = SequenceCurator.parse_fasta(input_fasta)

        valid_records = []
        invalid_records = []

        for header, seq in records:
            seq_clean = seq.replace(" ", "").replace("\n", "").replace("\r", "")
            length = len(seq_clean)
            mol_type = SequenceCurator.detect_molecule_type(seq_clean)

            reason = None
            if mol_type == "PROTEIN" and length <= min_protein_len:
                reason = f"Protein length <= {min_protein_len}"
            elif mol_type in ["DNA", "RNA"] and length <= min_na_len:
                reason = f"Nucleic acid length <= {min_na_len}"

            if reason:
                invalid_records.append((header, seq_clean, mol_type, length, reason))
            else:
                valid_records.append((header, seq_clean, mol_type, length))

        if output_valid_fasta:
            out_v = Path(output_valid_fasta)
            out_v.parent.mkdir(parents=True, exist_ok=True)
            with open(out_v, "w", encoding="utf-8") as f:
                for header, seq, _, _ in valid_records:
                    f.write(f">{header}\n{seq}\n")

        if output_invalid_fasta:
            out_inv = Path(output_invalid_fasta)
            out_inv.parent.mkdir(parents=True, exist_ok=True)
            with open(out_inv, "w", encoding="utf-8") as f:
                for header, seq, m_type, l_val, r_msg in invalid_records:
                    f.write(f">{header} | {m_type} | len={l_val} | {r_msg}\n{seq}\n")

        logger.info(f"FASTA Quality Control for {Path(input_fasta).name}: {len(valid_records)} valid, {len(invalid_records)} invalid/flagged entries")
        return {
            'valid_records': valid_records,
            'invalid_records': invalid_records,
            'valid_count': len(valid_records),
            'invalid_count': len(invalid_records)
        }

    @staticmethod
    def read_rcsb_csv(csv_path: Path) -> pd.DataFrame:
        """
        Intelligently reads RCSB PDB custom report CSV files (1-header or 2-header format).
        """
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(f"Input CSV file not found: {csv_path}")

        with open(csv_path, 'r', encoding='utf-8', errors='replace') as f:
            l1 = f.readline().strip()
            l2 = f.readline().strip()

        if any(k in l2 for k in ['Entry ID', 'Auth Asym ID', 'Entity ID', 'Accession Code', 'Refinement Resolution']):
            logger.info(f"Reading RCSB PDB CSV with 2-header rows: {csv_path.name}")
            df = pd.read_csv(csv_path, skiprows=1, header=0, low_memory=False)
        else:
            logger.info(f"Reading RCSB PDB CSV with standard 1-header row: {csv_path.name}")
            df = pd.read_csv(csv_path, header=0, low_memory=False)

        return df

    @staticmethod
    def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
        """
        Maps varying RCSB PDB column headers to standardized internal column names.
        """
        renames = {}
        existing_cols = list(df.columns)

        for std_name, aliases in COLUMN_ALIASES.items():
            for col in existing_cols:
                if col in aliases or col.strip() in aliases:
                    if col != std_name:
                        renames[col] = std_name
                    break

        return df.rename(columns=renames)

    @staticmethod
    def abbreviate_organism(val: str) -> str:
        """
        Abbreviates Source Organism strings (e.g. Escherichia coli -> E. coli, Homo sapiens -> H. sapiens).
        """
        if pd.isna(val) or not str(val).strip() or str(val).strip() in ['nan', 'None', 'Unknown']:
            return 'Unknown'

        parts = str(val).strip().split()
        if len(parts) == 1:
            return parts[0]
        return f"{parts[0][0].upper()}. {' '.join(parts[1:])}"

    @staticmethod
    def preprocess_df(df: pd.DataFrame) -> pd.DataFrame:
        """
        Cleans whitespace and forward-fills entry-level metadata across multi-entity rows.
        """
        df_clean = df.copy()

        for col in df_clean.select_dtypes(include=['object', 'string']).columns:
            df_clean[col] = df_clean[col].astype(str).str.strip()
            df_clean[col] = df_clean[col].replace({'nan': np.nan, 'None': np.nan, '': np.nan})

        ffill_cols = ['entry_id', 'uniprot_id', 'protein_name', 'resolution', 'organism', 'expression_host']
        for col in ffill_cols:
            if col in df_clean.columns:
                df_clean[col] = df_clean[col].ffill()

        return df_clean

    @staticmethod
    def merge_metadata_report(seq_df: pd.DataFrame, metadata_csv_path: Path) -> pd.DataFrame:
        """
        Merges additional structure metadata (e.g. from TFNRDv1.0_custom_report_pdb.csv) on entry_id.
        """
        meta_df_raw = SequenceCurator.read_rcsb_csv(metadata_csv_path)
        meta_df = SequenceCurator.normalize_columns(meta_df_raw)

        if 'entry_id' in meta_df.columns:
            meta_df['entry_id'] = meta_df['entry_id'].ffill()

        merge_cols = ['entry_id']
        meta_feature_cols = [c for c in ['resolution', 'organism', 'expression_host', 'protein_name'] if c in meta_df.columns and c not in seq_df.columns]

        if meta_feature_cols:
            meta_sub = meta_df[merge_cols + meta_feature_cols].drop_duplicates(subset=['entry_id'])
            seq_df = seq_df.merge(meta_sub, on='entry_id', how='left')
            logger.info(f"Merged metadata columns {meta_feature_cols} from {Path(metadata_csv_path).name}")

        return seq_df

    @staticmethod
    def filter_complete_complexes(df: pd.DataFrame, min_protein_len: int = 30, min_na_len: int = 5) -> pd.DataFrame:
        """
        Filters Entry IDs to keep complexes containing BOTH a qualifying protein (length >= min_protein_len)
        AND a qualifying nucleic acid (DNA or RNA, length >= min_na_len).
        """
        if 'macromolecule_type' not in df.columns or 'sequence_length' not in df.columns:
            return df

        df_copy = df.copy()
        df_copy['seq_len_num'] = pd.to_numeric(df_copy['sequence_length'], errors='coerce').fillna(0)

        c_poly = (df_copy['macromolecule_type'] == 'polypeptide(L)') & (df_copy['seq_len_num'] >= min_protein_len)
        c_rna = (df_copy['macromolecule_type'] == 'polyribonucleotide') & (df_copy['seq_len_num'] >= min_na_len)
        c_dna = (df_copy['macromolecule_type'] == 'polydeoxyribonucleotide') & (df_copy['seq_len_num'] >= min_na_len)

        df_filtered = df_copy[c_poly | c_rna | c_dna].copy()

        poly_ids = set(df_filtered[df_filtered['macromolecule_type'] == 'polypeptide(L)']['entry_id'].unique())
        na_ids = set(df_filtered[df_filtered['macromolecule_type'].isin(['polyribonucleotide', 'polydeoxyribonucleotide'])]['entry_id'].unique())
        valid_ids = poly_ids & na_ids

        if not valid_ids:
            logger.warning("No Entry IDs contain both qualifying protein and nucleic acid macromolecule types.")
            return df.iloc[0:0].copy()

        result = df[df['entry_id'].isin(valid_ids)].copy().reset_index(drop=True)
        logger.info(f"Filtered {len(valid_ids)} complete complex Entry IDs containing both protein >= {min_protein_len} and NA >= {min_na_len}")
        return result

    @staticmethod
    def combine_asym_ids(df: pd.DataFrame, condition_dict: dict = None, chain_sep: str = '') -> pd.DataFrame:
        """
        Combines chain IDs (auth_asym_id / asym_id) for identical entity sequences within each PDB entry.
        """
        df_target = df.copy()

        if condition_dict:
            mask = pd.Series([True] * len(df_target), index=df_target.index)
            for col, val in condition_dict.items():
                col_name = col
                if col not in df_target.columns:
                    for std_name, aliases in COLUMN_ALIASES.items():
                        if col in aliases:
                            col_name = std_name
                            break

                if col_name in df_target.columns:
                    mask &= (df_target[col_name].astype(str).str.lower() == str(val).lower())

            df_matching = df_target[mask].copy()
            df_non_matching = df_target[~mask].copy()
        else:
            df_matching = df_target
            df_non_matching = pd.DataFrame()

        if df_matching.empty:
            return df_target

        group_cols = []
        for c in ['entry_id', 'sequence', 'macromolecule_type', 'sequence_length']:
            if c in df_matching.columns:
                group_cols.append(c)

        if not group_cols:
            group_cols = [c for c in ['entry_id', 'entity_id'] if c in df_matching.columns]

        chain_col = 'auth_asym_id' if 'auth_asym_id' in df_matching.columns else ('asym_id' if 'asym_id' in df_matching.columns else None)

        def _merge_chains(series):
            valid = sorted(set(str(x).strip() for x in series.dropna() if str(x).strip() not in ['nan', 'None', '']))
            return chain_sep.join(valid)

        agg_dict = {}
        if chain_col:
            agg_dict[chain_col] = _merge_chains

        for col in df_matching.columns:
            if col not in group_cols and col not in agg_dict:
                agg_dict[col] = 'first'

        df_filtered = df_matching.groupby(group_cols, as_index=False, dropna=False).agg(agg_dict)

        if chain_col and chain_col in df_filtered.columns:
            df_filtered[chain_col] = df_filtered[chain_col].astype(str).str.replace(' ', '')

        result = pd.concat([df_non_matching, df_filtered], ignore_index=True)
        if 'entry_id' in result.columns:
            result = result.sort_values(['entry_id']).reset_index(drop=True)

        return result

    @staticmethod
    def filter_macromolecule_type(df: pd.DataFrame, mol_type: str = 'polypeptide(L)') -> pd.DataFrame:
        """
        Filters DataFrame for specified macromolecule type.
        """
        if 'macromolecule_type' not in df.columns:
            return df

        if not mol_type or mol_type.lower() == 'all':
            return df

        mask = df['macromolecule_type'].astype(str).str.lower() == mol_type.lower()
        return df[mask].copy()

    @staticmethod
    def build_complex_table(df: pd.DataFrame) -> pd.DataFrame:
        """
        Constructs a one-row-per-Entry-ID summary table of TF-DNA / TF-RNA complexes.
        Columns: Entry ID, Chains, protein_name, protein_length, DNA_name, DNA_length, Refinement Resolution (Å), Source Organism
        """
        df_copy = df.copy()
        df_copy['seq_len_str'] = pd.to_numeric(df_copy['sequence_length'], errors='coerce').astype('Int64').astype(str)

        res_col = 'resolution' if 'resolution' in df_copy.columns else df_copy.columns[0]
        org_col = 'organism' if 'organism' in df_copy.columns else df_copy.columns[0]
        name_col = 'protein_name' if 'protein_name' in df_copy.columns else 'macromolecule_type'
        chain_col = 'auth_asym_id' if 'auth_asym_id' in df_copy.columns else ('asym_id' if 'asym_id' in df_copy.columns else df_copy.columns[0])

        base_info = df_copy.groupby('entry_id').agg({
            res_col: 'first',
            org_col: 'first'
        })

        def aggregate_type(mac_type, prefix):
            subset = df_copy[df_copy['macromolecule_type'] == mac_type]
            if subset.empty:
                return pd.DataFrame()
            return subset.groupby('entry_id').agg({
                name_col: lambda x: '+'.join(sorted(set(str(v) for v in x.dropna() if str(v).strip() not in ['nan', 'None', '']))),
                'seq_len_str': lambda x: '+'.join(sorted(set(str(v) for v in x.dropna() if str(v).strip() not in ['nan', 'None', '']))),
                chain_col: lambda x: ''.join(sorted(set(str(v) for v in x.dropna() if str(v).strip() not in ['nan', 'None', ''])))
            }).rename(columns={
                name_col: f'{prefix}_name',
                'seq_len_str': f'{prefix}_length',
                chain_col: f'{prefix}_chain'
            })

        protein = aggregate_type('polypeptide(L)', 'protein')
        dna = aggregate_type('polydeoxyribonucleotide', 'DNA')
        rna = aggregate_type('polyribonucleotide', 'RNA')

        output_df = base_info.join([protein, dna, rna]).reset_index()

        chain_cols = [c for c in ['protein_chain', 'DNA_chain', 'RNA_chain'] if c in output_df.columns]
        output_df['Chains'] = output_df[chain_cols].fillna('').agg(':'.join, axis=1).str.strip(':')

        if org_col in output_df.columns:
            output_df[org_col] = output_df[org_col].apply(SequenceCurator.abbreviate_organism)

        final_rename = {
            'entry_id': 'Entry ID',
            res_col: 'Refinement Resolution (Å)',
            org_col: 'Source Organism'
        }
        output_df = output_df.rename(columns=final_rename)
        return output_df

    @staticmethod
    def export_complex_table(df: pd.DataFrame, output_file: Path) -> pd.DataFrame:
        """
        Exports the 1-row-per-Entry-ID complex summary table to CSV.
        """
        tbl = SequenceCurator.build_complex_table(df)
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        tbl.to_csv(output_file, index=False)
        logger.info(f"Successfully exported complex table ({len(tbl)} entries) to {output_file}")
        return tbl

    @staticmethod
    def export_fasta(df: pd.DataFrame, output_file: Path, mol_type: str = 'polypeptide(L)', header_format: str = None, line_width: int = 0) -> int:
        """
        Exports filtered sequence dataset to FASTA format.
        """
        df_sub = SequenceCurator.filter_macromolecule_type(df, mol_type)

        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        count = 0
        with open(output_file, 'w', encoding='utf-8') as f:
            for idx, row in df_sub.iterrows():
                seq = str(row.get('sequence', '')).strip().replace('\n', '').replace('\r', '')
                if not seq or seq.lower() in ['nan', 'none']:
                    continue

                entry_id = str(row.get('entry_id', 'NA'))
                chain_id = str(row.get('auth_asym_id', row.get('asym_id', 'NA')))
                length = str(row.get('sequence_length', len(seq)))
                if length.endswith('.0'):
                    length = length[:-2]

                uniprot_id = str(row.get('uniprot_id', 'NA'))

                if header_format:
                    try:
                        header = header_format.format(
                            entry_id=entry_id,
                            chain_id=chain_id,
                            length=length,
                            uniprot_id=uniprot_id,
                            **row.to_dict()
                        )
                        if not header.startswith('>'):
                            header = '>' + header
                    except KeyError:
                        header = f">{entry_id}:{chain_id} | {length}"
                else:
                    header = f">{entry_id}:{chain_id} | {length}"

                f.write(f"{header}\n")

                if line_width > 0:
                    for i in range(0, len(seq), line_width):
                        f.write(f"{seq[i:i+line_width]}\n")
                else:
                    f.write(f"{seq}\n")

                count += 1

        logger.info(f"Successfully wrote {count} FASTA sequence entries to {output_file}")
        return count

    @staticmethod
    def export_csv(df: pd.DataFrame, output_file: Path, mol_type: str = 'polypeptide(L)') -> pd.DataFrame:
        """
        Exports cleaned sequence dataset to CSV file.
        """
        df_sub = SequenceCurator.filter_macromolecule_type(df, mol_type)

        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        export_rename = {
            'entry_id': 'Entry ID',
            'entity_id': 'Entity ID',
            'auth_asym_id': 'Auth Asym ID',
            'asym_id': 'Asym ID',
            'uniprot_id': 'Accession Code(s)',
            'sequence': 'Sequence',
            'sequence_length': 'Polymer Entity Sequence Length',
            'macromolecule_type': 'Entity Macromolecule Type',
            'resolution': 'Refinement Resolution (Å)',
            'organism': 'Source Organism',
            'expression_host': 'Expression Host'
        }

        df_out = df_sub.rename(columns={k: v for k, v in export_rename.items() if k in df_sub.columns})
        df_out.to_csv(output_file, index=False)
        logger.info(f"Successfully exported {len(df_out)} rows to CSV: {output_file}")
        return df_out


# Package Top-Level API Wrappers
def process_pdb_custom_report(
    input_csv_path: Path,
    metadata_csv_path: Path = None,
    output_fasta_path: Path = None,
    output_csv_path: Path = None,
    complex_table_path: Path = None,
    filter_complexes: bool = False,
    min_protein_len: int = 30,
    min_na_len: int = 5,
    mol_type: str = 'polypeptide(L)',
    chain_sep: str = '',
    line_width: int = 0
) -> pd.DataFrame:
    """
    Master pipeline function to read RCSB PDB custom report CSV, forward-fill metadata,
    merge structure metadata, filter complete complexes, combine chain IDs,
    and write output FASTA, CSV, and complex table datasets.
    """
    df_raw = SequenceCurator.read_rcsb_csv(input_csv_path)
    df = SequenceCurator.normalize_columns(df_raw)
    df = SequenceCurator.preprocess_df(df)

    if metadata_csv_path and Path(metadata_csv_path).exists():
        df = SequenceCurator.merge_metadata_report(df, metadata_csv_path)

    if filter_complexes:
        df = SequenceCurator.filter_complete_complexes(df, min_protein_len=min_protein_len, min_na_len=min_na_len)

    conditions = [
        {'macromolecule_type': 'polypeptide(L)'},
        {'macromolecule_type': 'polydeoxyribonucleotide'},
        {'macromolecule_type': 'polyribonucleotide'}
    ]
    for cond in conditions:
        df = SequenceCurator.combine_asym_ids(df, cond, chain_sep=chain_sep)

    if output_fasta_path:
        SequenceCurator.export_fasta(df, output_fasta_path, mol_type=mol_type, line_width=line_width)

    if output_csv_path:
        SequenceCurator.export_csv(df, output_csv_path, mol_type=mol_type)

    if complex_table_path:
        SequenceCurator.export_complex_table(df, complex_table_path)

    return df


def main():
    parser = argparse.ArgumentParser(description="SequenceCurator: Curation, Quality Checking, FASTA & Complex Table Generator.")
    parser.add_argument("-i", "--input", default=None, help="Input RCSB PDB Custom Sequence CSV file path")
    parser.add_argument("-b", "--metadata-csv", default=None, help="Optional structure metadata CSV file path to merge")
    parser.add_argument("-o", "--output-fasta", default=None, help="Output FASTA file path")
    parser.add_argument("-c", "--output-csv", default=None, help="Output CSV file path")
    parser.add_argument("-t", "--complex-table-csv", default=None, help="Output Complex Table CSV file path")
    parser.add_argument("-f", "--filter-complexes", action="store_true", help="Filter complete protein-NA complexes")
    parser.add_argument("--check-fasta", default=None, help="FASTA file path to quality check by sequence length")
    parser.add_argument("--min-protein-len", type=int, default=30, help="Minimum protein sequence length cutoff (default: 30)")
    parser.add_argument("--min-na-len", type=int, default=5, help="Minimum nucleic acid sequence length cutoff (default: 5)")
    parser.add_argument("-m", "--mol-type", default="polypeptide(L)", help="Macromolecule type filter (default: polypeptide(L))")
    parser.add_argument("-s", "--chain-sep", default="", help="Separator for combined chain IDs (default: '')")
    parser.add_argument("-w", "--line-width", type=int, default=0, help="FASTA sequence line wrap width (0 for single line)")

    args = parser.parse_args()

    # Mode 1: Quality check existing FASTA file
    if args.check_fasta:
        logger.info(f"Quality checking FASTA file: {args.check_fasta}")
        res = SequenceCurator.filter_fasta_by_length(
            args.check_fasta,
            min_protein_len=args.min_protein_len,
            min_na_len=args.min_na_len
        )
        print(f"\n=== FASTA Quality Check Summary ===")
        print(f"Valid Sequences:   {res['valid_count']}")
        print(f"Flagged Sequences: {res['invalid_count']}")
        if res['invalid_records']:
            print("\nFlagged Entries:")
            for h, _, m_type, l_val, r_msg in res['invalid_records'][:10]:
                print(f"  {h[:40]:<40} | {m_type:<8} | len={l_val:<5} | {r_msg}")
        return

    # Mode 2: Process RCSB Custom Report CSV
    input_file = Path(args.input) if args.input else INPUT_DIR / "TFNRDv1.0_custom_report_seq.csv"
    metadata_file = Path(args.metadata_csv) if args.metadata_csv else SCRIPT_DIR.parent.parent / "swarnava_TF_work" / "shrikant_script" / "TFNRDv1.0_custom_report_pdb.csv"

    output_fasta = Path(args.output_fasta) if args.output_fasta else OUTPUT_DIR / "TFNRDv1.0_Sequence_367_oneline.fasta"
    output_csv = Path(args.output_csv) if args.output_csv else OUTPUT_DIR / "TFNRDV1.0_Sequence_367.csv"
    complex_table_csv = Path(args.complex_table_csv) if args.complex_table_csv else OUTPUT_DIR / "TFNRDv1.0_complex_table.csv"

    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)

    logger.info("Executing SequenceCurator Pipeline...")
    process_pdb_custom_report(
        input_csv_path=input_file,
        metadata_csv_path=metadata_file if metadata_file.exists() else None,
        output_fasta_path=output_fasta,
        output_csv_path=output_csv,
        complex_table_path=complex_table_csv if args.complex_table_csv or args.filter_complexes else None,
        filter_complexes=args.filter_complexes,
        min_protein_len=args.min_protein_len,
        min_na_len=args.min_na_len,
        mol_type=args.mol_type,
        chain_sep=args.chain_sep,
        line_width=args.line_width
    )
    logger.info("SequenceCurator Pipeline Completed Successfully!")


if __name__ == "__main__":
    main()
