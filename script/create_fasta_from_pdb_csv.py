"""
create_fasta_from_pdb_csv.py
----------------------------
Backward-compatibility wrapper importing from sequence_curator.py.
Primary module location: script/sequence_curator.py.
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from sequence_curator import (
    SequenceCurator,
    read_rcsb_csv,
    normalize_columns,
    abbreviate_organism,
    combine_asym_ids,
    filter_complete_complexes,
    build_complex_table,
    write_fasta,
    process_pdb_custom_report,
    main
)

__all__ = [
    "SequenceCurator",
    "read_rcsb_csv",
    "normalize_columns",
    "abbreviate_organism",
    "combine_asym_ids",
    "filter_complete_complexes",
    "build_complex_table",
    "write_fasta",
    "process_pdb_custom_report",
    "main"
]

if __name__ == "__main__":
    main()