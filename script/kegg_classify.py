#!/usr/bin/env python3
"""
Classify KEGG disease annotations into major disease classes.

- Supports KEGG DISEASE entries (e.g., H00224) via CATEGORY parsing.
- Supports KEGG Human Disease pathway IDs (e.g., hsa05010) via BRITE hierarchy br:br08902.
- Input: CSV with a column of KEGG IDs (default: 'kegg_id').
- Output: CSV with kegg_id, disease_name, major_class, subclass, source.

Usage:
  python kegg_classify.py --in input.csv --out output.csv --col kegg_id
"""

import argparse
import csv
import re
import sys
import time
from collections import defaultdict
from typing import Dict, Tuple, Optional

import pandas as pd
import requests

KEGG_BASE = "https://rest.kegg.jp"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "KEGG-Classifier/1.0"})

# --------- Helpers ---------

def http_get(url: str, max_retries: int = 4, backoff: float = 0.8) -> str:
    for i in range(max_retries):
        try:
            r = SESSION.get(url, timeout=30)
            if r.status_code == 200 and r.text:
                return r.text
            # Retry on non-200 or empty body
        except requests.RequestException:
            pass
        time.sleep(backoff * (2 ** i))
    raise RuntimeError(f"Failed to GET after {max_retries} retries: {url}")

def is_pathway_id(kegg_id: str) -> bool:
    # hsaXXXXX (human pathways)
    return bool(re.fullmatch(r"hsa\d{5}", kegg_id))

def is_disease_id(kegg_id: str) -> bool:
    # HXXXXX (KEGG DISEASE entries)
    return bool(re.fullmatch(r"H\d{5}", kegg_id))

# --------- BRITE br:br08902 parsing (Human Diseases pathways) ---------

def fetch_brite_human_diseases() -> str:
    # KEGG BRITE hierarchy for Human Diseases pathways
    # Contains top-level (A), subclass (B), and pathway entries with hsa codes under (C)
    url = f"{KEGG_BASE}/get/br:br08902"
    return http_get(url)

def parse_brite_hierarchy(br_text: str) -> Dict[str, Tuple[str, Optional[str]]]:
    """
    Parse br:br08902 text.
    Returns mapping: pathway_id -> (Major_Class_A, Subclass_B or None)
    """
    mapping: Dict[str, Tuple[str, Optional[str]]] = {}
    major = None
    sub = None
    for raw in br_text.splitlines():
        line = raw.rstrip()
        if not line:
            continue
        # Lines start with 'A', 'B', 'C' followed by two spaces and content in KEGG 'flat' text
        if line.startswith("A "):
            major = line[2:].strip()
            sub = None
        elif line.startswith("B "):
            sub = line[2:].strip()
        elif line.startswith("C "):
            # C entries hold one or more pathway IDs and names, e.g.:
            # C    hsa05010  Alzheimer disease
            parts = line[2:].strip().split()
            if len(parts) >= 2:
                pid = parts[0]
                if re.fullmatch(r"hsa\d{5}", pid):
                    mapping[pid] = (major or "Unknown", sub)
    return mapping

# --------- Name resolvers ---------

_name_cache: Dict[str, str] = {}
_class_cache: Dict[str, Tuple[str, Optional[str], str]] = {}  # id -> (major, subclass, source)

def get_pathway_name(hsa_id: str) -> Optional[str]:
    if hsa_id in _name_cache:
        return _name_cache[hsa_id]
    # KEGG list returns: <entry>\t<name>
    # Example: hsa05010\tAlzheimer disease
    url = f"{KEGG_BASE}/list/{hsa_id}"
    try:
        txt = http_get(url)
    except RuntimeError:
        return None
    line = txt.strip().splitlines()
    if not line:
        return None
    fields = line[0].split("\t")
    if len(fields) >= 2:
        name = fields[1].strip()
        _name_cache[hsa_id] = name
        return name
    return None

def parse_kegg_flat(text: str) -> Dict[str, list]:
    """
    Parse a KEGG flat file (simple field:value lines).
    Returns dict of field -> list of values (for multi-line fields).
    """
    data = defaultdict(list)
    current_key = None
    for raw in text.splitlines():
        if not raw.strip():
            continue
        if len(raw) >= 12 and raw[:12].strip().isalpha() and raw[12] == ' ':
            # New field starts (KEY at cols 1-12)
            current_key = raw[:12].strip()
            value = raw[12:].strip()
            if value:
                data[current_key].append(value)
        else:
            # Continuation line
            if current_key is not None:
                data[current_key].append(raw.strip())
    return data

def get_disease_name_and_category(h_id: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    For KEGG DISEASE entry Hxxxxx:
      - NAME field for disease name (first item)
      - CATEGORY field like "Neurodegenerative disease; Alzheimer disease"
        -> major class before semicolon, subclass after semicolon (first pair)
    """
    url = f"{KEGG_BASE}/get/{h_id}"
    try:
        txt = http_get(url)
    except RuntimeError:
        return None, None, None
    data = parse_kegg_flat(txt)
    # NAME
    name = None
    if "NAME" in data and data["NAME"]:
        # First NAME line; strip trailing semicolons
        name = re.sub(r"\s*;$", "", data["NAME"][0]).strip()
    # CATEGORY
    major = None
    subclass = None
    if "CATEGORY" in data and data["CATEGORY"]:
        # CATEGORY lines can repeat; take the first, split on ';'
        first_cat = data["CATEGORY"][0]
        parts = [p.strip() for p in first_cat.split(";")]
        if parts:
            major = parts[0] or None
        if len(parts) > 1:
            subclass = parts[1] or None
    return name, major, subclass

# --------- Main mapping function ---------

def classify_ids(ids: pd.Series, br_map: Dict[str, Tuple[str, Optional[str]]]) -> pd.DataFrame:
    rows = []
    for kid in ids.astype(str).str.strip():
        if not kid:
            rows.append((kid, None, None, None, None))
            continue

        disease_name = None
        major = None
        subclass = None
        source = None

        try:
            if is_pathway_id(kid):
                # Pathway-based classification via BRITE
                source = "pathway"
                disease_name = get_pathway_name(kid)
                if kid in br_map:
                    major, subclass = br_map[kid][0], br_map[kid][1]
                else:
                    # Fallback: unknown in br08902
                    major, subclass = None, None

            elif is_disease_id(kid):
                # KEGG DISEASE entry classification via CATEGORY
                source = "disease"
                disease_name, major, subclass = get_disease_name_and_category(kid)

            else:
                # Try best-effort generic lookups:
                # (a) Try as pathway
                if is_pathway_id(kid.lower()):
                    kid = kid.lower()
                    source = "pathway"
                    disease_name = get_pathway_name(kid)
                    if kid in br_map:
                        major, subclass = br_map[kid][0], br_map[kid][1]
                # (b) Try as disease
                elif is_disease_id(kid.upper()):
                    kid = kid.upper()
                    source = "disease"
                    disease_name, major, subclass = get_disease_name_and_category(kid)
                else:
                    # Final attempt: get name only (no class)
                    # works for many KEGG entries
                    url = f"{KEGG_BASE}/list/{kid}"
                    try:
                        txt = http_get(url)
                        parts = txt.strip().split("\t")
                        if len(parts) >= 2:
                            disease_name = parts[1].strip()
                    except Exception:
                        pass
                    source = "unknown"

        except Exception:
            # Keep row with Nones; do not crash entire run
            pass

        rows.append((kid, disease_name, major, subclass, source))

        # Gentle pacing to be nice to KEGG server
        time.sleep(0.15)

    return pd.DataFrame(rows, columns=["kegg_id", "disease_name", "major_class", "subclass", "source"])

# --------- CLI ---------

def main():
    ap = argparse.ArgumentParser(description="Categorize KEGG disease annotations into major classes.")
    ap.add_argument("--in", dest="inp", required=True, help="Input CSV path")
    ap.add_argument("--out", dest="out", required=True, help="Output CSV path")
    ap.add_argument("--col", dest="col", default="kegg_id", help="Column name containing KEGG IDs (default: kegg_id)")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    if args.col not in df.columns:
        print(f"ERROR: Column '{args.col}' not found in {args.inp}. Columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Build BRITE map for Human Diseases pathways
    try:
        br_txt = fetch_brite_human_diseases()
        br_map = parse_brite_hierarchy(br_txt)
    except Exception as e:
        print(f"WARNING: Failed to fetch/parse br:br08902; pathway major classes may be missing. ({e})", file=sys.stderr)
        br_map = {}

    out_df = classify_ids(df[args.col], br_map)
    # Preserve input order by merging back if needed
    merged = df.copy()
    merged = merged.merge(out_df, how="left", left_on=args.col, right_on="kegg_id")
    # Drop duplicate key column if present
    if "kegg_id_y" in merged.columns or "kegg_id_x" in merged.columns:
        merged = merged.drop(columns=[c for c in ["kegg_id_y"] if c in merged.columns])
        merged = merged.rename(columns={"kegg_id_x": "kegg_id"}) if "kegg_id_x" in merged.columns else merged

    merged.to_csv(args.out, index=False)
    print(f"Saved: {args.out}")

if __name__ == "__main__":
    main()
