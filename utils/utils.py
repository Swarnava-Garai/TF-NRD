import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Long JSON Parser with more details
def parse_features_long_TFDNA(json_file):
    rows = []

    with open(json_file, 'r') as f:
        data = json.load(f)

    results = data.get('results', [])
    logging.info(f"Total mapped entries found: {len(results)}")

    for item in results:
        entry = item.get('to', {})   # 🔥 key fix: unwrap "to"
        source_id = item.get('from') # original query ID

        if not entry:
            continue

        # -----------------------------
        # Basic Info
        # -----------------------------
        uid = entry.get('primaryAccession')

        protein_name = (
            entry.get('proteinDescription', {})
                 .get('recommendedName', {})
                 .get('fullName', {})
                 .get('value')
        )

        # -----------------------------
        # Organism Info
        # -----------------------------
        org_name = entry.get('organism', {})
        org_sci_name = org_name.get('scientificName')
        org_comm_name = org_name.get('commonName')

        # -----------------------------
        # Sequence
        # -----------------------------

        seq = entry.get('sequence', {})
        sequence = seq.get('value')
        seq_len = seq.get('length')

        # -----------------------------
        # Gene names
        # -----------------------------
        # genes = entry.get('genes', [])
        # if genes:
        #     gene_name = genes[0].get('geneName', {}).get('value')
        # else:
        #     gene_name = None

        # -----------------------------
        # Gene names (correct)
        # -----------------------------
        genes = entry.get('genes', [])

        gene_names = [
            g.get('geneName', {}).get('value')
            for g in genes
            if g.get('geneName')
        ]
        
        # -----------------------------
        # FUNCTION comments
        # -----------------------------
        comments = entry.get('comments', [])
        functions = []

        for c in comments:
            if c.get('commentType') == 'FUNCTION':
                for text in c.get('texts', []):
                    functions.append(text.get('value'))
        
        # -----------------------------
        # Subcellular location
        # -----------------------------
        locations = []

        for c in comments:
            if c.get('commentType') == 'SUBCELLULAR LOCATION':
                for loc in c.get('subcellularLocations', []):
                    location = loc.get('location', {}).get('value')
                    if location:
                        locations.append(location)

        # -----------------------------
        # Feature counts
        # -----------------------------

        for feat in entry.get('features', []):
            if feat.get('type') in ['Motif', 'Domain', 'Region']:

                loc = feat.get('location', {})
                start = loc.get('start', {}).get('value')
                end = loc.get('end', {}).get('value')

                rows.append({
                    'uniprot_id': uid,
                    'gene_name' : "; ".join(gene_names) if gene_names else None,
                    'protein_name': protein_name,
                    'org_sci_name': org_sci_name,
                    'org_comm_name': org_comm_name,
                    'sequence': sequence,
                    'sequence_length': seq_len,
                    'feature_type': feat.get('type'),
                    'feature_name': feat.get('description'),
                    'start': start,
                    'end': end,
                    'subcellular_location' : "; ".join(locations) if locations else None,
                    'function': " ".join(functions) if functions else None
                })

    return pd.DataFrame(rows)

# Long JSON Parser with more details
def parse_features_long_RBP(json_file):

    rows = []

    with open(json_file, 'r') as f:
        data = json.load(f)

    entries = data.get('results', data)
    logging.info(f"Total entries found: {len(entries)}")

    for entry in entries:
        # -----------------------------
        # Basic Info
        # -----------------------------
        uid = entry.get('primaryAccession')

        protein_name = (
            entry.get('proteinDescription', {})
                 .get('recommendedName', {})
                 .get('fullName', {})
                 .get('value')
        )

        # -----------------------------
        # Organism Name
        # -----------------------------
        org_name = entry.get('organism', {})
        org_sci_name = org_name.get('scientificName')
        org_comm_name = org_name.get('commonName')

        # -----------------------------
        # Sequence
        # -----------------------------

        seq = entry.get('sequence', {})
        sequence = seq.get('value')
        seq_len = seq.get('length')

        # -----------------------------
        # Gene names
        # -----------------------------
        # genes = entry.get('genes', [])
        # if genes:
        #     gene_name = genes[0].get('geneName', {}).get('value')
        # else:
        #     gene_name = None

        # -----------------------------
        # Gene names (correct)
        # -----------------------------
        genes = entry.get('genes', [])

        gene_names = [
            g.get('geneName', {}).get('value')
            for g in genes
            if g.get('geneName')
        ]
        
        # -----------------------------
        # FUNCTION comments
        # -----------------------------
        comments = entry.get('comments', [])
        functions = []

        for c in comments:
            if c.get('commentType') == 'FUNCTION':
                for text in c.get('texts', []):
                    functions.append(text.get('value'))
        
        # -----------------------------
        # Subcellular location
        # -----------------------------
        locations = []

        for c in comments:
            if c.get('commentType') == 'SUBCELLULAR LOCATION':
                for loc in c.get('subcellularLocations', []):
                    location = loc.get('location', {}).get('value')
                    if location:
                        locations.append(location)

        # -----------------------------
        # Feature counts
        # -----------------------------

        for feat in entry.get('features', []):
            if feat.get('type') in ['Motif', 'Domain', 'Region']:

                loc = feat.get('location', {})
                start = loc.get('start', {}).get('value')
                end = loc.get('end', {}).get('value')

                rows.append({
                    'uniprot_id': uid,
                    'gene_name' : "; ".join(gene_names) if gene_names else None,
                    'protein_name': protein_name,
                    'org_sci_name': org_sci_name,
                    'org_comm_name': org_comm_name,
                    'sequence': sequence,
                    'sequence_length': seq_len,
                    'feature_type': feat.get('type'),
                    'feature_name': feat.get('description'),
                    'start': start,
                    'end': end,
                    'subcellular_location' : "; ".join(locations) if locations else None,
                    'function': " ".join(functions) if functions else None
                })

    return pd.DataFrame(rows)


#Parser for PDB mapping
def extract_pdb_mapping(entry):
    mappings = []

    uid = entry.get('primaryAccession')

    for ref in entry.get('uniProtKBCrossReferences', []):
        if ref.get('database') == 'PDB':

            pdb_id = ref.get('id')

            chains_info = None
            for prop in ref.get('properties', []):
                if prop.get('key') == 'Chains':
                    chains_info = prop.get('value')

            if chains_info:
                for chain_entry in chains_info.split(','):

                    chain_entry = chain_entry.strip()

                    # Skip malformed entries
                    if '=' not in chain_entry:
                        continue

                    chain, region = chain_entry.split('=', 1)

                    # Handle missing or invalid regions
                    if '-' not in region:
                        continue

                    try:
                        start, end = region.split('-')
                        start = int(start)
                        end = int(end)
                    except:
                        continue

                    mappings.append({
                        'uniprot_id': uid,
                        'pdb_id': pdb_id,
                        'chain': chain.strip(),
                        'pdb_start': start,
                        'pdb_end': end
                    })

    return mappings

def build_pdb_df_TFDNA(input_data):
    import json
    import pandas as pd
    import logging

    rows = []

    # ============================
    # Load input
    # ============================
    if isinstance(input_data, str):
        with open(input_data) as f:
            data = json.load(f)
    else:
        data = input_data

    results = data.get('results', [])
    logging.info(f"Total mapped entries: {len(results)}")

    # ============================
    # Parse mappings
    # ============================
    for item in results:
        entries = item.get('to', [])
        source_id = item.get('from')

        # 🔥 Handle string case
        if isinstance(entries, str):
            logging.warning(f"Skipping string entry: {entries}")
            continue

        # 🔥 Handle dict case
        if isinstance(entries, dict):
            entries = [entries]

        for entry in entries:
            if not isinstance(entry, dict):
                continue

            mapped = extract_pdb_mapping(entry)

            for m in mapped:
                m['query_id'] = source_id

            rows.extend(mapped)

    return pd.DataFrame(rows)

def map_features_to_pdb(feature_df, pdb_df):
    merged = feature_df.merge(pdb_df, on='uniprot_id')

    # Overlap condition
    mapped = merged[
        (merged['start'] <= merged['pdb_end']) &
        (merged['end'] >= merged['pdb_start'])
    ].copy()

    return mapped

def build_structure_hierarchy(feature_df, pdb_df):
    results = []

    for _, pdb_row in pdb_df.iterrows():
        uid = pdb_row['uniprot_id']
        pdb_id = pdb_row['pdb_id']
        p_start = pdb_row['pdb_start']
        p_end = pdb_row['pdb_end']

        # Get features for this protein
        feats = feature_df[feature_df['uniprot_id'] == uid]

        # Filter features inside this PDB region
        feats_in_pdb = feats[
            (feats['start'] <= p_end) &
            (feats['end'] >= p_start)
        ].copy()

        # Separate by type
        domains = feats_in_pdb[feats_in_pdb['feature_type'] == 'Domain']
        motifs = feats_in_pdb[feats_in_pdb['feature_type'] == 'Motif']

        for _, dom in domains.iterrows():
            d_start, d_end = dom['start'], dom['end']

            # Find motifs inside this domain
            motifs_in_domain = motifs[
                (motifs['start'] >= d_start) &
                (motifs['end'] <= d_end)
            ]

            results.append({
                'uniprot_id': uid,
                'pdb_id': pdb_id,
                'pdb_region': f"{p_start}-{p_end}",
                'domain_name': dom['feature_name'],
                'domain_start': d_start,
                'domain_end': d_end,
                'motifs': list(motifs_in_domain['feature_name']),
                'motif_positions': list(zip(
                    motifs_in_domain['start'],
                    motifs_in_domain['end']
                ))
            })

    return pd.DataFrame(results)

def plot_architecture_advanced(protein_id, protein_df, feature_df, pdb_df, output_file="architecture.svg"):

    # =========================================================
    # 1. LOAD DATA
    # =========================================================
    prot = protein_df[protein_df['uniprot_id'] == protein_id].iloc[0]
    seq_len = prot['sequence_length']

    features = feature_df[feature_df['uniprot_id'] == protein_id]
    pdbs = pdb_df[pdb_df['uniprot_id'] == protein_id].copy()

    logging.info(f"Plotting architecture for {protein_id}")
    logging.info(f"Sequence length: {seq_len}")
    logging.info(f"Features: {len(features)}, PDBs: {len(pdbs)}")

    if pdbs.empty:
        logging.warning("No PDB structures found")
        return

    # =========================================================
    # 2. STACK PDB TRACKS (avoid overlap)
    # =========================================================
    pdbs = pdbs.sort_values('pdb_start')
    tracks = []

    pdbs['track'] = 0
    for i, row in pdbs.iterrows():
        placed = False
        for t, track in enumerate(tracks):
            if row['pdb_start'] > track[-1]:
                track.append(row['pdb_end'])
                pdbs.at[i, 'track'] = t
                placed = True
                break

        if not placed:
            tracks.append([row['pdb_end']])
            pdbs.at[i, 'track'] = len(tracks) - 1

    n_tracks = pdbs['track'].max() + 1

    # =========================================================
    # 3. DEFINE LAYOUT (dynamic spacing)
    # =========================================================
    y_seq = 5
    y_pdb_top = 4
    pdb_spacing = 0.35
    y_pdb_bottom = y_pdb_top - (n_tracks - 1) * pdb_spacing

    y_domain = y_pdb_bottom - 1.0
    y_motif  = y_domain - 1.0
    y_region = y_motif - 1.0

    # =========================================================
    # 4. INITIALIZE FIGURE
    # =========================================================
    fig, ax = plt.subplots(figsize=(15, 5))

    # =========================================================
    # 5. SEQUENCE TRACK
    # =========================================================
    ax.hlines(y_seq, 0, seq_len, linewidth=8, color='black', alpha=0.7)

    ax.text(
        seq_len / 2, y_seq + 0.4,
        f"UniProt AA Length: {seq_len} AA",
        fontsize=14,
        ha='center',
        fontweight='bold'
    )

    # =========================================================
    # 6. PDB TRACKS
    # =========================================================
    for _, row in pdbs.iterrows():
        y = y_pdb_top - row['track'] * pdb_spacing
        length = row['pdb_end'] - row['pdb_start']

        ax.hlines(
            y,
            row['pdb_start'],
            row['pdb_end'],
            linewidth=6,
            color='navy',
            alpha=0.8
        )

        # ---- Label (PDB + chain)
        if length > 5:
            mid = (row['pdb_start'] + row['pdb_end']) / 2
            label = f"{row['pdb_id']}:{row.get('chain','')}".strip(':')
            y_range = y + 0.15 if length > 30 else y - 0.25

            ax.text(
                mid,
                y + 0.15,
                label,
                ha='right',
                fontsize=9,
                rotation=0,
                fontweight='bold'
            )

        # ---- Range label (cleaner than start/end separately)
        if length > 5:
            # shifted above to avoid overlap with PDB label
            y_range = y + 0.15 if length > 150 else y - 0.25
            # Shifted right to avoid overlap with PDB label
            if length > 150:
                mid = (row['pdb_start'] + row['pdb_end']) / 2 + 10
            else:
                mid = (row['pdb_start'] + row['pdb_end']) / 2
            
            ax.text(
                mid,
                y_range,
                f"{row['pdb_start']}-{row['pdb_end']}",
                ha='left',
                fontsize=9
            )
# =========================================================
# 7. DOMAIN TRACK
# =========================================================
    domains = features[features['feature_type'] == 'Domain']

    for _, row in domains.iterrows():
        start, end = row['start'], row['end']
        length = end - start
        mid = (start + end) / 2

        # ---- Domain bar
        ax.hlines(
            y_domain,
            start,
            end,
            linewidth=6,
            color='#2ca02c'
        )

        # ---- Domain name (above)
        if length > 15:
            ax.text(
                mid,
                y_domain + 0.3,
                row['feature_name'],
                ha='center',
                fontsize=10,
                fontweight='bold'
            )

        # ---- Domain range (centered)
        if length > 30:
            ax.text(
            mid,
            y_domain - 0.25,
            f"{start}-{end}",
            ha='center',
            fontsize=9,
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')
        )

    # =========================================================
    # 8. MOTIF TRACK
    # =========================================================
    motifs = features[features['feature_type'] == 'Motif']

    last_x = -100

    for _, row in motifs.sort_values('start').iterrows():
        mid = (row['start'] + row['end']) / 2

        if mid - last_x > 20:  # spacing threshold
            ax.text(mid, y_motif - 0.25, f"{row['start']}-{row['end']}",
                    ha='center', fontsize=7)
            last_x = mid

            # ---- Motif bar
            ax.hlines(
                y_motif,
                start,
                end,
                linewidth=3,
                color='#d62728'
            )

            # ---- Motif name (only if reasonable size)
            if length >= 6:
                ax.text(
                    mid,
                    y_motif + 0.25,
                    row['feature_name'],
                    ha='center',
                    fontsize=8,
                    rotation=90
                )

            # ---- Motif range (key addition)
            # Show for all motifs, but small font
            ax.text(
                mid,
                y_motif - 0.25,
                f"{start}-{end}",
                fontsize=7,
                ha='center',
                bbox=dict(facecolor='white', alpha=0.5, edgecolor='none')
            )
    # =========================================================
    # 9. REGION TRACK
    # =========================================================
    regions = features[features['feature_type'] == 'Region']

    for _, row in regions.iterrows():
        ax.hlines(
            y_region,
            row['start'],
            row['end'],
            linewidth=3,
            color='#7f7f7f'
        )

    # =========================================================
    # 10. AXIS FORMATTING
    # =========================================================
    ax.set_yticks([
        y_seq,
        (y_pdb_top + y_pdb_bottom) / 2,
        y_domain,
        y_motif,
        y_region
    ])

    ax.set_yticklabels([
        f"UniProt ID:\n {protein_id}",
        "PDB structures",
        "Domains",
        "Motifs",
        "Regions"
    ], fontsize=12, fontweight='bold')

    ax.set_xlabel("Sequence Position", fontsize=14, fontweight='bold')
    ax.set_xlim(-5, seq_len + 5)
    ax.set_ylim(y_region - 0.5, y_seq + 0.6)

    # Clean axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()