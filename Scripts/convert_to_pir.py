import pandas as pd
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio import AlignIO

def merge_df_for_pir(df_light, df_heavy):
    """# Combine heavy and light chain dataframes into one for .pir creation."""
    merged_df = pd.merge(df_light, df_heavy, on=["pdb", "template", "H_isotype_clean", "HC_species"])
    merged_df

    # Filter redundant columns (making a copy of df)
    df_combined = merged_df[[
        "pdb", "template", "Lchain", "Hchain", "L_chain_first_residue", 
        "H_chain_first_residue", "L_chain_last_residue", "H_chain_last_residue"
    ]].copy()

    # Add empty columns to new copy
    df_combined["light_chain_is_first"] = None
    df_combined["start_point"] = None
    df_combined["start_letter"] = None
    df_combined["gapped_seq_1"] = None
    df_combined["end_point"] = None
    df_combined["end_letter"] = None
    df_combined["gapped_seq_2"] = None

    return df_combined


def cif_parse(pdb: str, cif_fname: str, filepath: str):
    """Parses .cif and determines correct chain order for MODELLER."""
    from Bio.PDB.MMCIFParser import MMCIFParser
    import os

    # Set full path
    full_cif_path = os.path.join(filepath, cif_fname)

    # Initialise parser
    parser = MMCIFParser(QUIET=True)

    # Load structure
    structure = parser.get_structure(pdb, full_cif_path)

    # Extract first model (usually there's only one)
    model = next(structure.get_models())

    # Get chains in order
    cif_chain_order = [chain.id for chain in model]

    print(cif_chain_order)

    return cif_chain_order


def relevant_chains(df, v_cif_chain_order, c_cif_chain_order, v_template, c_template):
    """Extracts relevant chain info from dataframe"""

    # Using iterrows for stable behaviour in pandas
    for idx, row in df.iterrows():
        if row["pdb"] == v_template:
            l_chain_letter = row["Lchain"]
            h_chain_letter = row["Hchain"]

            v_relevant_chains = [chain_id for chain_id in v_cif_chain_order if chain_id in (h_chain_letter, l_chain_letter)]
            print(v_relevant_chains)

            if v_relevant_chains[0] == l_chain_letter:
                # assign True/False to light_chain_is_first column
                df.at[idx, "light_chain_is_first"] = True
                df.at[idx, "start_point"] = int(row["L_chain_first_residue"])
                df.at[idx, "start_letter"] = row["Lchain"]
                df.at[idx, "end_point"] = int(row["H_chain_last_residue"])
                df.at[idx, "end_letter"] = row["Hchain"]
            else:
                df.at[idx, "light_chain_is_first"] = False
                df.at[idx, "start_point"] = int(row["H_chain_first_residue"])
                df.at[idx, "start_letter"] = row["Hchain"]
                df.at[idx, "end_point"] = int(row["L_chain_last_residue"])
                df.at[idx, "end_letter"] = row["Lchain"]

        if row["pdb"] == c_template:
            l_chain_letter = row["Lchain"]
            h_chain_letter = row["Hchain"]

            c_relevant_chains = [chain_id for chain_id in c_cif_chain_order if chain_id in (h_chain_letter, l_chain_letter)]
            print(c_relevant_chains)

            if c_relevant_chains[0] == l_chain_letter:
                df.at[idx, "light_chain_is_first"] = True
                df.at[idx, "start_point"] = int(row["L_chain_first_residue"])
                df.at[idx, "start_letter"] = row["Lchain"]
                df.at[idx, "end_point"] = int(row["H_chain_last_residue"])
                df.at[idx, "end_letter"] = row["Hchain"]
            else:
                df.at[idx, "light_chain_is_first"] = False
                df.at[idx, "start_point"] = int(row["H_chain_first_residue"])
                df.at[idx, "start_letter"] = row["Hchain"]
                df.at[idx, "end_point"] = int(row["L_chain_last_residue"])
                df.at[idx, "end_letter"] = row["Lchain"]
    return df


def extract_gapped_seqs(df, msa_filepath, light_clustal_fname, heavy_clustal_fname):
    """Parses .aln-clustal files with Bio.AlignIO and extracts the gapped sequences."""
    
    # Read clustal files with nested function
    def load_gapped_seqs(filename):
        full_clustal_path = os.path.join(msa_filepath, filename)
        alignment = AlignIO.read(full_clustal_path, "clustal")

        # Dictionary of alignments using a for loop to return {ID: seq...}
        return {record.id.split("|")[0].lower(): str(record.seq) for record in alignment}

    # Load both alignments
    light_seqs = load_gapped_seqs(light_clustal_fname)
    heavy_seqs = load_gapped_seqs(heavy_clustal_fname)
    
    # Assign sequences to gapped_seq_1 and gapped_seq_2 if 
    # light_chain_is_first is True, otherwise 2 then 1
    for idx, row in df.iterrows():
        pdb = row["pdb"].lower()
        if row["light_chain_is_first"] is True:
            df.at[idx, "gapped_seq_1"] = light_seqs.get(pdb, "").replace("\n", "")
            df.at[idx, "gapped_seq_2"] = heavy_seqs.get(pdb, "").replace("\n", "")
        elif row["light_chain_is_first"] is False:
            df.at[idx, "gapped_seq_1"] = heavy_seqs.get(pdb, "").replace("\n", "")
            df.at[idx, "gapped_seq_2"] = light_seqs.get(pdb, "").replace("\n", "")

    # Handle Fab_hybrid row manually by finding the target row
    hybrid_idx = df[df["pdb"].str.lower() == "fab_hybrid"].index
    if not hybrid_idx.empty:
        idx = hybrid_idx[0]
        df.at[idx, "gapped_seq_1"] = light_seqs.get("fab_hybrid", "").replace("\n", "")
        df.at[idx, "gapped_seq_2"] = heavy_seqs.get("fab_hybrid", "").replace("\n", "")

    return df

def write_modeller_pir(
        df: pd.DataFrame, 
        out_path: str,
        v_template: str,
        c_template: str,
        isotype_label: str
        ):
    """Writes a .pir file from a df containing gapped sequences"""
    
    with open(out_path, "w") as pir:
        for _, row in df.iterrows():
            if row["pdb"] == v_template:
                
                v_pir_header = (
                    f">P1;{row['pdb']}\n"
                    f"structureX:{row['pdb']}:{row["start_point"]}:{row["start_letter"]}:{row["end_point"]}:{row["end_letter"]}"
                    f":variable_template:::\n"
                )
                v_gapped_sequence = (f"{row['gapped_seq_1']}/{row['gapped_seq_2']}*\n")

            if row["pdb"] == c_template:

                c_pir_header = (
                    f">P1;{row['pdb']}\n"
                    f"structureX:{row['pdb']}:{row["start_point"]}:{row["start_letter"]}:{row["end_point"]}:{row["end_letter"]}"
                    f":constant_template_{isotype_label}:::\n"
                )
                c_gapped_sequence = (f"{row['gapped_seq_1']}/{row['gapped_seq_2']}*\n")

            if row["template"] == "target":

                target_pir_header = (
                    f">P1;{row['pdb']}_{isotype_label}_target\n"
                    f"sequence:{row['pdb']}_{isotype_label}_target::.::.:hybrid_{isotype_label}_target:::\n"
                )
                target_gapped_sequence = (f"{row['gapped_seq_1']}/{row['gapped_seq_2']}*")

        pir.write(v_pir_header + v_gapped_sequence)
        pir.write(c_pir_header + c_gapped_sequence)
        pir.write(target_pir_header + target_gapped_sequence)

