# This script runs the whole pipeline when executed through MODELLER

# === User Inputs ===
import os
import my_run_info
from datetime import datetime

# Create a date-time stamp for the run that can be added to file names
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Set the project root using absolute path and adding ".." to 
# move up 1 level from .../Modelling-Pipeline/Scripts -> .../Modelling-Pipeline
project_root = os.path.abspath(os.path.join(os.getcwd(), ".."))

# List of required directories
required_dirs = [
    "VCAb_data",
    "fasta_sequences",
    "alignments",
    "pir_files",
    "atom_files",
    "models",
    "pickles"
]

# Create each directory if it doesn't exist
for dir_name in required_dirs:
    os.makedirs(os.path.join(project_root, dir_name), exist_ok=True)

# Loop inputs until matching files are found
def prompt_for_existing_file(prompt, suffix="cif", search_dir="../VCAb_data"):
    while True:
        pdb_id = input(prompt).strip().lower()
        file_path = os.path.join(search_dir, f"{pdb_id}{suffix}")
        if os.path.isfile(file_path):
            return pdb_id
        print(f"File '{file_path}' not found. Please try again.")

# Define the templates of interest
v_template = prompt_for_existing_file(
    "Enter PDB code for variable region template (e.g., 1n8z): "
    )
c_template = prompt_for_existing_file(
    "Enter PDB code for Fab constant region template (e.g., 3m8o): "
    )
#c_Fc_template = prompt_for_existing_file(
#    "Enter PDB code for Fc constant region template (e.g., XXXX): "
#    ) # Potentially deactivate unless my_run_info states True
#predicted_hinge_template = prompt_for_existing_file(
#    "Enter filename for predicted hinge template (e.g., AF_NNNN_hinge): "
#    ) # Potentially deactivate unless my_run_info states True



# === Sequence Preparation ===
import prepare_sequences

# Load VCAb into dataframe
df_raw = prepare_sequences.load_VCAb()

# Refine & filter VCAb dataframe
df_refined = prepare_sequences.refine_VCAb(df_raw)

# Visualise the filtered VCAb dataframe
print(f"Rows, columns: {df_refined.shape}")
df_refined.head(5)

# Filter only rows where matching user-entered pdb is True and create a copy of the dataframe. 
# NOTE: `df_matches` is an independant copy of `df_refined`. Use df_matches from this point on
df_matches = df_refined[
    df_refined["pdb"].isin([v_template, c_template])
].copy()

# Add a 'template' annotation column
df_matches.insert(loc=1, column="template", value=None)
df_matches.loc[df_matches["pdb"] == v_template, "template"] = "v_template"
df_matches.loc[df_matches["pdb"] == c_template, "template"] = "c_template"

df_matches

# Extract Immunoglobulin isotype label and store as variable
# NOTE: add this to the write FASTA function
for idx, row in df_matches.iterrows():
    if row["pdb"] == c_template:
        isotype_label = str(row["H_isotype_clean"])
        
print(f"Constant region template is {isotype_label}.")

# Do you wish to create FASTA files?
    # User input [y/n]. If n, skip to next section
    
    # Extract sequences and trim overlapping regions

    # Write heavy and light chain FASTA files

# Uses * to unpack the tuple from zip function
cl_zip, ch_zip = prepare_sequences.get_constant_region(
    *prepare_sequences.zip_template_cif(df_matches, c_template)
    )
vl_zip, vh_zip = prepare_sequences.get_variable_region(
    *prepare_sequences.zip_template_cif(df_matches, v_template)
    )

# Make recombinant antibody Fab region sequences
recombinant_seq_light, recombinant_seq_heavy = prepare_sequences.make_recombinant_seqs(vl_zip, vh_zip, cl_zip, ch_zip)

# Create 1 dataframe for heavy chains and 1 for light chains
df_heavy = prepare_sequences.make_df_heavies(df_matches, recombinant_seq_heavy)
df_light = prepare_sequences.make_df_lights(df_matches, recombinant_seq_light)

df_light
df_heavy

# Write 1 heavy and 1 light fasta file for submission to clustal aligner
prepare_sequences.write_fastas(df_light, df_heavy, 
                               my_run_info.fasta_out_dir, 
                               df_matches,
                               timestamp
                               )

# === Clustal Alignment ===
# NOTE: Not currently set up

proceed_clustal = input("""Pause here:
    1. Create heavy and light alignment files with clustal 
    2. Rename files appropriately
    3. Press 'y' to continue: """)
if proceed_clustal.lower() != 'y':
    print("Exiting... Please run again when ready.")
    exit()

# === .pir File Creation ===

import convert_to_pir

# Combine the light and heavy information into one table for .pir
df_combined = convert_to_pir.merge_df_for_pir(df_light, df_heavy)

# Determine the order that chains appear in the cif file
v_cif_chain_order = convert_to_pir.cif_parse(v_template, f"{v_template}.cif", my_run_info.cif_dir)
c_cif_chain_order = convert_to_pir.cif_parse(c_template, f"{c_template}.cif", my_run_info.cif_dir)

# Extract relevant chain info and add to table
df_combined_populated = convert_to_pir.relevant_chains(
    df_combined, 
    v_cif_chain_order, 
    c_cif_chain_order, 
    v_template, 
    c_template
)

# Parses .aln-clustal files with Bio.AlignIO and extracts the gapped sequences.
df_for_pir = convert_to_pir.extract_gapped_seqs(
    df_combined_populated,
    my_run_info.clustal_out_dir,
    light_clustal_fname=f"Fab_alignment_{isotype_label}_light.aln-clustal", # replace this with automation
    heavy_clustal_fname=f"Fab_alignment_{isotype_label}_heavy.aln-clustal" # replace this with automation
)

# Write .pir file populated to the standard required by MODELLER
# The allowed format:
#   >P1;3m8o
#   structure:pdb_file:.:.:.:.::::
#   seq1---/seq2---*

convert_to_pir.write_modeller_pir(
    df_for_pir, 
    f"../pir_files/pir_alignment_{isotype_label}_{timestamp}.pir",
    v_template,
    c_template,
    isotype_label
    )