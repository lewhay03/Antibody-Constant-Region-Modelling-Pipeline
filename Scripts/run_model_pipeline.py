# This script runs the whole pipeline when executed through MODELLER

# === User Inputs ===
import os

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

# Define the templates of interest
v_template = input("Enter PDB code for variable region template (e.g., 1n8z): ").strip().lower()
c_template = input("Enter PDB code for Fab constant region template (e.g., 3m8o): ").strip().lower()
c_Fc_template = input("Enter PDB code for Fc constant region template (e.g., XXXX): ").strip().lower() # Potentially deactivate unless my_run_info states True
predicted_hinge_template = input("Enter filename for predicted hinge template (e.g., AF_NNNN_hinge): ") # Potentially deactivate unless my_run_info states True

# Check the pdb code matches a file name in the folder
if os.path.exists(f"../atom_files/{v_template}.cif") & os.path.exists(f"../atom_files/{c_template}.cif"): # fix this conditional and make it loop back or terminate if no value entered
    print("Files found.")
else:
    print(f"File for {v_template} not found")

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
cl_zip, ch_zip = prepare_sequences.get_constant_region(*prepare_sequences.zip_template_cif(df_matches, c_template))
vl_zip, vh_zip = prepare_sequences.get_variable_region(*prepare_sequences.zip_template_cif(df_matches, v_template))

# Make recombinant antibody Fab region
recombinant_seq_light, recombinant_seq_heavy = prepare_sequences.make_recombinant_seqs(vl_zip, vh_zip, cl_zip, ch_zip)

# === Clustal Alignment ===
# NOTE: Not currently set up

# === .pir File Creation ===
