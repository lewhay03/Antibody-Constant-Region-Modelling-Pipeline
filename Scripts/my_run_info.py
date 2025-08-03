# my_run_info
# The purpose of this file is to centralise information about templates, name convensions
# filepaths and outputs.
# Change these as required to guide MODELLER to the correct files.

# my_run_info.py
# This file should be edited by the user before each run of the pipeline

# Template and target identifiers - NOTE move to run_model_pipeline.py
template_v = "1n8z"
template_c = "3m8o"
target_name = "target_IgA1"

# Chain IDs (from the structure files) (these should be extracted after runtime input)
template_v_chain_heavy = "B"
template_v_chain_light = "A"
template_c_chain_heavy = "H"
template_c_chain_light = "L"

# VCAb variables
VCAb_dir = "../VCAb_data/VCAb.csv"
# insert a kappa/lambda light chain preference here

# Expected final VH residues in VH-to-CH1 boundary
vh_boundary = "LVTVSS"

# Directories
fasta_out_dir = "../fasta_sequences"
clustal_out_dir = "../alignments"
pir_out_dir = "../pir_files"
cif_dir = "../atom_files"
models_out_dir = "../models"
pickle_out_dir = "../pickles"

# Output file names (will be suffixed with a timestamp) - change this
FASTA_heavy_out_file = "IgA1_heavy_chain_out.fasta"
FASTA_light_out_file = "IgA1_light_chain_out.fasta"
clustal_heavy_file = "clustal_alignment_IgA1_heavy.clustal"
clustal_light_file = "clustal_alignment_IgA1_light.clustal"
pir_combined_file = "alignment_IgA1_thegoal.pir"

# Utility options
append_timestamp_to_outputs = True
auto_hybridise_template = True
make_alignment_automatically = False  # if True, clustal is triggered automatically (not yet implemented)
create_pir_file = True
store_pickle_output = True

# Use this in the pipeline to pull information from above
from my_run_info import models_out_dir, template_v, template_c, target_name