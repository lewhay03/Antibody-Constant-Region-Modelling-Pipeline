import pandas as pd
import my_run_info

def load_VCAb():
    """Load the raw VCAb.csv database (from my_run_info.py) into a pandas DataFrame."""
    return pd.read_csv(my_run_info.VCAb_dir)

def refine_VCAb(df):
    """Clean and reduce the VCAb DataFrame to important and clean columns only.
    By default, only keeps entries containing kappa light chains 
    and human heavy chains."""

    # Generate clean, readable isotype columns
    # Simplify the naming conventions in certain columns using str.extract(r"regex_expression")
    # This captures the matched word at the start of string, stopping at the first non-word character "("
    df["H_isotype_clean"] = df["Htype"].str.extract(r"^(\w+)")
    df["L_isotype_clean"] = df["Ltype"].str.extract(r"^(\w+)")

    # Keep just the key columns (especially the "_clean" ones)
    columns_to_keep = [
        "pdb","Hchain", "Lchain", "H_seq", "L_seq", "H_coordinate_seq", "L_coordinate_seq", 
        "H_PDB_numbering", "L_PDB_numbering", "pdb_H_VC_Boundary", "pdb_L_VC_Boundary", 
        "title", "release_date","method", "resolution", "carbohydrate", "HC_species", 
        "H_isotype_clean", "L_isotype_clean", "HC_coordinate_seq", 
        "LC_coordinate_seq", "HV_seq", "LV_seq","disulfide_bond"
    ]
    df_keep = df[columns_to_keep]

    # Filter down to rows of interest with conditions. Using wrapped conditions with & 
    # `[(cond_1) & (cond_2) & ...etc]`.
    df_refined_entries = df_keep[
        (df_keep["HC_species"] == "homo_sapiens") & 
        (df_keep["L_isotype_clean"].str.contains("kappa", case=False))
    ]
    return df_refined_entries

# Can be called separately for each template
def zip_template_cif(df: pd.DataFrame, pdb_code: str):
    """Takes in a filtered DataFrame and PDB code. Zips residues with PDB coordinates to their respective PDB numbering, 
    only residues appearing in the x-ray structure are included. Outputs a list of tuples for residues and their numbering."""
    
    # Using iterrows for stable behaviour in pandas
    for idx, row in df.iterrows():

        if row["pdb"] == pdb_code:
            
            # Parse .cif Heavy Chain 
            heavy_cif_residue_list = list(row["H_coordinate_seq"])
            heavy_cif_numbering_list = list(map(int, row["H_PDB_numbering"].split(",")))
            heavy_cif_VC_boundary = int(row["pdb_H_VC_Boundary"])
        
            # Parse .cif Light Chain
            light_cif_residue_list = list(row["L_coordinate_seq"])
            light_cif_numbering_list = list(map(int, row["L_PDB_numbering"].split(",")))
            light_cif_VC_boundary = int(row["pdb_L_VC_Boundary"])

            # Create heavy/light chain zips
            heavy_zip = (list(zip(heavy_cif_residue_list, heavy_cif_numbering_list)))
            light_zip = (list(zip(light_cif_residue_list, light_cif_numbering_list)))
            
            print(f"Zipped light chain length: {len(light_zip)}")
            print(f"Zipped heavy chain length: {len(heavy_zip)}")
            print(f"Light chain VC boundary {light_cif_VC_boundary}th residue")
            print(f"Heavy chain VC boundary {heavy_cif_VC_boundary}th residue")

    return (light_zip, heavy_zip, light_cif_VC_boundary, heavy_cif_VC_boundary)

def get_constant_region(light_zip, heavy_zip, light_boundary, heavy_boundary):
    """Extracts the tuple after the Variable-Constant boundary within a heavy/light chain zip."""

    # Expression format: `(return this for (item1, item2) in iterable if condition)`
    # expression - what we want to return
    # for ... in ... - the loop part
    # if ... - optional condition

    # Find index of the boundary tuple in each chain
    light_boundary_index = next(i for i, (_, num) in enumerate(light_zip) if num == light_boundary)
    heavy_boundary_index = next(i for i, (_, num) in enumerate(heavy_zip) if num == heavy_boundary)

    #print(f"Light Boundary tuple index: {light_boundary_index}")
    #print(f"Heavy Boundary tuple index: {heavy_boundary_index}")

    # Get all tuples after the boundary
    light_post_boundary = light_zip[light_boundary_index + 1:]
    heavy_post_boundary = heavy_zip[heavy_boundary_index + 1:]  

    return (light_post_boundary, heavy_post_boundary)

def get_variable_region(light_zip, heavy_zip, light_boundary, heavy_boundary):
    """Extracts the tuple before the Variable-Constant boundary within a heavy/light chain zip."""

    # Find index of the boundary tuple in each chain
    light_boundary_index = next(i for i, (_, num) in enumerate(light_zip) if num == light_boundary)
    heavy_boundary_index = next(i for i, (_, num) in enumerate(heavy_zip) if num == heavy_boundary)

    #print(f"Light Boundary tuple index: {light_boundary_index}")
    #print(f"Heavy Boundary tuple index: {heavy_boundary_index}")

    # Get all tuples BEFORE AND INCLUDING the boundary
    light_pre_boundary = light_zip[:light_boundary_index + 1]
    heavy_pre_boundary = heavy_zip[:heavy_boundary_index + 1]

    return (light_pre_boundary, heavy_pre_boundary)

# The recombinant antibody must only combine V sequence of v_template + C sequence of c_template
def make_recombinant_seqs(vl_zip, vh_zip, cl_zip, ch_zip):
    """Combines the variable and constant regions of 2 different template sequences.
    Takes 'zipped' (residue-position sequence) chains, joins VL to CL and VH to CH.
    Outputs a tuple of two strings (light_sequence, heavy_sequence)."""

    # Make variable region sequence
    vl_str = ''.join(res for res, _ in vl_zip)
    vh_str = ''.join(res for res, _ in vh_zip)

    # Make constant region sequence
    cl_str = ''.join(res for res, _ in cl_zip)
    ch_str = ''.join(res for res, _ in ch_zip)

    recombinant_seq_light = vl_str + cl_str
    recombinant_seq_heavy = vh_str + ch_str

    print(f"Recombinant Sequence (light): {recombinant_seq_light}")
    print(f"Recombinant Sequence (heavy): {recombinant_seq_heavy}")

    return (recombinant_seq_light, recombinant_seq_heavy)

def make_df_heavies(df_filtered: pd.DataFrame) -> pd.DataFrame:
    """Creates a new dataframe for heavy chain metadata from an already-filtered dataframe."""

    # Call the recombinant sequence function and assign variables
    _, recombinant_seq_heavy = make_recombinant_seqs(vl_zip, vh_zip, cl_zip, ch_zip)

    df_heavies = pd.DataFrame({
    "pdb": df_filtered["pdb"],
    "template": df_filtered["template"],
    "Hchain": df_filtered["Hchain"],
    "H_isotype_clean": df_filtered["H_isotype_clean"],
    "HC_species": df_filtered["HC_species"],
    "H_coordinate_seq": df_filtered["H_coordinate_seq"],
    })

    # Clean Hchain to keep only first character
    df_heavies["Hchain"] = df_heavies["Hchain"].str[0]
    
    first_residue_nos = []
    last_residue_nos = []

    # Using iterrows for stable behaviour in pandas
    for idx, row in df_filtered.iterrows():

        # Split the string of numbering to a list and extract start/end values
        residue_numbering_list = list(map(int, row["H_PDB_numbering"].split(",")))
        last_residue_nos.append(residue_numbering_list[-1])
        first_residue_nos.append(residue_numbering_list[0])
    
    #Add new columns to new df
    df_heavies["H_chain_first_residue"] = first_residue_nos
    df_heavies["H_chain_last_residue"] = last_residue_nos

    # Append recombinant hybrid row to the new dataframe
    hybrid_row = {
        "pdb": "Fab_hybrid",
        "template": "target",
        
        # Assigns the constant region's isotype as the isotype for the target
        "H_isotype_clean": df_filtered.loc[df_filtered["template"] == "c_template", "H_isotype_clean"].values[0],

        # Assigns the species
        "HC_species": df_filtered.loc[df_filtered["template"] == "c_template", "HC_species"].values[0],
        
        # Add the recombinant sequence to the new row dict
        "H_coordinate_seq": recombinant_seq_heavy
    }
    
    df_heavies = pd.concat([df_heavies, pd.DataFrame([hybrid_row])], ignore_index=True)

    return df_heavies

df_heavy = make_df_heavies(df_filtered)
df_heavy