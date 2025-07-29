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
        "pdb","Hchain", "Lchain", "H_coordinate_seq", "L_coordinate_seq", 
        "title", "release_date","method", "resolution", "carbohydrate", 
        "HC_species", "H_isotype_clean", "L_isotype_clean", "HC_coordinate_seq", 
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

