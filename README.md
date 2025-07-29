# Antibody-Constant-Region-Modelling-Pipeline-Research-Project

This repository contains code and data for building and evaluating structural models of antibody constant regions 
from different isotypes, combined with the trastuzumab variable region (VH).

## Project Goals
- Generate and compare Fab models with the 9 human isotype CH1 domains in the first instance
- Generate and compare full antibody models with hinge regions modelled and decorated
- Use MODELLER, Clustal Omega, Biopython, and PDB data
- Develop a modular and reproducible modeling pipeline in Python


## Current Status
- VCAb database downloaded
- Isotype Fab template structures selected
- FASTA and PDB templates collected
- Workspace and repo initialized
- Initial script underway
- Directory structure finalised

## Next Steps
- Automate MSA and PIR generation
- Build batch model generation for 9 isotype hybrids
- Add scoring and model evaluation

## Prerequisites
- MODELLER 10.6 or later installed on machine
- Python environment of 3.9 or later is set up
- Set working directory to the location of this file
- Visit https://fraternalilab.cs.ucl.ac.uk/VCAb/ and download a .csv copy of "all entries in VCAb" from 
the download tab, saving it within "./VCAb_data"
- Save a copy of your desired .PDBx/.cif template files to "./PDB_data" with their 4-digit PDB code 
as the file name

## Instructions for use
- Open modeller exe to bring up commandline interface
- Type `python run_model_pipeline.py` and Enter
- When prompted input your template information
- The pipeline generates all necessary files for MODELLER which will generate time-stamped models