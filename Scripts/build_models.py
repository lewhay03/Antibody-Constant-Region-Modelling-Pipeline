# Comparative modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the AutoModel class

log.verbose()    # request verbose output

# space reserved here for defining custom parameters `class MyModel(AutoModel):`

env = Environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = AutoModel(env,
              alnfile  = '../PIR_files/alignment_IgA1_thegoal.pir', # alignment filename (NOTE: relative filepath)
              knowns   = ('1n8z', '3m8o'),                          # codes of the templates
              sequence = 'target_IgA1',                             # code of the target (used by MODELLER to name files)
              assess_methods=(assess.DOPE, assess.GA341))           # assessment methods
a.starting_model= 1                 # index of the first model
a.ending_model  = 5                 # index of the last model
                                    # (determines how many models to calculate)

import os
os.chdir("../Models")               # change output directory (this may need changing back to ../Scripts each run)
a.set_output_model_format("MMCIF")  # request mmCIF rather than PDB outputs
a.name = "VC_Hybrid_IgA1_Model_Fab" # change the output filename prefix
a.initial_malign3d = True           # superpose the 3d structures (good for multi chain)
a.make()                            # do the actual comparative modeling

# Get a list of all successfully built models from a.outputs
ok_models = [x for x in a.outputs if x['failure'] is None]

# Rank the models by DOPE score
key = 'DOPE score'
ok_models.sort(key=lambda a: a[key])

# Top Model
top_model = ok_models[0]
print("Top model: %s (DOPE score %.3f)" % (top_model['name'], top_model[key]))

# Preserve the model output (a.output) with pickle
import pickle
this_temp_v = "1n8z"
this_temp_c = "3m8o"
this_target = "target_IgA1"
pickle_counter = 0


# Loop until file is written
while True:
    
    # Set this run's pickle file name
    this_run_pickle_filename = f"model_outputs_{this_temp_v}_{this_temp_c}_{this_target}{pickle_counter}.pkl"
    
    # Avoid overwriting if the name already exists using a counter
    if not os.path.exists(this_run_pickle_filename):
        break
    print("File name already exists. Adding a suffix digit.")
    pickle_counter =+ 1
    

# Write binary pickle file
with open(this_run_pickle_filename, "wb") as f:
    pickle.dump(a.outputs, f)