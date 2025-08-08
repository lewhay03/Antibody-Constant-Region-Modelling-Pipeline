# Addition of restraints to the default ones
from modeller import *
from modeller.automodel import *    # Load the AutoModel class

log.verbose()
env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

class MyModel(AutoModel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
#       Add some restraints from a file:
#       rsr.append(file='my_rsrs1.rsr')

#       Residues 20 through 30 should be an alpha helix:
        rsr.add(secondary_structure.Alpha(self.residue_range('20:A', '30:A')))
#       Two beta-strands:
        rsr.add(secondary_structure.Strand(self.residue_range('1:A', '6:A')))
        rsr.add(secondary_structure.Strand(self.residue_range('9:A', '14:A')))
#       An anti-parallel sheet composed of the two strands:
        rsr.add(secondary_structure.Sheet(at['N:1:A'], at['O:14:A'],
                                          sheet_h_bonds=-5))
#       Use the following instead for a *parallel* sheet:
#       rsr.add(secondary_structure.Sheet(at['N:1:A'], at['O:9:A'],
#                                         sheet_h_bonds=5))

#       Restrain the specified CA-CA distance to 10 angstroms (st. dev.=0.1)
#       Use a harmonic potential and X-Y distance group.
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CA:547:B'], # Pro 329 + 218 = 547
                                                         at['CA:1178:G']), # Ser 298 +662 +218 = 516
                               mean=5.0, stdev=0.1))
    
    def special_patches(self, aln):
        # A disulfide between residues 444 and 1106 in chain B & D:
        self.patch(residue_type='DISU', residues=(self.residues['444:B'],
                                                  self.residues['1106:D']))
        # A disulfide between residues 447 and 1109 in chain B & D:
        self.patch(residue_type='DISU', residues=(self.residues['444:B'],
                                                  self.residues['1109:D']))

a = MyModel(env,
    alnfile  = '../pir_files/pir-alignment-5dk3-monomer-with-AF-hinge-dimertest2.pir', # alignment filename NOTE: import this from other script
    knowns   = ('IgG4ModellerHingeTemplate'),                          # codes of the templates
    sequence = 'IgG4FullAbTarget',                             # code of the target (used by MODELLER to name files)
    assess_methods=(assess.DOPE, assess.GA341)      # assessment methods
    )              # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 5                 # index of the last model
                                    # (determines how many models to calculate)

import os
os.chdir("../models")               # change output directory (this may need changing back to ../Scripts each run)
a.set_output_model_format("MMCIF")  # request mmCIF rather than PDB outputs

a.initial_malign3d = True           # superpose the 3d structures before modelling

a.make()                            # do comparative modeling
