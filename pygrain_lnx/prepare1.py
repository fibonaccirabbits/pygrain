#converts a prot.pdb(atomic) to pro_cg.pdb (coarse grained)
#also creates other files: .itp, .top, .ssd, .gro
#uses martini and insane

import pygrain

cg = pygrain.Cg()
cg.prepare_martini_insane()

print '\n DONE!!! \n'
