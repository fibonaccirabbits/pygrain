#performs minimizations and equilibrations.
#uses GROMACS
#uses the new gmx syntax

import pygrain
cg = pygrain.Cg()
cg.gmx_minim_eq()

print '\n DONE!!! \n'
