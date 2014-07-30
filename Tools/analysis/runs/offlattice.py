#!/usr/bin/env python
import math
import os
import numpy as np
import commands

# Note: must be run from within your HemeLB repository.
def run_fab(cmd):
    print commands.getoutput('fab %s' % cmd)

def get_base_name(Re, Vx, x_rotation, y_rotation):
    if (x_rotation, y_rotation) == (0,0):
        base_name = 'cylRe' + ('%g' % Re) + 'VxSz' + ('%.1e' % Vx)
    else:
        base_name =  'cyl' + ('X%iY%i' % (x_rotation, y_rotation)) +'Re' + ('%g' % Re) + 'VxSz' + ('%.1e' % Vx)            
    return base_name

lattices = [15,19,27]
kernels = ['Lbgk', 'Mrt', 'EntA', 'EntC']
# Exclude JYG
bcs = ['Sbb','Fin','Gzs']

resolution = 8e-5
Reynolds = [1]

delta_rotation = 15

for lattice in lattices:
    for kernel in kernels:
        # Skip the non-existent MRT on D3Q27
        if lattice == 27 and kernel == 'Mrt':
            continue
        for bc in bcs:
            # Don't bother rebuilding, already did for Poiseuille and Womersley runs

            for x_rotation in range(delta_rotation,45 + delta_rotation,delta_rotation):
                for y_rotation in range(0, x_rotation + delta_rotation, delta_rotation):
                    for Reynold in Reynolds:
                        base_name = get_base_name(Reynold,resolution,x_rotation,y_rotation)
                        run_fab('hector%d%s%s hemelb:%s,cores=512,images=0'
                            ',snapshots=0,wall_time="1:0:0"' % (lattice,kernel,bc,base_name))
                    pass
                pass
            pass
        run_fab('hector wait_complete')
    pass
pass
