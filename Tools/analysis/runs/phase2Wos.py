#!/usr/bin/env python
import commands

# Note: must be run from within your HemeLB repository.
def run_fab(cmd):
    print commands.getoutput('fab %s' % cmd)

def get_base_name(Re, Vx, Wo):
    base_name = 'cylPhase2Re%gWo%gVx%.1e' % (Re, Wo, Vx)
    return base_name

lattices = [15,19]
kernels = ['Lbgk', 'Mrt']
bcs = ['Sbb', 'Fin', 'Gzs', 'Jyg']

# FOR NOW, leave off the 192-site diameters, in case we're unhappy with these others
resolutions = [12*2**2]
diameter = 8e-3
voxel_sizes = [diameter / x for x in resolutions]

ReWoCombos = [(30,4), (100,8), (300,12)]

for lattice in lattices:
    for kernel in kernels:
        for bc in bcs:
            for (Re,Wo) in ReWoCombos:
                for delta_x in voxel_sizes:
                    base_name = get_base_name(Re,delta_x,Wo)
                    run_fab('hector%d%s%s hemelb:%s,cores=256,images=0'
                       ',snapshots=0,wall_time="12:0:0",label="better_sound_time_ratio_and_longer"' % (lattice,kernel,bc,base_name))
            pass
        run_fab('hector%d%s%s wait_complete' % (lattice,kernel,bc))
    pass
pass
