#!/usr/bin/env python
import commands

# Note: must be run from within your HemeLB repository.
def run_fab(cmd):
    print commands.getoutput('fab %s' % cmd)

def get_base_name(Re, Vx, Wo):
    base_name = 'cylPhase2Re%gWo%gVx%.1e' % (Re, Wo, Vx)
    return base_name

lattices = [15,19]
kernels = ['Lbgk']
bcs = ['Sbb']

# Just do the 192-site diameters
resolutions = [12*2**4]
diameter = 8e-3
voxel_sizes = [diameter / x for x in resolutions]

Reynolds = [1, 30]
Wos = [0]

for lattice in lattices:
    for kernel in kernels:
        # Skip the non-existent MRT on D3Q27
        if lattice == 27 and kernel == 'Mrt':
            continue
        for bc in bcs:
            # No rebuild
            for Reynold in Reynolds:
                for Wo in Wos:
                    for delta_x in voxel_sizes:
                        base_name = get_base_name(Reynold,delta_x,Wo)
                        run_fab('hector%d%s%s hemelb:%s,cores=2048,images=0'
                           ',snapshots=0,wall_time="3:0:0"' % (lattice,kernel,bc,base_name))
            pass
        pass
    pass
pass
