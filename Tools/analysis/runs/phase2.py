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

# Only do finest in a different script
resolutions = [12*2**n for n in range(0, 5)]
diameter = 8e-3
voxel_sizes = [diameter / x for x in resolutions]

Reynolds = [100]
Wos = [0]

for lattice in lattices:
    for kernel in kernels:
        # Skip the non-existent MRT on D3Q27
        if lattice == 27 and kernel == 'Mrt':
            continue
        for bc in bcs:
            for Reynold in Reynolds:
                for Wo in Wos:
                    for delta_x in voxel_sizes:
                        base_name = get_base_name(Reynold,delta_x,Wo)
                        run_fab('hector%d%s%s hemelb:%s,cores=2048,images=0'
                           ',snapshots=0,wall_time="1:0:0"' % (lattice,kernel,bc,base_name))
            pass
            run_fab('hector%d%s%s wait_complete' % (lattice,kernel,bc))
        pass
    pass
pass
