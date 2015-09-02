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
bcs = ['Fin', 'Gzs', 'Jyg']

# Only do finest in a different script
resolutions = [12*2**n for n in range(0, 5)]
cores = [128, 128, 256, 2048, 16384]
diameter = 8e-3
voxel_sizes = [diameter / x for x in resolutions]

Reynolds = [1, 30, 100]
Wos = [0]

for lattice in lattices:
    for kernel in kernels:
        # Skip the non-existent MRT on D3Q27
        if lattice == 27 and kernel == 'Mrt':
            continue
        for bc in bcs:
            for Reynold in Reynolds:
                for Wo in Wos:
                    for (delta_x,core_count) in zip(voxel_sizes, cores):
                        base_name = get_base_name(Reynold,delta_x,Wo)
                        run_fab('hector%d%s%s hemelb:%s,cores=%i,images=0'
                           ',snapshots=0,wall_time="3:0:0"' % (lattice,kernel,bc,base_name,core_count))
                run_fab('hector%d%s%s wait_complete' % (lattice,kernel,bc))
            pass
        pass
    pass
pass
