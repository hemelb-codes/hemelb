#!/usr/bin/env python
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

resolutions = [8e-5 * 100 / x for x in range(30,100,10)]
Reynolds = [10, 100]
x_rot=0
y_rot=0

for lattice in lattices:
    for kernel in kernels:
        # Skip the non-existent MRT on D3Q27
        if lattice == 27 and kernel == 'Mrt':
            continue
        for bc in bcs:
            # Don't bother rebuilding, already did for Poiseuille and Womersley runs
            for Reynold in Reynolds:
                for resolution in resolutions:
                    base_name = get_base_name(Reynold,resolution,x_rot,y_rot)
                    run_fab('hector%d%s%s hemelb:%s,cores=512,images=0'
                       ',snapshots=0,wall_time="1:0:0"' % (lattice,kernel,bc,base_name))
                run_fab('hector%d%s%s wait_complete' % (lattice,kernel,bc))
            pass
        pass
    pass
pass
