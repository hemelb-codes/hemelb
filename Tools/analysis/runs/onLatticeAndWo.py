#!/usr/bin/env python
import commands
import os
import math

# Note: must be run from within your HemeLB repository.
def run_fab(cmd):
    print commands.getoutput('fab %s' % cmd)

def GetBaseName(Re, Vx, Wo = None):
    if (Wo == None):
        return 'cylRe' + ('%g' % Re) + 'VxSz' + ('%.1e' % Vx)
    return 'cylRe' + ('%g' % Re) + 'VxSz' + ('%.1e' % Vx) + 'Wo' + ('%g' % Wo)

lattices = [15,19,27]
kernels = ['Lbgk', 'Mrt', 'EntA', 'EntC']
bcs = ['Sbb','Fin','Gzs','Jyg']

resolutions=[8e-5]
Reynolds=[0.1,1,10,100,1000]
womersley_list=[4, 8, 12, 16, 20]

for lattice in lattices:
    for kernel in kernels:
        # Skip the non-existent MRT on D3Q27
        if lattice == 27 and kernel == 'Mrt':
            continue

        for bc in bcs:
            # Build first then submit all the runs
            run_fab('hector%d%s%s sync require_recopy batch_build:no_streaklines,HEMELB_BUILD_UNITTESTS=OFF,no_debug,HEMELB_LOG_LEVEL=Debug wait_complete' %
                (lattice,kernel,bc))

            for resolution in resolutions:
                for Reynold in Reynolds:
                    time = "1:0:0" if Reynold == 1000 else "0:20:0"

                    # On-lattice, Poiseuille runs
                    base_name = GetBaseName(Reynold,resolution)
                    run_fab('hector%d%s%s hemelb:%s,cores=2048,images=0'
                         ',snapshots=0,wall_time="%s"' % (lattice,kernel,bc,base_name,time))

                    # Only do Wo for 10 and 100
                    if Reynold not in [10,100]:
                        continue

                    # On-lattice, Womersley runs
                    for womersley in womersley_list:
                        base_name = GetBaseName(Reynold, resolution, womersley)
                        run_fab('hector%d%s%s hemelb:%s,cores=2048,images=0'
                                ',snapshots=0,wall_time="%s"' % (lattice,kernel,bc,base_name,time))

                    # Ensure the several jobs just queued are finished before we continue
                    run_fab('hector%d%s%s wait_complete' % (lattice,kernel,bc))
            pass
        pass
    pass
pass
