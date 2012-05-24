#!/usr/bin/env python
# encoding: utf-8

import unittest
import subprocess
import os
import shutil
import tempfile
from hemeTools.parsers.extraction import ExtractedProperty

class TestPoiseuilleFlowTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.viscosity = 4e-3 # HemeLB's default dynamic viscosity (Pa s)
        self.pressure_diff = 16 * 133.3223874 # 16mmHg in Pa
        self.pipe_length = 6e-2 # 60mm in m

        self.temp_dir = tempfile.mkdtemp("_HemeLB_RegressionTest")
        shutil.copy("../../build/hemelb", self.temp_dir)
        shutil.copy("resources/poiseuille_flow_test.gmy", self.temp_dir)
        shutil.copy("resources/poiseuille_flow_test.xml", self.temp_dir)
        os.chdir(self.temp_dir)


    def test_run_simulation_and_check_output_created(self):
        try:
            subprocess.call("mpirun -np 4 hemelb -in poiseuille_flow_test.xml -i 10 -s 1", shell=True)
        except OSError, e:
            print >>sys.stderr, "Call to HemeLB failed:", e
            
        # Make sure the .dat files have been created
        os.path.isfile("results/Extracted/velocity_40mm_in.dat")
        os.path.isfile("results/Extracted/shear_stress_40mm_in.dat")

    def compute_analytical_velocity(self, site_data):
        x_coord = site_data[0]
        dist_centre = abs(self.centre - x_coord)
        analytical_vel = (1/(4*self.viscosity)) * (self.pressure_diff/self.pipe_length) * (pow(self.radius,2) - pow(dist_centre,2))
        return analytical_vel

    def test_velocity_profile(self):
        filename = "results/Extracted/velocity_40mm_in.dat"
        propFile = ExtractedProperty(filename)

        print '# Dump of file "{}"'.format(filename)
        print '# File has {} sites.'.format(propFile.siteCount)
        print '# File has {} fields:'.format(propFile.fieldCount)
        for name, xdrType, memType, length, offset in propFile._fieldSpec:
            print '#     "{0}", length {1}'.format(name, length)
        print '# Geometry origin = {} m'.format(propFile.originMetres)
        print '# Voxel size = {} m'.format(propFile.voxelSizeMetres)

        header  = '# '+ ', '.join(name for name, xdrType, memType, length, offset in propFile._fieldSpec)
        print header

        # Property extraction files could have info for more than one time step depending on the frequency requested
        for t in propFile.times:
            sites_along_line = propFile.GetByTimeStep(t)
            print "# Timestep {:d}".format(t)
            
            # Create a list of tuples (x_coordinate, z_velocity) for all the sites along the line
            velocities_along_line = []
            [velocities_along_line.append((site.position[0], site.velocity_40mm_in[2])) for site in sites_along_line]

            # Work out pipe radius
            max_coord = max(velocities_along_line)[0]
            min_coord = min(velocities_along_line)[0]

            # With bounce back, the actual wall is half a lattice site away
            self.radius = (max_coord - min_coord)/2 + propFile.voxelSizeMetres/2
            self.centre = (max_coord + min_coord)/2

            # Compare velocities along the line with analytical solution
            analytical_solutions = map(self.compute_analytical_velocity, velocities_along_line)

            # Check the results of the comparison
            [self.assertAlmostEqual(analytical, computed, delta=1e-3, msg="Velocity {0} differs from analytical solution {1} at site with z coordinate {2}".format(computed,analytical,z_coord))
             for (analytical, (z_coord,computed)) in zip(analytical_solutions, velocities_along_line)]

    def test_shear_stress_profile(self):
        pass
    
