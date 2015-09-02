#!/usr/bin/env python
import os.path
import numpy as np
import quantities as pq
from xml.etree import ElementTree
import shutil
import glob
import time

import _fixpath

from .HemeLbParameters import HemeLbParameters
from pvs.helpers import IndentXml, cd
from pvs.simplify import simplify

import pdb

class CylinderGeneration(HemeLbParameters):
    def __init__(self, reynolds_number, diameter_voxels, tau):
        self.ReynoldsNumber = reynolds_number * pq.dimensionless
        self._DiameterVoxels = diameter_voxels
        self._Tau = tau
        return

    @property
    def Tau(self):
        return self._Tau
    
    @property
    def VoxelSize(self):
        return self.Diameter / self._DiameterVoxels
    
    @property
    @simplify
    def TimeStep(self):
        kin_visc_lat = (self.Tau - 0.5) / 3.
        return kin_visc_lat * self.VoxelSize**2 / self.KinematicViscosity
    
    @property
    def Diameter(self):
        return 0.008 * pq.metre

    @property
    def Radius(self):
        return 0.5 * self.Diameter
    
    @property
    def Length(self):
        return 0.032 * pq.metre

    @property
    @simplify
    def MaximumVelocity(self):
        return self.ReynoldsNumber * self.KinematicViscosity / self.Diameter

    @property
    @simplify
    def MachNumber(self):
        return self.MaximumVelocity / self.SpeedOfSound

    @property
    @simplify
    def PressureGradient(self):
        return 4 * self.MaximumVelocity * self.DynamicViscosity / self.Radius**2

    @property
    @simplify
    def PressureDifference(self):
        return self.PressureGradient * self.Length

    @property
    @simplify
    def MomentumDiffusionTime(self):
        return self.Diameter**2 / self.KinematicViscosity

    @property
    @simplify
    def SoundTime(self):
        return self.Length / self.SpeedOfSound
    
    @property
    def NumberOfTimeSteps(self):
        nSteps = self.SimulationTime / self.TimeStep
        return int(nSteps.simplified)

    @property
    def TotalSites(self):
        return self.ScaleToLatticeUnits(np.pi * self.Radius**2 * self.Length)

    @property
    def Axis(self):
        approx = np.array((-0.299, 0.382, 0.874))
        return pq.dimensionless * approx/np.sqrt(np.dot(approx,approx))
    
    _OutputUnitsMap = {'pressure': 'mmHg'}

    @property
    def InletPressure(self):
        return self.ReferencePressure + self.PressureDifference
    
    @property
    def OutletPressure(self):
        return self.ReferencePressure 

    @property
    def InletPosition(self):
        return 0.5 * self.Length * self.Axis
    
    @property
    def OutletPosition(self):
        return -0.5 * self.Length * self.Axis
        
    @property
    def OutputGeometryFile(self):
        return 'config.gmy'
    @property
    def OutputXmlFile(self):
        return 'config.xml'
    @property
    def OutputSummaryFile(self):
        return 'config.smy'
    
    def Generate(self, outdir):
        def ConvertPressure(pMean):
            pMean_mmHg = float(pMean / self.mmHg)
            return Vector.Vector(pMean_mmHg, 0., 0.)
        from HemeLbSetupTool.Model import OutputGeneration, Vector

        rundir = os.path.join(outdir, self.RunName)
        os.mkdir(rundir)
        with cd(rundir):
            gen = OutputGeneration.CylinderGenerator(
                self.OutputGeometryFile,
                self.OutputXmlFile,
                self.ScaleForOutput(self.VoxelSize),
                self.ScaleForOutput(self.Axis),
                self.ScaleForOutput(self.Length),
                self.ScaleForOutput(self.Radius),
                InletPressure=ConvertPressure(self.InletPressure),
                OutletPressure=ConvertPressure(self.OutletPressure)
                )
            gen.Execute()
            self.FixXml()
            self.WriteSummary()
            self.WriteSubmitScript()
            
    @property
    def SubmitScript(self):
        return 'sub_for_exe.sh'
    
    def WriteSubmitScript(self):
        script = '''#!/bin/bash
cd $HEMELB
{0.SubmitCommand},executable=$1
'''.format(self)
        with file(self.SubmitScript, 'w') as f:
            f.write(script)
        os.chmod(self.SubmitScript, 0755)
        return
            
    def FixXml(self):
        tree = ElementTree.parse(self.OutputXmlFile)
        root = tree.getroot()
        simulation = root.find('simulation')

        simulation.attrib.pop('cycles', None)
        simulation.attrib.pop('cyclesteps', None)
        simulation.attrib['steps'] = str(self.NumberOfTimeSteps)
        simulation.attrib['step_length'] = str(self.ScaleForOutput(self.TimeStep))
        simulation.attrib['extra_warmup_steps'] = str(self.WarmUpSteps)
        
        prop = root.find('properties')
        while prop is not None:
            root.remove(prop)
            prop = root.find('properties')

        self.MakeInitialConditionsElement(tree.getroot())
        self.MakeExtractionElement(tree.getroot())
        self.RewriteInlets(root)
        
        IndentXml(tree.getroot())
        tree.write(self.OutputXmlFile)
        return

    def MakeInitialConditionsElement(self, root):
        initConds = ElementTree.SubElement(root, 'initialconditions')
        pressure = ElementTree.SubElement(initConds, 'pressure')
        uniform = ElementTree.SubElement(
            pressure, 'uniform',
            {'value': str(0.0),
             'units': 'mmHg'}
            )
    
    def MakeExtractionElement(self, root):
        def VecToDict(v):
            v = self.ScaleForOutput(v)
            return dict((d, str(v[i])) for i, d in enumerate(('x', 'y', 'z')))
        
        extr = ElementTree.SubElement(root, 'properties')
        plane = ElementTree.SubElement(
            extr, 'propertyoutput',
            {'frequency': str(self.ExtractionInterval),
             'file': 'plane.xtr'}
            )
        geom = ElementTree.SubElement(plane, 'planegeometry')
        ElementTree.SubElement(geom, 'point', VecToDict(pq.metre * (0,0,0)))
        ElementTree.SubElement(geom, 'normal', VecToDict(self.Axis))
        ElementTree.SubElement(geom, 'field', type='velocity')

        line = ElementTree.SubElement(
            extr, 'propertyoutput',
            {'frequency': str(self.ExtractionInterval),
             'file': 'line.xtr'}
             )
        geom = ElementTree.SubElement(line, 'linegeometry')
        ElementTree.SubElement(geom, 'point', VecToDict(self.InletPosition))
        ElementTree.SubElement(geom, 'point', VecToDict(self.OutletPosition))
        ElementTree.SubElement(geom, 'field', type='pressure')        
        return

    @property
    def HectorCores(self):
        basic = int(self.TotalSites / 5000) + 1
        if basic < 32:
            return max(basic, 4)
        nodes = basic / 32
        return nodes * 32

    @property
    def HectorSeconds(self):
        sups = self.TotalSites * self.NumberOfTimeSteps
        time = sups / (1e6 * self.HectorCores)
        return 2 * time + 300

    @property
    def HectorWallTime(self):
        return time.strftime('%H:%M:%S', time.gmtime(self.HectorSeconds))

    @property
    def SubmitCommand(self):
        return 'fab hector hemelb:config={0.RunName},images=0,cores={0.HectorCores},wall_time={0.HectorWallTime}'.format(self)
    
    @property
    def RunInfo(self):
        lines = []
        for attr in ['ReynoldsNumber', 'Diameter', 'HectorCores', 'HectorSeconds']:
            val = getattr(self, attr)
            if isinstance(val, pq.Quantity):
                val = self.ScaleToLatticeUnits(val)
            
            lines.append(attr + ' = ' + str(val))
        
        return '\n'.join(lines)
