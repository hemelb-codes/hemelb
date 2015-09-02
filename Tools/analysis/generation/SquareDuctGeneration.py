#!/usr/bin/env python
import os.path
import numpy as np
import quantities as pq
from xml.etree import ElementTree
import shutil
import glob

import _fixpath
from .HemeLbParameters import HemeLbParameters
from pvs.helpers import IndentXml, cd

class SquareDuctParameters(HemeLbParameters):
    def __init__(self, reynolds_number, side_length_voxels, tau):
        self.ReynoldsNumber = reynolds_number * pq.dimensionless
        self._SideLengthVoxels = side_length_voxels * pq.dimensionless
        self._Tau = tau
        return

    @property
    def Tau(self):
        return self._Tau
    
    @property
    def VoxelSize(self):
        return self.SideLength / self._SideLengthVoxels
    
    @property
    def TimeStep(self):
        kin_visc_lat = (self.Tau - 0.5) / 3.
        return kin_visc_lat * self.VoxelSize**2 / self.KinematicViscosity
    
    @property
    def SideLength(self):
        return 0.008 * pq.metre

    @property
    def DuctLength(self):
        return 0.016 * pq.metre

    @property
    def MaximumVelocity(self):
        return (self.ReynoldsNumber * self.KinematicViscosity / self.SideLength).simplified

    @property
    def MachNumber(self):
        return (self.MaximumVelocity / self.SpeedOfSound).simplified

    @property
    def PressureGradient(self):
        return (self.MaximumVelocity * self.DynamicViscosity / (self.SideLength**2 * 0.0736714)).simplified

    @property
    def PressureDifference(self):
        return self.PressureGradient * self.DuctLength

    @property
    def MomentumDiffusionTime(self):
        return self.SideLength**2 / self.KinematicViscosity

    @property
    def SimulationTime(self):
        return 4 * self.MomentumDiffusionTime
    
    @property
    def NumberOfTimeSteps(self):
        nSteps = self.SimulationTime / self.TimeStep
        return int(nSteps.simplified)

    @property
    def TotalSites(self):
        return self.ScaleToLatticeUnits(self.SideLength**2 * self.DuctLength)
    
    pass

class SquareDuctGeneration(SquareDuctParameters):
    _OutputUnitsMap = {'pressure': 'mmHg'}

    @property
    def InletPressure(self):
        return self.ReferencePressure + self.PressureDifference
    
    @property
    def OutletPressure(self):
        return self.ReferencePressure
    
    @property
    def RunName(self):
        return 'SquareDuct_Re{}_L{}'.format(self.ScaleForOutput(self.ReynoldsNumber),
                                            int(self.ScaleToLatticeUnits(self.SideLength)))
    
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
            gen = OutputGeneration.SquareDuctGenerator(
                self.OutputGeometryFile,
                self.OutputXmlFile,
                float((self.VoxelSize / pq.metre).simplified),
                2,
                int(np.rint((self.DuctLength / self.VoxelSize).simplified)),
                int(np.rint((self.SideLength / self.VoxelSize).simplified)),
                InletPressure=ConvertPressure(self.InletPressure),
                OutletPressure=ConvertPressure(self.OutletPressure)
                )
            gen.Execute()
            self.FixXml()
            self.WriteSummary()

    def WriteSummary(self):
        import yaml
        data = {}
        data['Type'] = 'SquareDuct'

        for attr in ['ReynoldsNumber', '_SideLengthVoxels', 'Tau']:
            val = getattr(self, attr)
            if isinstance(val, pq.Quantity):
                val = self.ScaleForOutput(val)
            data[attr] = val

        with file(self.OutputSummaryFile, 'w') as f:
            yaml.dump(data, f)
        return
    
    @classmethod
    def LoadFromSummary(cls, summary):
        import yaml
        with file(summary) as f:
            data = yaml.load(f)

        assert data['Type'] == 'SquareDuct'
        
        args = []
        for attr in ['ReynoldsNumber', '_SideLengthVoxels', 'Tau']:
            args.append(data[attr])
        return cls(*args)
    
    # def DuplicateXml(self):
    #     xmls = glob.glob('SquareDuct_Re*_L{}/config.xml'.format(int(self.ScaleToLatticeUnits(self.SideLength))))
    #     if xmls[0] != self.OutputXmlFile:
    #         shutil.copyfile(xmls[0], self.OutputXmlFile)
    #     return
    
    def FixXml(self):
        tree = ElementTree.parse(self.OutputXmlFile)
        root = tree.getroot()
        simulation = root.find('simulation')

        simulation.attrib.pop('cycles', None)
        simulation.attrib.pop('cyclesteps', None)
        simulation.attrib['steps'] = str(self.NumberOfTimeSteps)
        simulation.attrib['step_length'] = str(self.ScaleForOutput(self.TimeStep))

        prop = root.find('properties')
        while prop is not None:
            root.remove(prop)
            prop = root.find('properties')
        
        self.MakeExtractionElement(tree.getroot())
        
        IndentXml(tree.getroot())
        tree.write(self.OutputXmlFile)
        return
        
    def MakeExtractionElement(self, root):
        extr = ElementTree.SubElement(root, 'properties')
        plane = ElementTree.SubElement(
            extr, 'propertyoutput',
            {'frequency': str(int(self.ScaleToLatticeUnits(self.MomentumDiffusionTime/4))),
             'file': 'plane.xtr'}
            )
        geom = ElementTree.SubElement(plane, 'planegeometry')
        ElementTree.SubElement(geom, 'point',
                               {'x': '0',
                                'y': '0',
                                'z': '0'})
        ElementTree.SubElement(geom, 'normal',
                               {'x': '0',
                                'y': '0',
                                'z': '1'})
        ElementTree.SubElement(geom, 'field', type='velocity')

        line = ElementTree.SubElement(
            extr, 'propertyoutput',
            {'frequency': str(int(self.ScaleToLatticeUnits(self.MomentumDiffusionTime/4))),
             'file': 'line.xtr'}
             )
        geom = ElementTree.SubElement(line, 'linegeometry')
        ElementTree.SubElement(geom, 'point',
                               {'x': '0',
                                'y': '0',
                                'z': str(self.ScaleForOutput(-0.5 * self.DuctLength))})
        ElementTree.SubElement(geom, 'point',
                               {'x': '0',
                                'y': '0',
                                'z': str(self.ScaleForOutput(0.5 * self.DuctLength))})
        ElementTree.SubElement(geom, 'field', type='pressure')
        
        return

    @property
    def HectorCores(self):
        return int(self.TotalSites / 5000)

    @property
    def HectorSeconds(self):
        sups = self.TotalSites * self.NumberOfTimeSteps
        time = sups / (1e6 * self.HectorCores)
        return 2 * time

    
if __name__ == "__main__":
    import pdb
    import sys
    outdir = sys.argv[1]
    params = []
    for re, tau in ((1., 1.),
                    (30., 0.75),
                    (100., 0.6)):
        for size in (12, 24, 48,96, 192):
            gen = SquareDuctGeneration(re, size, tau)
            gen.Generate(outdir)

    # fields = ['ReynoldsNumber', 'Tau', 'SideLength', 'MachNumber', 'VoxelSize']
    # def Print():
    #     print ' & '.join(fields) + r' \\'

    #     for p in params:
    #         print ' & '.join(str(getattr(p, f)) for f in fields) + r' \\'
    # Print()
