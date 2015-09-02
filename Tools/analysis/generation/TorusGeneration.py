import os.path
import contextlib
import quantities as pq
import numpy as np
from xml.etree import ElementTree

import _fixpath

from .TorusParameters import TorusParameters, IterRunParams
from .WriteStlTorus import WriteStlTorus
from pvs.helpers import IndentXml

from HemeLbSetupTool.Model.Vector import Vector
from HemeLbSetupTool.Model.Iolets import Inlet, Outlet
from HemeLbSetupTool.Model.Profile import Profile, millimetre


class TorusGeneration(TorusParameters):
    @property
    def WarmUpSteps(self):
        return int(self.ScaleToLatticeUnits(self.SoundTime))

    @property
    def InletPressure(self):
        ans = np.zeros((3,), dtype=self.ReferencePressure.dtype) * self.ReferencePressure.units
        ans[0] = self.ReferencePressure + self.PressureDifference
        return ans
    
    @property
    def OutletPressure(self):
        ans = np.zeros((3,), dtype=self.ReferencePressure.dtype) * self.ReferencePressure.units
        ans[0] = self.ReferencePressure
        return ans

    StlFileUnit = millimetre
    
    _OutputUnitsMap = dict(**TorusParameters._OutputUnitsMap)
    _OutputUnitsMap[pq.metre.dimensionality.simplified] = pq.millimetre
                
    def MakeVector(self, array):
        assert array.shape == (3,)
        return Vector(float(x) for x in array)
    
    @property
    def Inlet(self):
        return Inlet(Name="Inlet",
                     Centre=self.MakeVector(self.ScaleForOutput(self.InletCentre)),
                     Normal=self.MakeVector(self.InletNormal),
                     Pressure=self.MakeVector(self.ScaleForOutput(self.InletPressure)),
                     Radius=1.5*self.ScaleForOutput(self.CrossSectionRadius))
    @property
    def Outlet(self):
        return Outlet(Name="Outlet",
                      Centre=self.MakeVector(self.ScaleForOutput(self.OutletCentre)),
                      Normal=self.MakeVector(self.OutletNormal),
                      Pressure=self.MakeVector(self.ScaleForOutput(self.OutletPressure)),
                      Radius=1.5*self.ScaleForOutput(self.CrossSectionRadius))

    @property
    def Profile(self):
        p = Profile()
        p.StlFile = self.StlFileName
        p.StlFileUnitId = Profile._UnitChoices.index(self.StlFileUnit)
        p.Iolets.append(self.Inlet)
        p.Iolets.append(self.Outlet)
        p.SeedPoint = self.MakeVector(self.ScaleForOutput(self.InletCentre))
        p.VoxelSize = self.ScaleForOutput(self.VoxelSize)
        p.Steps = self.NumberOfTimeSteps
        p.OutputGeometryFile = self.GmyFileName
        p.OutputXmlFile = self.OutputXmlFile
        return p

    @property
    def ShouldWriteStl(self):
        return not os.path.exists(self.StlFileName)
    
    @property
    def WriteStlParameters(self):
        """c - radius from the centre to the middle of the ring
        a - radius of the cross-section of the ring

        n_c - number of sections around the torus
        n_a - number of sections around the cross section

        filename - name of file to write
        """
        def dedim(q):
            q = q.simplified
            assert q.dimensionality == pq.dimensionless
            return float(q)
        
        return (dedim(self.RingRadius / pq.mm),
                dedim(self.CrossSectionRadius / pq.mm),

                int(dedim(2 * np.pi * (self.RingRadius + self.CrossSectionRadius) / self.VoxelSize)),
                int(dedim(2 * np.pi * self.CrossSectionRadius / self.VoxelSize)),

                self.StlFileName)

    def Generate(self, outDir):
        assert self.ScaleToLatticeUnits(self.PressureDifference) < 0.05, \
            ("Pressure difference too large (%f)" % self.ScaleToLatticeUnits(self.PressureDifference))
        
        with cd(outDir):
            if os.path.exists(self.RunName):
                print "Warning: not regenerating '%s'" % os.path.abspath(self.RunName)
                return
            
            # Create directory
            os.mkdir(self.RunName)

            # Create the STL if needed
            if self.ShouldWriteStl:
                WriteStlTorus(*self.WriteStlParameters)

            # Create the Profile
            p = self.Profile
            # Write it
            p.Save(self.ProfileFileName)

            self.WriteSummary()

            # Actually create the geometry file
            p.Generate()

            self.FixXml()

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

    @property
    def ExtractionCentreToroidal(self):
        return np.array([0.,0., self.PipeAngle / 2])
    
    @property
    def ExtractionCentre(self):
        return self.T2C.TransformPoint(self.ExtractionCentreToroidal)

    @property
    def ExtractionNormal(self):
        return self.T2C.TransformVector(self.ExtractionCentreToroidal, [0., 0., -1.]*pq.dimensionless)
    
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
        ElementTree.SubElement(geom, 'point', VecToDict(self.ExtractionCentre / pq.metre))
        ElementTree.SubElement(geom, 'normal', VecToDict(self.ExtractionNormal))
        ElementTree.SubElement(geom, 'field', type='velocity')

    def RewriteInlets(self, root):
        ioType = 'inlets'
        iolets = root.find(ioType)
        direction = {'inlets': 1.,
                     'outlets': -1.}[ioType]
        for iolet in iolets:
            pressure = iolet.find('pressure')
            iolet.remove(pressure)
            vel = ElementTree.SubElement(iolet, 'velocity',
                                         {'radius': str(self.ScaleToLatticeUnits(self.CrossSectionRadius)),
                                          'maximum': str(direction *
                                                         self.ScaleToLatticeUnits(self.MaximumVelocity))})
        return

@contextlib.contextmanager
def cd(path):
    start = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(start)
    

if __name__ == "__main__":
    import sys
    outpath = sys.argv[1]

    for params in IterRunParams():
        tGen = TorusGeneration(*params)
        print tGen.GetParametersReport()
        tGen.Generate(outpath)
        print ''
            
