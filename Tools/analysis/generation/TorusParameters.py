import os.path
import numpy as np
import quantities as pq

import _fixpath
from pvs.simplify import simplify
from pvs.coordinates import toroidal
from .HemeLbParameters import HemeLbParameters
from analysis.memo import memo_property

class TorusParameters(HemeLbParameters):
    """Encapsulate all parameters about the geometry of the torus, how
    much of it to use for simulation, the discretisation of the domain
    to a HemeLB gmy and run options.
   
    """
    RunTypeName = 'Torus'
    
    def __init__(self, curvatureRatio, deanNumber, crossSectionRadius_m, voxelSize_m):
        self.CurvatureRatio = curvatureRatio * pq.dimensionless
        self.DeanNumber = deanNumber * pq.dimensionless
        self._VoxelSize = voxelSize_m * pq.metre
        self.CrossSectionRadius = crossSectionRadius_m * pq.metre
        
        # Fixed by vtkParamatricTorus
        self.Origin = [0., 0., 0.] * pq.metre
        # Chosen as relevant to blood flow
        
        return

    @classmethod
    def LoadFromSummary(cls, summary):
        import yaml
        with file(summary) as f:
            data = yaml.load(f)
            
        assert data['type'] == 'torus'
        
        args = []
        for attr in ['CurvatureRatio', 'DeanNumber', 'CrossSectionRadius', 'VoxelSize']:
            args.append(data[attr])
        return cls(*args)

    def WriteSummary(self):
        import yaml
        
        data = {}
        data['type'] = 'torus'
        with self.SIunits():
            for attr in ['CurvatureRatio', 'DeanNumber', 'CrossSectionRadius', 'VoxelSize']:
                val = getattr(self, attr)
                if isinstance(val, pq.Quantity):
                    val = self.ScaleForOutput(val)
                data[attr] = val
        
        with file(self.OutputSummaryFile, 'w') as f:
            yaml.dump(data, f)
        return
    
    @property
    def Tau(self):
        return 0.55

    @property
    @simplify
    def TimeStep(self):
        kin_visc_lat = (self.Tau - 0.5) / 3.
        return kin_visc_lat * self.VoxelSize**2 / self.KinematicViscosity

    @property
    def VoxelSize(self):
        return self._VoxelSize
    
    @property
    def CrossSectionDiameter(self):
        return 2 * self.CrossSectionRadius
    
    @property
    @simplify
    def ReynoldsNumber(self):
        """Given the definition of the Dean number (D = 4 Re Sqrt[2
        delta]), compute the Reynolds number based on the maximum
        axial velocity of the flow that would be driven by the same
        steady pressure gradient in a straight pipe with the same
        radius as the cross section.
        """
        return self.DeanNumber * (32. * self.CurvatureRatio)**-0.5
    
    @property
    def RingRadius(self):
        return self.CrossSectionRadius / self.CurvatureRatio
    
    @property
    @simplify
    def MomentumDiffusionTime(self):
        """Time for momentum to diffuse from the surface to the centre
        of the pipe.
        """
        return self.CrossSectionRadius**2 / self.KinematicViscosity
    
    @property
    @simplify
    def MaximumVelocity(self):
        return self.ReynoldsNumber * self.KinematicViscosity / self.CrossSectionDiameter
    
    @property
    @simplify
    def MachNumber(self):
        return self.MaximumVelocity / self.SpeedOfSound
    
    @property
    @simplify
    def EntranceLength(self):
        return 0.06 * self.ReynoldsNumber * self.CrossSectionDiameter
        
    @property
    @simplify
    def PipeLength(self):
        # greatest of (4 entrance lengths and half the torus)
        ans = max(4. * self.EntranceLength, 0.9*np.pi*self.RingRadius)
        assert ans < 2 * np.pi * self.RingRadius, \
            "Pipe length must be less than the torus circumference"
        return ans
    
    @property
    def PipeAngle(self):
        return self.PipeLength / self.RingRadius
    
    @property
    @simplify
    def PressureDifference(self):
        return 4. * self.MaximumVelocity * self.DynamicViscosity * self.PipeLength / \
            self.CrossSectionRadius**2

    @property
    @simplify
    def DensityDifference(self):
        return self.PressureDifference / self.SpeedOfSound**2
    
    @property
    @simplify
    def TimeStep(self):
        """Choose the timestep such that the maximum velocity == 0.05 c_s
        """
        return 0.05 * self.VoxelSize / \
            (np.sqrt(3.) * self.MaximumVelocity)
    
    @property
    @simplify
    def SoundTime(self):
        return 2. * np.sqrt(3) * self.PipeLength * self.TimeStep / self.VoxelSize
    
    @property
    @simplify
    def SimulationTime(self):
        return self.NumberOfTimeSteps * self.TimeStep

    @property
    @simplify
    def ExtractionInterval(self):
        return int(self.ScaleToLatticeUnits(max(self.SoundTime, self.MomentumDiffusionTime)))
    
    @property
    def NumberOfTimeSteps(self):
        return 4 * self.ExtractionInterval

    @property
    def TotalSites(self):
        return self.ScaleToLatticeUnits(np.pi * self.CrossSectionRadius**2 * self.PipeLength)

    @memo_property
    def T2C(self):
        return toroidal.ToroidalToCartesianTransformer(self.CrossSectionRadius, self.RingRadius)

    @property
    def InletCentreToroidal(self):
        return np.array([0.,0.,0.])
    
    @property
    def OutletCentreToroidal(self):
        return np.array([0.,0.,self.PipeAngle])
    
    @property
    def InletCentre(self):
        return self.T2C.TransformPoint(self.InletCentreToroidal)
    
    @property
    def OutletCentre(self):
        return self.T2C.TransformPoint(self.OutletCentreToroidal)

    @property
    def InletNormal(self):
        return self.T2C.TransformVector(self.InletCentreToroidal, [0., 0., 1.])
    
    @property
    def OutletNormal(self):
        return self.T2C.TransformVector(self.OutletCentreToroidal, [0., 0., -1.])
    
    @property
    def StlFileName(self):
        return 'Torus_Radius{:g}_Delta{:g}.stl'.format(float(self.CrossSectionRadius),
                                                       float(self.CurvatureRatio))
    
    @property
    def RunName(self):
        return 'Torus_Dean{:g}_Delta{:g}_Voxel{:g}'.format(float(self.DeanNumber),
                                                           float(self.CurvatureRatio),
                                                           float(self.VoxelSize))
    @property
    def ProfileFileName(self):
        return os.path.join(self.RunName, 'config.pro')
    
    @property
    def OutputSummaryFile(self):
        return os.path.join(self.RunName, 'config.smy')
    
    @property
    def GmyFileName(self):
        return os.path.join(self.RunName, 'config.gmy')

    @property
    def OutputXmlFile(self):
        return os.path.join(self.RunName, 'config.xml')
    
    pass

    _ReportedParameters = [
        'RunName', 'CrossSectionRadius', 'CurvatureRatio', 'RingRadius',
        'PipeLength', 'EntranceLength', 'ReynoldsNumber', 'DeanNumber',
        'MachNumber', 'MomentumDiffusionTime', 'SoundTime',
        'TotalSites', 'NumberOfTimeSteps', 'Tau', 'DensityDifference']
    
    def GetParametersReport(self):
        return '\n'.join(self._FormatParamReport(attr)
                         for attr in self._ReportedParameters)
    
    def _FormatParamReport(self, attr):
        val = getattr(self, attr)
        if isinstance(val, pq.Quantity):
            val = self.ScaleToLatticeUnits(val)
        
        return '{0} = {1}'.format(attr, val)

def IterRunParams():
    curvatureRatio = 0.1
    crossSectionRadius_m = 24e-3
    voxelSize_m = 1e-3
    for deanNumber in (200.,):
        yield (curvatureRatio, deanNumber, crossSectionRadius_m, voxelSize_m)
    
if __name__ == "__main__":
    
    reports = []
    for pars in IterRunParams():
        tp = TorusParameters(*pars)
        reports.append(tp.GetParametersReport())
    
    print '\n\n'.join(reports)
