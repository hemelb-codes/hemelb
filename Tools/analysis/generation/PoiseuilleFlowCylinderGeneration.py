#!/usr/bin/env python
import os.path
import numpy as np
import quantities as pq
from xml.etree import ElementTree
import shutil
import glob
import time

import _fixpath
from .CylinderGeneration import CylinderGeneration

import pdb

class PoiseuilleFlowCylinderGeneration(CylinderGeneration):
    def __init__(self, reynolds_number, diameter_voxels, tau):
        self.ReynoldsNumber = reynolds_number * pq.dimensionless
        self._DiameterVoxels = diameter_voxels
        self._Tau = tau
        return
    
    @property
    def SimulationTime(self):
        return 4 * self.MomentumDiffusionTime
    
    @property
    def NumberOfTimeSteps(self):
        nSteps = self.SimulationTime / self.TimeStep
        return int(nSteps.simplified)

    @property
    def WarmUpSteps(self):
        return int(self.ScaleToLatticeUnits(self.SoundTime))
    
    @property
    def TotalSites(self):
        return self.ScaleToLatticeUnits(np.pi * self.Radius**2 * self.Length)

    @property
    def Axis(self):
        approx = np.array((-0.299, 0.382, 0.874))
        return pq.dimensionless * approx/np.sqrt(np.dot(approx,approx))
    
    @property
    def RunName(self):
        return 'Cylinder_Re{}_D{}'.format(self.ScaleForOutput(self.ReynoldsNumber),
                                          int(self.ScaleToLatticeUnits(self.Diameter)))
                
    def WriteSummary(self):
        import yaml
        data = {}
        data['type'] = 'cylinder'
        
        end1 = self.ScaleForOutput(0.5 * self.Length * self.Axis)
        end2 = -1 * end1
        
        data['ends'] = [self.ScaleForOutput(self.InletPosition).tolist(),
                        self.ScaleForOutput(self.OutletPosition).tolist()]
        data['radius'] = self.ScaleForOutput(self.Radius)
        data['flow_type'] = 'poiseuille'

        for attr in ['ReynoldsNumber', '_DiameterVoxels', 'Tau']:
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

        assert data['type'] == 'cylinder'
        
        args = []
        for attr in ['ReynoldsNumber', '_DiameterVoxels', 'Tau']:
            args.append(data[attr])
        return cls(*args)
        
    def RewriteInlets(self, root):
        ioType = 'inlets'
        iolets = root.find(ioType)
        direction = {'inlets': 1.,
                     'outlets': -1.}[ioType]
        for iolet in iolets:
            pressure = iolet.find('pressure')
            iolet.remove(pressure)
            vel = ElementTree.SubElement(iolet, 'velocity',
                                         {'radius': str(self.ScaleToLatticeUnits(self.Radius)),
                                          'maximum': str(direction *
                                                         self.ScaleToLatticeUnits(self.MaximumVelocity))})
        return

    @property
    def ExtractionInterval(self):
        return int(self.ScaleToLatticeUnits(self.MomentumDiffusionTime/4))
    
    pass

if __name__ == "__main__":
    import sys
    outdir = sys.argv[1]
    for re, tau in ((1., 1.),
                    (30., 0.75),
                    (100., 0.6),
                    (300., 0.55)):
        for size in (12, 24, 48, 96, 192):
            gen = PoiseuilleFlowCylinderGeneration(re, size, tau)
            gen.Generate(outdir)
