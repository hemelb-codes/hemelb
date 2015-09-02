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
from CylinderGeneration import CylinderGeneration

class WomersleyFlowCylinderGeneration(CylinderGeneration):
    def __init__(self, reynolds_number, target_alpha, diameter_voxels, tau):
        self.ReynoldsNumber = reynolds_number * pq.dimensionless
        self._DiameterVoxels = diameter_voxels
        self._Tau = tau
        
        # We need to nudge this such that QUARTER of the period is an exact number of timesteps
        self.WomersleyNumber = target_alpha * pq.dimensionless
        # Using the target Wo, compute the period in lattice units and round to an int
        quarterPeriodSteps = int(np.round(self.ScaleToLatticeUnits(self.Period) / 4.))
        periodSteps = 4 * quarterPeriodSteps
        # Convert it back to a physical quantity
        period = periodSteps * self.TimeStep
        # Compute the new Wo (bearing in mind that pq doesn't handle sqrt properly)
        woSq = ((2 * np.pi * self.Radius**2) / (self.KinematicViscosity * period)).simplified
        assert woSq.units == pq.dimensionless
        # Set the Wo we'll actually use
        self.WomersleyNumber = np.sqrt(woSq)
        return
    
    @property
    def SimulationTime(self):
        return 4 * self.MomentumDiffusionTime
    
    @property
    def ExtractionInterval(self):
        return int(np.round(self.ScaleToLatticeUnits(self.Period/4)))

    @property
    def WarmUpSteps(self):
        return 0
    
    @property
    @simplify
    def Period(self):
        return (2 * np.pi * self.Radius**2) / (self.KinematicViscosity * self.WomersleyNumber**2)
    
    @property
    def RunName(self):
        return 'Womersley_Re{}_D{}'.format(self.ScaleForOutput(self.ReynoldsNumber),
                                           int(np.round(self.ScaleToLatticeUnits(self.Diameter))))
    
    def WriteSummary(self):
        import yaml
        data = {}
        data['type'] = 'cylinder'
        
        end1 = self.ScaleForOutput(0.5 * self.Length * self.Axis)
        end2 = -1 * end1
        
        data['ends'] = [self.ScaleForOutput(self.InletPosition).tolist(),
                        self.ScaleForOutput(self.OutletPosition).tolist()]
        data['radius'] = self.ScaleForOutput(self.Radius)
        data['flow_type'] = 'womersley'

        for attr in ['ReynoldsNumber', 'WomersleyNumber', '_DiameterVoxels', 'Tau']:
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
        assert data['flow_type'] == 'womersley'
        
        args = []
        for attr in ['ReynoldsNumber', 'WomersleyNumber', '_DiameterVoxels', 'Tau']:
            args.append(data[attr])
        return cls(*args)
        
    def RewriteInlets(self, root):
        ioType = 'inlets'
        iolets = root.find(ioType)
        direction = {'inlets': 1.,
                     'outlets': -1.}[ioType]
        for iolet in iolets:
                # <inlet>
                #   <womersley_velocity> 
                #     <pressure_gradient_amplitude value="2.5" units="lattice"/>
                #     <period value="0.5" units="lattice"/>
                #     <womersley_number value="2" units="dimensionless"/>
                #     <radius value="0.10" units="lattice"/>
                #     <normal x="0.0" y="0.0" z="1.0" />
                #     <position x="0.0" y="0.0" z="-0.05" />
                #   </womersley_velocity>
                # </inlet>
            pressure = iolet.find('pressure')
            iolet.remove(pressure)
            normal = iolet.find('normal')
            iolet.remove(normal)
            position = iolet.find('position')
            iolet.remove(position)
            
            vel = ElementTree.SubElement(iolet, 'womersley_velocity')
            ElementTree.SubElement(
                vel, 'pressure_gradient_amplitude',
                {'value': str(self.ScaleToLatticeUnits(self.PressureGradient)),
                 'units': 'lattice'}
                )
            ElementTree.SubElement(
                vel, 'period',
                {'value': str(int(np.round(self.ScaleToLatticeUnits(self.Period)))),
                 'units': 'lattice'}
                )
            ElementTree.SubElement(
                vel, 'womersley_number',
                {'value': str(self.ScaleToLatticeUnits(self.WomersleyNumber)),
                 'units': 'dimensionless'}
                )
            ElementTree.SubElement(
                vel, 'radius',
                {'value': str(self.ScaleToLatticeUnits(self.Radius)),
                 'units': 'lattice'}
                 )
            ElementTree.SubElement(vel, 'normal', normal.attrib)
            ElementTree.SubElement(vel, 'position', position.attrib)
        return
    pass

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument('--table', action="store_true")
    p.add_argument('outdir', nargs='?')
    args = p.parse_args()
    
    for re, wo, tau in (
            (30., 4., 0.620),
            (100., 8., 0.565),
            (300., 12., 0.531)):
        
        for size in (48, ):
            gen = WomersleyFlowCylinderGeneration(re, wo, size, tau)
            def scale(gen, attr):
                q = getattr(gen, attr)
                if isinstance(q, pq.Quantity):
                    return gen.ScaleToLatticeUnits(q)
                else:
                    return q
                
            if args.table:
                for attr in ['ReynoldsNumber', 'WomersleyNumber', 'MachNumber', 'Tau', 'PressureDifference', 'Period', 'MomentumDiffusionTime', 'SoundTime']:
                    print attr, ' ', scale(gen, attr)

                print ''
            else:
                gen.Generate(args.outdir)
            
