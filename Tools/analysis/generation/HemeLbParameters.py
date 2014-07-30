import abc
import quantities as pq
import numpy as np
from contextlib import contextmanager

import _fixpath
from pvs.simplify import simplify

class HemeLbParameters(object):
    __metaclass__ = abc.ABCMeta
    
    mmHg = 133.3223874 * pq.pascal
    ReferencePressure = 0.0 * mmHg
    
    Density = 1e3 * pq.kilogram / pq.metre**3
    DynamicViscosity = 0.004 * pq.pascal*pq.second
    KinematicViscosity = DynamicViscosity / Density
    
    @property
    @simplify
    def SpeedOfSound(self):
        return self.VoxelSize / (np.sqrt(3.) * self.TimeStep)

    @abc.abstractproperty
    def VoxelSize(self):
        pass
    
    @abc.abstractproperty
    def TimeStep(self):
        pass
    
    @property
    @simplify
    def UnitMass(self):
        """Choose this such that the density is one in lattice units.
        """
        return self.Density * self.VoxelSize**3

    @property
    def Tau(self):
        return self.ScaleToLatticeUnits(self.KinematicViscosity / self.SpeedOfSound**2) + 0.5
    
    @property
    def _LatticeUnitsMap(self):
        return {
            pq.second: self.TimeStep,
            pq.metre: self.VoxelSize,
            pq.kilogram: self.UnitMass
            }
    
    def ScaleToLatticeUnits(self, q):
        # Check it's a quantity
        assert isinstance(q, pq.Quantity)
        dims = q.dimensionality.simplified

        scaleFactor = pq.dimensionless
        for dim, power in dims.iteritems():
            scaleFactor = scaleFactor * self._LatticeUnitsMap[dim]**power
                    
        # Scale
        scaled = (q / scaleFactor).simplified
        assert scaled.units == pq.dimensionless
        
        # Now forget about quantities
        scaled = scaled.view(np.ndarray)

        # If it's a scalar, return a float
        if q.shape == ():
            return float(scaled)
        # otherwise, the np array
        return scaled
    
    _OutputUnitsMap = {}
    
    @contextmanager
    def SIunits(self):
        self._OutputUnitsMap = {}
        yield
        del self._OutputUnitsMap
        

    def ScaleForOutput(self, q):
        # Check it's a quantity
        assert isinstance(q, pq.Quantity)
        try:
            # Have we specified output units for quantities with these dimensions?
            scaleFactor = self._OutputUnitsMap[q.dimensionality.simplified]
        except KeyError:
            # No, so use SI base units instead
            scaleFactor = q.units.simplified.units
            pass
        
        # Scale
        scaled = (q / scaleFactor).simplified
        assert scaled.units == pq.dimensionless
        
        # Now forget about quantities
        scaled = scaled.view(np.ndarray)

        # If it's a scalar, return a float
        if q.shape == ():
            return float(scaled)
        # otherwise, the np array
        return scaled
    
    pass

