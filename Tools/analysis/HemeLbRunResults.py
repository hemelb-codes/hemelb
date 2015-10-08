import os.path
import warnings
import numpy as np
import quantities as pq

from hemeTools.parsers.extraction import ExtractedProperty

from pvs.simplify import simplify
from .memo import memo_property
from .cache import cache_property

import pdb

class HemeLbRunResults(object):
    @classmethod
    def LoadFromSummary(cls, run):
        smy = os.path.join(run, 'config.smy')
        inst = super(HemeLbRunResults, cls).LoadFromSummary(smy)
        inst.RunDirectory = run
        inst.ParseRunName()
        return inst
    
    def ParseRunName(self):
        """Parse the run name into attributes.
        name = "${config}_${executable}"
        where config is dealt with by the subclass
        and executatble = "revision_velocityset_collision_wallbc"
        """
        runtypename, reynolds, dia, revision, vs, co, bc = os.path.basename(self.RunDirectory).split('_')
        
        assert runtypename == self.RunTypeName

        assert reynolds[0:2] == 'Re'
        assert (float(reynolds[2:]) - self.ReynoldsNumber)**2 < 1e-6

        assert dia[0] == 'D'
        assert int(dia[1:]) == self.ScaleToLatticeUnits(self.Diameter)

        assert int(revision, 16) > 0
        self.RevisionHash = revision

        assert vs in ('d3q15', 'd3q15i', 'd3q19', 'd3q27')            
        self.VelocitySet = vs

        assert co in ('lbgk', 'mrt')
        self.CollisionOperator = co

        assert bc in ('sbb', 'bfl', 'gzs', 'jy')
        self.BoundaryCondition = bc

        return
    
    @property
    def VsCo(self):
        return self.VelocitySet + '+' + self.CollisionOperator

    @memo_property
    def Plane(self):
        planeFileName = os.path.join(self.RunDirectory, 'results', 'Extracted', 'plane59.xtr')
        return ExtractedProperty(planeFileName)

        return ans

    @memo_property
    def Line(self):
        planeFileName = os.path.join(self.RunDirectory, 'results', 'Extracted', 'line.xtr')
        return ExtractedProperty(planeFileName)
    
    @memo_property
    def ConvergedVelocityPlane(self):
        t = self.ConvergedTimeStep
        return self.GetVelocityPlaneByTimeStep(t)
    
    def GetVelocityPlaneByTimeStep(self, t):
        return self.Plane.GetByTimeStep(t)
    
    @memo_property
    def ConvergedPressureLine(self):
        t = self.ConvergedTimeStep
        lineData = self.Line.GetByTimeStep(t)
        z = np.einsum('ij,j', lineData.position, self.Axis)
        ordering = z.argsort()
        return lineData[ordering]

    def _IterTimeVelocityPairs(self):
        def getV(i):
            return self.Plane.GetByIndex(i).velocity * (pq.metre / pq.second)
        
        times = self.Plane.times
        nT = len(times)
        queue = []
        for i, t in enumerate(times):
            queue.append((times[i], getV(i)))
            if i < self.ConvergenceTestInterval:
                continue
            
            yield queue[0], queue[self.ConvergenceTestInterval]
            queue.pop(0)
            
    @cache_property
    def ConvergedTimeStep(self):
        pdb.set_trace()
        for (t1, u1), (t2, u2) in self._IterTimeVelocityPairs():
            delta = self.CalculateDeltaU(u1, u2)
            if delta < 1e-7:
                return t1
            continue
        warnings.warn("%s not converged: delta = %e" % (self.RunDirectory, delta))
        return t1
    
    @memo_property
    @simplify
    def VelocityConvergenceError(self):
        ultimate = self.GetVelocityPlaneByTimeStep(self.Plane.times[self.ConvergedTimeStep])
        penultimate = self.GetVelocityPlaneByTimeStep(
            self.Plane.times[self.ConvergedTimeStep + self.ConvergenceTestInterval]
            )
        vDiffNormSq = np.sum((ultimate.velocity - penultimate.velocity)**2, axis=-1)
        return np.sqrt(np.mean(vDiffNormSq)) * (pq.metre / pq.second) / self.MaximumVelocity

    @simplify
    def CalculateDeltaU(self, u1, u2):
        """Calculate the maximum difference between two fields. (Eq 14)
        """
        vDiffNormSq = np.sum((u1 - u2)**2, axis=-1)
        return np.max(vDiffNormSq)**0.5 / self.MaximumVelocity

