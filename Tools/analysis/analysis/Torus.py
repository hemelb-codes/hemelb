import numpy as np
import quantities as pq
import os.path
import Gnuplot
import warnings
from xml.etree import ElementTree
from glob import glob

from hemeTools.parsers.extraction import ExtractedProperty

import _fixpath
from generation.TorusGeneration import TorusGeneration
from .memo import memo_property
from .cache import cache_property
from .HemeLbRunResults import HemeLbRunResults

import pdb
        
class TorusResults(HemeLbRunResults, TorusGeneration):
    ConvergenceTestInterval = 1
    
    def ParseRunName(self):
        """Parse the run name into attributes.
        name = "${config}_${executable}"
        where config is dealt with by the subclass
        and executatble = "revision_velocityset_collision_wallbc"
        """
        # orus_Dean10_Delta0.1_Voxel0.001_5ff4e11dc563_d3q15_lbgk_sbb
        runtypename, dean, delta, voxel, revision, vs, co, bc = os.path.basename(self.RunDirectory).split('_')
        
        assert runtypename == self.RunTypeName

        assert dean[0:4] == 'Dean'
        assert (float(dean[4:]) - self.DeanNumber)**2 < 1e-6

        assert delta[0:5] == 'Delta'
        assert (float(delta[5:]) - self.CurvatureRatio)**2 < 1e-6
        
        assert voxel[0:5] == 'Voxel'
        assert (float(voxel[5:])*pq.metre - self.VoxelSize)**2 < 1e-12 * pq.metre**2

        assert int(revision, 16) > 0
        self.RevisionHash = revision

        assert vs in ('d3q15', 'd3q15i', 'd3q19', 'd3q27')            
        self.VelocitySet = vs

        assert co in ('lbgk', 'mrt')
        self.CollisionOperator = co

        assert bc in ('sbb', 'bfl', 'gzs', 'jy')
        self.BoundaryCondition = bc

        return
    
    def TransformOutput(self):
        plane = self.ConvergedVelocityPlane
        
        posCart = plane.position * pq.metre
        velCart = plane.velocity * (pq.metre/pq.second)

        c2t = self.T2C.GetInverse()

        posTor = c2t.TransformPoint(posCart)
        velTor = c2t.TransformVector(posCart, velCart)
        return posTor, velTor
    
    pass

if __name__ == "__main__":
    import sys
    fab_results_dir = sys.argv[1]

    for run in glob(os.path.join(fab_results_dir, 'Torus_Dean200_Delta0.1_Voxel0.001_5ff4e11dc563_*')):
        res = TorusResults.LoadFromSummary(run)
        posTor, velTor = res.TransformOutput()
        newPosTor = posTor.copy()
        newPosTor[:,2] = 0.
        res.T2C.TransformPoint(newPosTor)
        posCart = res.T2C.TransformPoint(newPosTor)
        velCart = res.T2C.TransformVector(newPosTor, velTor)
        
        for var in ('pos', 'vel'):
            for coord in ('Cart', 'Tor'):
                np.savetxt(os.path.join(res.RunDirectory, var+coord+'.dat'), locals()[var+coord])
            
    # scale = 1e5 * pq.second/pq.metre

    # def Circle(r, **opts):
    #     return Gnuplot.Func('{r}*cos(t), {r}*sin(t)'.format(r=r), **opts)

    
    # def Axial(r, theta, dean, delta):
    #     rSq = r**2
    #     w00 = (1.0 - rSq) / 4.
    #     w01 = -3.0 * r * (1. - rSq) * np.cos(theta) / 16.
    #     w02 = (1. - rSq) * (-3. + 11. * rSq + 10. * rSq * np.cos(2.*theta)) / 128.

    #     w0 = w00 + delta * w01 + delta**2 * w02
        
    #     return dean * w0
    
    # g = Gnuplot.Gnuplot()
    # g('set parametric')
    # g.plot(Gnuplot.Data(
    #     (posCart[:,0] - res.RingRadius)/res.CrossSectionRadius,
    #     posCart[:,2]/res.CrossSectionRadius,
    #     scale*velCart[:,0],
    #     scale*velCart[:,2],
    #     **{'with': 'vectors'}), Circle(1., title=None))

    
