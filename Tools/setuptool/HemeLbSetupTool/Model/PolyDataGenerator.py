import numpy as np
from vtk import vtkTransform, vtkTransformFilter
from .Vector import Vector
from .GeometryGenerator import GeometryGenerator
from .Clipper import Clipper

import Generation
import pdb

np.seterr(divide='ignore')

class PolyDataGenerator(GeometryGenerator):

    def __init__(self, profile):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        GeometryGenerator.__init__(self)
        self._profile = profile
        self.generator = Generation.PolyDataGenerator()
        self._SetCommonGeneratorProperties()
        self.generator.SetSeedPointWorking(
            profile.SeedPoint.x / profile.VoxelSize,
            profile.SeedPoint.y / profile.VoxelSize,
            profile.SeedPoint.z / profile.VoxelSize)

        # This will create the pipeline for the clipped surface
        clipper = Clipper(profile)

        # Scale by the voxel size
        trans = vtkTransform()
        scale = 1. / profile.VoxelSize
        trans.Scale(scale, scale, scale)

        transformer = vtkTransformFilter()
        transformer.SetTransform(trans)
        transformer.SetInputConnection(
            clipper.ClippedSurfaceSource.GetOutputPort())

        # Uncomment this an insert the output path to debug pipeline construction
        # write = StageWriter('/Users/rupert/working/compare/aneurysm').WriteOutput
        # i = 0
        # for alg in getpipeline(transformer):
        #     print i
        #     i += 1
        #     print alg
        #     write(alg)

        transformer.Update()
        self.ClippedSurface = transformer.GetOutput()
        self.generator.SetClippedSurface(self.ClippedSurface)
        
        originWorking, nSites = self._ComputeOriginWorking()
        self.generator.SetOriginWorking(*(float(x) for x in originWorking))
        self.generator.SetSiteCounts(*(int(x) for x in nSites))
        self.OriginMetres = Vector(originWorking * self.VoxelSizeMetres) 
        return
    
    def __getattr__(self, attr):
        """Delegate unknown attribute access to profile object.
        """
        return getattr(self._profile, attr)
    
    def _ComputeOriginWorking(self):
        """
        Here we are setting the location of our domain's origin in the input
        space and the number of sites along each axis. Sites will all have
        positions of:
               Origin + Index * VoxelSize,
        where:
               0 <= Index[i] < nSites[i]
        
        We also require that there be at least one solid site outside the fluid
        sites. For the case of axis-aligned faces which are an integer number
        of VoxelSizes apart (e.g. synthetic datasets!) this can cause numerical
        issues for the classifier if all the points that are "outside" are very
        close to the surface so we further require that these sites are a
        little further from the bounding box of the PolyData.
        """
        SurfaceBoundsWorking = self.ClippedSurface.GetBounds()
        OriginWorking = np.zeros(3, dtype=float)
        nSites = np.zeros(3, dtype=np.uint)
        
        for i in xrange(3):
            # Bounds of the vtkPolyData
            mins = SurfaceBoundsWorking[2 * i]
            maxs = SurfaceBoundsWorking[2 * i + 1]
            size = maxs - mins
            
            nSites[i] = int(size)
            # Since int() truncates, we have:
            #         0 < size/VoxelSize - nSites < 1.
            # Hence we need nSites + 1 links and therefore nSites + 2 sites
            nSites[i] += 2
            
            # The extra distance from size to the distance from x[0] to x[nSites -1]
            extra = (nSites[i] - 1) - size

            # To avoid numerical problems with the classifier, ensure that the
            # sites just outside the fluid region are at least 1% of a VoxelSize
            # away.
         
            if (extra < 0.01):
                # They weren't, so add one to the # sites and recalculate extra
                nSites[i] += 1
                extra = (nSites[i] - 1) - size
                pass

            # Now ensure this extra space is equally balanced before & after the
            # fluid region with the placement of the first site.
            OriginWorking[i] = mins - 0.5 * extra
            continue
        
        return OriginWorking, nSites
    pass


