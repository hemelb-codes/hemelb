import numpy as np
from vtk import vtkTransform, vtkTransformFilter
from vtk.util import numpy_support

from .Vector import Vector
#from .GeometryGenerator import GeometryGenerator
from .Clipper import Clipper

from SurfaceVoxeliser import SurfaceVoxeliser

import pdb

# class PolyDataGenerator(GeometryGenerator):
class PolyDataGenerator(object):

    def __init__(self, profile):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        #GeometryGenerator.__init__(self)
        self._profile = profile
        #self.generator = Generation.PolyDataGenerator()
        #self._SetCommonGeneratorProperties()
        #self.generator.SetSeedPointWorking(
        #    profile.SeedPoint.x / profile.VoxelSize,
        #    profile.SeedPoint.y / profile.VoxelSize,
        #    profile.SeedPoint.z / profile.VoxelSize)
        
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

        transformer.Update()
        self.ClippedSurface = transformer.GetOutput()
        #self.generator.SetClippedSurface(self.ClippedSurface)
        
        self._ComputeOriginWorking()
        #self.generator.SetOriginWorking(*(float(x) for x in originWorking))
        #self.generator.SetSiteCounts(*(int(x) for x in nSites))
        self.OriginMetres = Vector(self.OriginWorking * self.VoxelSizeMetres) 
        return
    
    def Execute(self):
        # Just override the base's exec method for now...
        
        # Get the surface data into numpy arrays, assuming they're all triangles.
        surf = self.ClippedSurface
        tridata = numpy_support.vtk_to_numpy(surf.GetPolys().GetData())
        tridata.shape = (surf.GetNumberOfPolys(), 4)
        assert np.all(tridata[:,0]==3), "Non triangle in surface!"
        triangles = tridata[:, 1:]
        points = numpy_support.vtk_to_numpy(surf.GetPoints().GetData())
        normals = self.ClippedSurface.GetCellData().GetNormals()
        normals = numpy_support.vtk_to_numpy(normals)

        # This guy will create an octree with nodes that represent the polydata
        # surface. They are tagged with those triangles that could intersect 
        # their 26-neighbourhood links        
        voxer = SurfaceVoxeliser(points, triangles, normals, self.TreeLevels, self.OriginWorking)
        voxer.Execute()
        
        # We aren't done yet, but that will do for today!
        voxer.Tree.Write(self.OutputGeometryFile)
        
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
        nSites = np.zeros(3, dtype=int)
        
        for i in xrange(3):
            # Bounds of the vtkPolyData
            bmin = SurfaceBoundsWorking[2 * i]
            bmax = SurfaceBoundsWorking[2 * i + 1]
            size = bmax - bmin
            
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
            
        # Get the biggest of the sides
        nSites = int(nSites.max())
        # Now round this up to the next power of two
        nBits = nSites.bit_length()
        if nSites > 2**(nBits -1):
            nSites = 2**nBits
            pass
        
        for i in xrange(3):
            # Now ensure this extra space is equally balanced before & after the
            # fluid region with the placement of the first site.
            bmin = SurfaceBoundsWorking[2 * i]
            bmax = SurfaceBoundsWorking[2 * i + 1]
            size = bmax - bmin
            extra = (nSites - 1) - size
            OriginWorking[i] = bmin - 0.5 * extra
            continue
        
        self.OriginWorking = OriginWorking
        self.CubeSize = nSites
        self.TreeLevels = nSites.bit_length() - 1
        return
    pass


