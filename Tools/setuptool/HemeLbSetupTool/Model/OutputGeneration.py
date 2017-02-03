#
# Copyright (C) University College London, 2007-2012, all rights reserved.
#
# This file is part of HemeLB and is CONFIDENTIAL. You may not work
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#

import numpy as np
import os.path

from vtk import vtkClipPolyData, vtkAppendPolyData, vtkPlane, vtkStripper, \
    vtkFeatureEdges, vtkPolyDataConnectivityFilter, vtkProgrammableFilter, \
    vtkTriangleFilter, vtkCleanPolyData, vtkIntArray, vtkPoints, vtkPolyData, \
    vtkCellArray, vtkTransform, vtkTransformFilter, vtkIdList, vtkPolyLine, \
    vtkXMLPolyDataWriter, vtkAlgorithm, vtkImplicitBoolean, vtkSphere, \
    vtkPolyDataNormals, vtkSTLWriter

from vmtk.vtkvmtk import vtkvmtkPolyDataBoundaryExtractor
from vmtk.vtkvmtk import vtkvmtkBoundaryReferenceSystems

#from CGAL.CGAL_Polyhedron_3 import Polyhedron_3

from .Iolets import Inlet, Outlet, Iolet
from .Vector import Vector
from .Profile import Profile, metre
from .XmlWriter import XmlWriter

import Generation
import pdb

np.seterr(divide='ignore')


def DVfromV(v):
    """Translate a Model.Vector.Vector to a Generation.DoubleVector.
    """
    return Generation.DoubleVector(v.x, v.y, v.z)


class GeometryGenerator(object):

    def __init__(self):
        self.skipNonIntersectingBlocks = False

    def _MakeIoletProxies(self):
        # Construct the Iolet structs
        nIn = 0
        nOut = 0
        ioletProxies = []
        for io in self._profile.Iolets:
            proxy = Generation.Iolet()

            proxy.Centre = DVfromV(io.Centre) / self._profile.VoxelSize
            proxy.Normal = DVfromV(io.Normal)
            proxy.Radius = io.Radius / self._profile.VoxelSize

            if isinstance(io, Inlet):
                io.Id = proxy.Id = nIn
                proxy.IsInlet = True
                nIn += 1
            elif isinstance(io, Outlet):
                io.Id = proxy.Id = nOut
                proxy.IsInlet = False
                nOut += 1
                pass
            ioletProxies.append(proxy)
            continue
        return ioletProxies

    def _SetCommonGeneratorProperties(self):
        self.generator.SetOutputGeometryFile(
            str(self._profile.OutputGeometryFile))
        # We need to keep a reference to this to make sure it's not GC'ed
        self.ioletProxies = self._MakeIoletProxies()
        self.generator.SetIolets(self.ioletProxies)
        return

    def Execute(self):
        """Forward this to the C++ implementation.
        """
        t = Timer()
        t.Start()
        self.generator.Execute(self.skipNonIntersectingBlocks)
        XmlWriter(self).Write()
        t.Stop()
        print "Setup time: %f s" % t.GetTime()
        return

    pass


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
            min = SurfaceBoundsWorking[2 * i]
            max = SurfaceBoundsWorking[2 * i + 1]
            size = max - min
            
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
            OriginWorking[i] = min - 0.5 * extra
            continue
        
        return OriginWorking, nSites
    pass


class CylinderGenerator(GeometryGenerator):

    def __init__(self, OutputGeometryFile, OutputXmlFile, VoxelSizeMetres,
                 Axis, LengthMetres, RadiusMetres,
                 InletPressure=None, OutletPressure=None):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        GeometryGenerator.__init__(self)
        self.Axis = Axis
        self.LengthMetres = LengthMetres
        self.RadiusMetres = RadiusMetres
        self.InletPressure = InletPressure
        self.OutletPressure = OutletPressure

        self._profile = Profile()
        self._profile.StlFileUnitId = Profile._UnitChoices.index(metre)
        self._profile.VoxelSize = VoxelSizeMetres
        self._profile.OutputGeometryFile = OutputGeometryFile
        self._profile.OutputXmlFile = OutputXmlFile
        self._MakeIolets()

        self.generator = Generation.CylinderGenerator()
        self._SetCommonGeneratorProperties()

        self.generator.SetCylinderLength(LengthMetres / VoxelSizeMetres)
        self.generator.SetCylinderRadius(RadiusMetres / VoxelSizeMetres)
        self.generator.SetCylinderCentre(Generation.DoubleVector(0., 0., 0.))
        self.generator.SetCylinderAxis(Generation.DoubleVector(*self.Axis))
        return

    def _MakeIolets(self):
        # Construct the Iolet structs
        inlet = Inlet()
        inlet.Centre = Vector(
            *(-0.5 * self.LengthMetres * n for n in self.Axis))
        inlet.Normal = Vector(*self.Axis)
        inlet.Radius = self.RadiusMetres
        if self.InletPressure is not None:
            inlet.Pressure = self.InletPressure
        self._profile.Iolets.append(inlet)

        outlet = Outlet()
        outlet.Centre = Vector(
            *(0.5 * self.LengthMetres * n for n in self.Axis))
        outlet.Normal = Vector(*(-n for n in self.Axis))
        outlet.Radius = self.RadiusMetres
        if self.OutletPressure is not None:
            outlet.Pressure = self.OutletPressure
        self._profile.Iolets.append(outlet)

        return

    pass

class SquareDuctGenerator(GeometryGenerator):

    def __init__(self, OutputGeometryFile, OutputXmlFile, VoxelSizeMetres,
                 OpenAxis, LengthVoxels, SideVoxels,
                 InletPressure=None, OutletPressure=None):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        GeometryGenerator.__init__(self)
        assert OpenAxis in (0, 1, 2)
        self.OpenAxis = OpenAxis
        self.LengthVoxels = LengthVoxels
        self.SideVoxels = SideVoxels
        self.Sizes = Generation.DoubleVector(SideVoxels, SideVoxels, SideVoxels)
        self.Sizes[OpenAxis] = LengthVoxels
        
        self.InletPressure = InletPressure
        self.OutletPressure = OutletPressure

        self._profile = Profile()
        self._profile.StlFileUnitId = Profile._UnitChoices.index(metre)
        self._profile.VoxelSize = VoxelSizeMetres
        self._profile.OutputGeometryFile = OutputGeometryFile
        self._profile.OutputXmlFile = OutputXmlFile
        self._MakeIolets()

        self.generator = Generation.SquareDuctGenerator()
        self._SetCommonGeneratorProperties()

        self.generator.SetOpenAxis(self.OpenAxis)
        lb = self.Sizes * -0.5
        self.generator.SetLowerBound(lb)
        ub = self.Sizes * 0.5
        self.generator.SetUpperBound(ub)
        return

    def _MakeIolets(self):
        # Construct the Iolet structs
        inlet = Inlet()
        c = [0., 0., 0.]
        c[self.OpenAxis] = -0.5 * self.LengthVoxels * self._profile.VoxelSize
        inlet.Centre = Vector(*c)
        
        n = Generation.DoubleVector()
        n[self.OpenAxis] = 1.
        inlet.Normal = Vector(n.x, n.y, n.z)
        
        inlet.Radius = self.SideVoxels * self._profile.VoxelSizeMetres
        if self.InletPressure is not None:
            inlet.Pressure = self.InletPressure
        self._profile.Iolets.append(inlet)

        outlet = Outlet()
        c = [0., 0., 0.]
        c[self.OpenAxis] = 0.5 * self.LengthVoxels * self._profile.VoxelSize
        outlet.Centre = Vector(*c)
        
        n = Generation.DoubleVector()
        n[self.OpenAxis] = -1.
        outlet.Normal = Vector(n.x, n.y, n.z)
        
        outlet.Radius = self.SideVoxels * self._profile.VoxelSizeMetres
        if self.OutletPressure is not None:
            outlet.Pressure = self.OutletPressure
        self._profile.Iolets.append(outlet)

        return

    pass

# TODO: organise this timer
import time


class Timer(object):

    def __init__(self):
        self._running = False

    def Start(self):
        assert not self._running
        self._running = True
        self._startTime = time.clock()
        return

    def Stop(self):
        assert self._running
        self._stopTime = time.clock()
        self._running = False
        return

    def GetTime(self):
        if self._running:
            return time.clock() - self._startTime
        return self._stopTime - self._startTime

    pass


class Clipper(object):
    """Clips the input STL file to the ROI and caps it.
    """

    def __init__(self, profile):
        """Create the generator for the supplied Model.Profile object.

        profile - the HemeLbSetupTool.Model.Profile object to use
        """
        # Copy the profile's _Args
        for k in profile._Args:
            setattr(self, k, getattr(profile, k))
            continue
        # Pull in the SurfaceSource too
        self.SurfaceSource = profile.StlReader

        self.ClippedSurfaceSource = self.ConstructClipPipeline()
        return

    def ConstructClipPipeline(self):
        """This constructs a VTK pipeline to clip the vtkPolyData read from
        the STL file against the Iolets. It also adds a scalar value to each
        polygon indicating which Iolet it represents, or -1 if it is just a
        wall, and computes polygon normals.

        Note that this does NOT EXECUTE the pipeline.
        """
        # seal any leaks first.
        closer = PolyDataCloser()
        closer.SetInputConnection(self.SurfaceSource.GetOutputPort())
        # Add the Iolet id -1 to all cells
        adder = IntegerAdder(Value=-1)
        adder.SetInputConnection(closer.GetOutputPort())
        # Have the name pdSource first point to the input, then loop
        # over IOlets, clipping and capping.
        pdSource = adder
        for i, iolet in enumerate(self.Iolets):
            capper = PolyDataClipCapAndLabeller(Value=i, Iolet=iolet,
                                                SeedPoint=self.SeedPoint)
            capper.SetInputConnection(pdSource.GetOutputPort())
            # Set the source of the next iteraction to the capped
            # surface producer.
            pdSource = capper
            continue

        # The following adds cell normals to the PolyData
        normer = vtkPolyDataNormals()
        normer.SetInputConnection(pdSource.GetOutputPort())
        normer.ComputeCellNormalsOn()
        normer.ComputePointNormalsOff()
        normer.SplittingOff()
        normer.ConsistencyOn()
        normer.AutoOrientNormalsOn()
        normer.NonManifoldTraversalOff()

        return normer

    pass


class IntegerAdder(vtkProgrammableFilter):
    """vtkFilter for adding an integer value to vtkPolyData's CellData. Use
    SetValue to set the value to be added, or supply the optional argument
    'Value' to the constructor.
    """

    def __init__(self, Value=-1):
        self.SetExecuteMethod(self._Execute)
        self.Value = Value
        return

    def SetValue(self, val):
        self.Value = val
        return

    def GetValue(self):
        return self.Value

    def _Execute(self, *args):
        """Run the filter.
        """
        # Get the input
        input = self.GetPolyDataInput()
        values = vtkIntArray()
        values.SetNumberOfTuples(input.GetNumberOfCells())
        values.FillComponent(0, self.Value)

        # Get the out PD
        out = self.GetPolyDataOutput()
        out.CopyStructure(input)
        out.GetCellData().SetScalars(values)

        return

    pass


class PolyDataClipCapAndLabeller(vtkProgrammableFilter):
    """vtkFilter for clipping and capping a vtkPolyData surface, and labeling
    the cap with an integer cell data value.
    """
    def __init__(self, Value=None, Iolet=None, SeedPoint=None):
        self.SetExecuteMethod(self._Execute)
        self.Value = Value
        self.Iolet = Iolet
        self.SeedPoint = (SeedPoint.x, SeedPoint.y, SeedPoint.z)
        return

    def SetValue(self, val):
        self.Value = val
        return

    def GetValue(self):
        return self.Value

    def SetIolet(self, val):
        self.Iolet = val
        return

    def GetIolet(self):
        return self.Iolet

    def _Clip(self, pd):
        # The plane implicit function will be >0 for all the points in the positive side
        # of the plane (i.e. x s.t. n.(x-o)>0, where n is the plane normal and o is the
        # plane origin).
        plane = vtkPlane()
        plane.SetOrigin(
            self.Iolet.Centre.x, self.Iolet.Centre.y, self.Iolet.Centre.z)
        plane.SetNormal(
            self.Iolet.Normal.x, self.Iolet.Normal.y, self.Iolet.Normal.z)

        # The sphere implicit function will be >0 for all the points outside
        # the sphere.
        sphere = vtkSphere()
        sphere.SetCenter(
            self.Iolet.Centre.x, self.Iolet.Centre.y, self.Iolet.Centre.z)
        sphere.SetRadius(self.Iolet.Radius)

        # The VTK_INTERSECTION operator takes the maximum value of all the registered
        # implicit functions. This will result in the function evaluating to >0 for all
        # the points outside the sphere plus those inside the sphere in the positive
        # side of the plane.
        clippingFunction = vtkImplicitBoolean()
        clippingFunction.AddFunction(plane)
        clippingFunction.AddFunction(sphere)
        clippingFunction.SetOperationTypeToIntersection()

        clipper = vtkClipPolyData()
        clipper.SetInput(pd)
        clipper.SetClipFunction(clippingFunction)

        # Filter to get part closest to seed point
        connectedRegionGetter = vtkPolyDataConnectivityFilter()
        connectedRegionGetter.SetExtractionModeToClosestPointRegion()
        connectedRegionGetter.SetClosestPoint(*self.SeedPoint)
        connectedRegionGetter.SetInputConnection(clipper.GetOutputPort())
        connectedRegionGetter.Update()
        return connectedRegionGetter.GetOutput()

    def _AddValue(self, pd):
        adder = IntegerAdder(Value=self.Value)
        adder.SetInput(pd)
        adder.Update()
        return adder.GetOutput()

    def _Execute(self, *args):
        assert isinstance(self.Value, int)
        assert isinstance(self.Iolet, Iolet)
        input = self.GetPolyDataInput()
        clipped = self._Clip(input)

        clipped.BuildLinks(0)

        newPoints = vtkPoints()
        newPoints.DeepCopy(clipped.GetPoints())

        newPolys = vtkCellArray()
        newPolys.DeepCopy(clipped.GetPolys())

        newData = vtkIntArray()
        newData.DeepCopy(clipped.GetCellData().GetScalars())

        boundaryExtractor = vtkvmtkPolyDataBoundaryExtractor()
        boundaryExtractor.SetInput(clipped)
        boundaryExtractor.Update()
        boundaries = boundaryExtractor.GetOutput()
        boundariesPointIdMap = boundaries.GetPointData().GetScalars()
        for i in xrange(boundaries.GetNumberOfCells()):
            boundary = vtkPolyLine.SafeDownCast(boundaries.GetCell(i))

            barycentre = [0., 0., 0.]
            vtkvmtkBoundaryReferenceSystems.ComputeBoundaryBarycenter(
                boundary.GetPoints(),
                barycentre)

            barycentreId = newPoints.InsertNextPoint(barycentre)

            numberOfBoundaryPoints = boundary.GetNumberOfPoints()
            trianglePoints = vtkIdList()
            trianglePoints.SetNumberOfIds(3)

            for j in xrange(numberOfBoundaryPoints):
                trianglePoints.SetId(0,
                                     boundariesPointIdMap.GetValue(boundary.GetPointId(j)))
                trianglePoints.SetId(1, barycentreId)
                trianglePoints.SetId(2,
                                     boundariesPointIdMap.GetValue(boundary.GetPointId((j + 1) % numberOfBoundaryPoints))
                                     )
                newPolys.InsertNextCell(trianglePoints)
                newData.InsertNextValue(self.Value)
                continue

            continue

        output = self.GetPolyDataOutput()
        output.SetPoints(newPoints)
        output.SetPolys(newPolys)
        output.GetCellData().SetScalars(newData)

        return

    pass


class PolyDataCloser(vtkProgrammableFilter):
    """Given an input vtkPolyData object, close any holes in the surface.
    """

    def __init__(self):
        self.SetExecuteMethod(self._Execute)

    def _Execute(self, *args):
        input = self.GetPolyDataInput()

        # Gets any edges of the mesh
        edger = vtkFeatureEdges()
        edger.BoundaryEdgesOn()
        edger.FeatureEdgesOff()
        edger.NonManifoldEdgesOff()
        edger.ManifoldEdgesOff()
        edger.SetInput(input)

        # Converts the edges to a polyline
        stripper = vtkStripper()
        stripper.SetInputConnection(edger.GetOutputPort())
        stripper.Update()

        # Change the polylines into polygons
        boundaryPoly = vtkPolyData()
        boundaryPoly.SetPoints(stripper.GetOutput().GetPoints())
        boundaryPoly.SetPolys(stripper.GetOutput().GetLines())

        # Triangulate
        tri = vtkTriangleFilter()
        tri.SetInput(boundaryPoly)
        tri.Update()

        # Join to the input
        merger = vtkAppendPolyData()
        merger.AddInput(input)
        merger.AddInput(tri.GetOutput())

        # Clean up by merging duplicate points
        cleaner = vtkCleanPolyData()
        cleaner.SetInputConnection(merger.GetOutputPort())
        cleaner.Update()

        # Set output
        output = self.GetPolyDataOutput()
        output.ShallowCopy(cleaner.GetOutput())
        return
    pass



# Debugging helpers


def getpreviousstage(algo):
    """Given a vtkAlgorithm, get the previous algorithm in its pipeline.
    """
    return algo.GetInputConnection(0, 0).GetProducer()


def getpipeline(last):
    """Return a list of all the algorithms in a pipeline, given the last one.
    """
    pipe = [last]
    while True:
        try:
            prev = getpreviousstage(pipe[0])
            pipe.insert(0, prev)
        except:
            break
        continue

    return pipe


class StageWriter(object):
    """For debugging, will easily let one write polydata.
    """

    def __init__(self, dir):
        if os.path.exists(dir):
            assert os.path.isdir(dir)
        else:
            os.mkdir(dir)
            pass
        self.dir = dir
        self.i = 0
        return

    def WriteOutput(self, stage, name=None):
        writer = vtkXMLPolyDataWriter()

        if isinstance(stage, vtkAlgorithm):
            # To ensure everything executes in order
            stage.Update()
            writer.SetInputConnection(stage.GetOutputPort())
        elif isinstance(stage, vtkPolyData):
            writer.SetInput(stage)
        else:
            raise ValueError(
                'Cannot cope with instances of "%s"' % type(stage))

        if name is None:
            fnPattern = '%02d.vtp'
        else:
            fnPattern = '%02d-' + name + '.vtp'
        filename = os.path.join(self.dir, fnPattern % self.i)
        self.i += 1
        writer.SetFileName(filename)
        writer.Write()
        return
    pass
