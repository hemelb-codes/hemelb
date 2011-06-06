import numpy as np
import os.path
from xml.etree.ElementTree import Element, SubElement, ElementTree

from vtk import vtkClipPolyData, vtkAppendPolyData, vtkPlane, \
    vtkStripper, vtkFeatureEdges, vtkPolyDataConnectivityFilter, \
    vtkPolyDataNormals, vtkProgrammableFilter, vtkDelaunay2D, \
    vtkTriangleFilter, vtkCleanPolyData, vtkIntArray, vtkOctreePointLocator, \
    vtkImplicitFunction, vtkPoints, vtkPolyData, vtkCellArray
from vtk import vtkXMLPolyDataWriter
from .Iolets import Inlet, Outlet
import Generation
import pdb

np.seterr(divide='ignore')
def DVfromV(v):
    """Translate a Model.Vector.Vector to a Generation.DoubleVector.
    """
    dv = Generation.DoubleVector()
    dv.x = v.x
    dv.y = v.y
    dv.z = v.z
    return dv

class ConfigGenerator(object):
    def __init__(self, profile):
        """Clip the STL and set attributes on the SWIG-proxied C++ ConfigGenerator object.
        """
#        pdb.set_trace()
        self.profile = profile

        self.generator = Generation.ConfigGenerator()
        self.generator.SetVoxelSize(profile.VoxelSize)
        self.generator.SetOutputConfigFile(str(profile.OutputConfigFile))
        self.generator.SetStressType(profile.StressType)

        # Construct the Iolet structs
        nIn = 0
        nOut = 0
        ioletProxies = []
        for io in profile.Iolets:
            proxy = Generation.Iolet()

            proxy.Centre = DVfromV(io.Centre)
            proxy.Normal = DVfromV(io.Normal)
            proxy.Radius = io.Radius

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
        # We need to keep a reference to this to make sure it's not GC'ed
        self.ioletProxies = ioletProxies
        self.generator.SetIolets(ioletProxies)

        self.generator.SetSeedPoint(profile.SeedPoint.x,
                                    profile.SeedPoint.y,
                                    profile.SeedPoint.z)

        # This will create the pipeline for the clipped surface and excute it
        clipper = Clipper(profile)
        clipper.ClippedSurfaceSource.Update()
        self.ClippedSurface = clipper.ClippedSurfaceSource.GetOutput()
        self.generator.SetClippedSurface(self.ClippedSurface)

        # Debugging stuff: output the clipped geometry
        writer = vtkXMLPolyDataWriter()
        stlPath, stl = os.path.split(profile.StlFile)
        stlBase, void = os.path.splitext(stl)
        writer.SetFileName(os.path.join(stlPath, 'clipped_' + stlBase + '.vtp'))
        writer.SetInput(self.ClippedSurface)
        writer.Write()

        # Create the locator
#        clipper.ClippedSurfaceSource.Update()
#        clippedSurface = clipper.ClippedSurfaceSource.GetOutput()

#        self.Locator = vtkOBBTree()
#        self.Locator.SetDataSet(clippedSurface)
#        self.Locator.BuildLocator()
        # Again, we need to make sure this is kept alive and not GC'ed
#        self.generator.SetLocator(self.Locator)
        return

    def Execute(self):
        """Forward this to the C++ implementation.
        """
        self.generator.Execute()
        XmlWriter(self.profile).Write()
        return
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
        self.SurfaceSource = profile.SurfaceSource

        self.ClippedSurfaceSource = self.ConstructClipPipeline()
        return

    def ConstructClipPipeline(self):
        """This constructs a VTK pipeline to clip the vtkPolyData read
        from the STL file against the Iolets. It also adds a scalar
        value to each polygon indicating which Iolet it represents, or
        -1 if it is just a wall, and computes polygon normals.
        
        Note that this does NOT EXECUTE the pipeline.
        """
        # TODO: should probably seal any leaks first.
        
        # Add the Iolet id -1 to all cells
        adder = IntegerAdder(Value= -1)
        adder.SetInputConnection(self.SurfaceSource.GetOutputPort())
        # Have the name pdSource first point to the input, then loop
        # over IOlets, clipping and capping.
        pdSource = adder
        for i, iolet in enumerate(self.Iolets):
            # TODO: switch from this simple ImplicitPlane (i.e.
            # infinite plane) clipping to excising a thin cuboid
            # from the Iolet and keeping the piece next to the
            # SeedPoint.

            # Create a plane to clip against.
            plane = vtkPlane()
            plane.SetOrigin(iolet.Centre.x, iolet.Centre.y, iolet.Centre.z)
            plane.SetNormal(iolet.Normal.x, iolet.Normal.y, iolet.Normal.z)

#            # Clips to give the surface we want
#            clipper = vtkClipPolyData()
#            clipper.SetClipFunction(plane)
#            clipper.SetInputConnection(pdSource.GetOutputPort())

            capper = PolyDataClipCapAndLabeller(Value=i, ClipFunction=plane)
            capper.SetInputConnection(pdSource.GetOutputPort())

            # Set the source of the next iteraction to the capped
            # surface producer.
            pdSource = capper
            continue

        # Adds cell normals to the PolyData
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
        if name is None:
            fnPattern = '%02d.vtp'
        else:
            fnPattern = name + '%02d.vtp'

        filename = os.path.join(self.dir, fnPattern % self.i)
        writer.SetInputConnection(stage.GetOutputPort())
        writer.SetFileName(filename)
        writer.Write()
        self.i += 1
        return
    pass
#ws = StageWriter('/Users/rupe/tmp/setuptooltest/stages').WriteOutput

class IntegerAdder(vtkProgrammableFilter):
    """vtkFilter for adding an integer value to vtkPolyData's
    CellData. Use SetValue to set the value to be added, or supply the
    optional argument 'Value' to the constructor.
    """
    def __init__(self, Value= -1):
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

class PolyDataCapAndLabeller1(vtkProgrammableFilter):
    """vtkFilter for capping a vtkPolyData surface, and labeling the
    cap with an integer cell data value.
    """
    def __init__(self, Value= -1):
        self.SetExecuteMethod(self._Execute)
        self.Value = Value
        return

    def SetValue(self, val):
        self.Value = val
        return
    def GetValue(self):
        return self.Value


    def _Execute(self, *args):
        # Setup the in/out PD
        input = self.GetPolyDataInput()

        output = self.GetPolyDataOutput()
        output.DeepCopy(input)
        outPoints = output.GetPoints()
        outCells = output.GetPolys()
        outCellData = output.GetCellData().GetScalars()

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

        lines = stripper.GetOutput()
        lineIndices = lines.GetLines().GetData()
        linePoints = lines.GetPoints()

        nCaps = lines.GetNumberOfLines()

        # Range of point indices for each strip, not including the last closing one
        pointIndexRanges = np.zeros((nCaps, 2), dtype=int)
        # Number of triangles that will need to be created for each cap
        nNewTris = np.zeros(nCaps, dtype=int)
        # Barycentres of faces
        barycentres = np.zeros((nCaps, 3), dtype=float)
        # Index of the point into the output vtkPoints 
        barycentreIds = np.zeros(nCaps, dtype=int)

        # This is locator we'll use to find the points of the boundaries 
        # the points list
        locator = vtkOctreePointLocator()
        locator.SetDataSet(input)
        locator.BuildLocator()

        # First count will be at index 0
        stripLenIdx = 0
        for iCap in xrange(nCaps):
            # The polyline is a number of points, followed by a list
            # of the indices of the points that make it up. The 
            # final index will be the same as the first one.
            stripLen = int(lineIndices.GetTuple1(stripLenIdx))

            nNewTris[iCap] = stripLen - 1
            # First point id after the length of the strip
            pointIndexRanges[iCap] = (stripLenIdx + 1, stripLenIdx + 1 + nNewTris[iCap])

            # Now compute barycentre of face
            # R = (Sum_i m_i r_i) / (Sum_i m_i)
            sum_miri = np.zeros(3, dtype=float)
            sum_mi = 0.

            # For efficiency, we will also construct the triangles
            # making up the face now. The barycentre will of course be
            # in the wrong place, but we'll fix that at the end
            barycentreIds[iCap] = outPoints.InsertNextPoint(0., 0., 0.)

            for iLine in xrange(*pointIndexRanges[iCap]):
                start = np.array(linePoints.GetPoint(int(lineIndices.GetTuple1(iLine))))
                end = np.array(linePoints.GetPoint(int(lineIndices.GetTuple1(iLine + 1))))
                ri = 0.5 * (start + end)
                # Use length as a proxy of mass
                mi = np.sqrt(np.dot(end - start, end - start))

                sum_mi += mi
                sum_miri += mi * ri

                # Create the triangle now
                outCells.InsertNextCell(3)
                startPointId = locator.FindClosestPoint(start)
                outCells.InsertCellPoint(startPointId)
                outCells.InsertCellPoint(barycentreIds[iCap])
                endPointId = locator.FindClosestPoint(end)
                outCells.InsertCellPoint(endPointId)

                outCellData.InsertNextTuple1(self.Value)

                continue

            barycentres[iCap] = sum_miri / sum_mi
            outPoints.SetPoint(barycentreIds[iCap], barycentres[iCap])

            stripLenIdx += stripLen + 1 # +1 for the count
            continue

        return

    pass

class PolyDataCapAndLabeller2(vtkProgrammableFilter):
    """vtkFilter for capping a vtkPolyData surface, and labeling the
    cap with an integer cell data value.
    """
    def __init__(self, Value= -1):
        self.SetExecuteMethod(self._Execute)
        self.Value = Value
        return

    def SetValue(self, val):
        self.Value = val
        return
    def GetValue(self):
        return self.Value


    def _Execute(self, *args):
        # Setup the in/out PD
        input = self.GetPolyDataInput()

        output = self.GetPolyDataOutput()
        output.DeepCopy(input)
        outPoints = output.GetPoints()
        outCells = output.GetPolys()
        outCellData = output.GetCellData().GetScalars()

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

        lines = stripper.GetOutput()
        lineIndices = lines.GetLines().GetData()
        linePoints = lines.GetPoints()

        nCaps = lines.GetNumberOfLines()

        # Range of point indices for each strip, not including the last closing one
        pointIndexRanges = np.zeros((nCaps, 2), dtype=int)
        # Number of triangles that will need to be created for each cap
        nNewTris = np.zeros(nCaps, dtype=int)
        # Barycentres of faces
        barycentres = np.zeros((nCaps, 3), dtype=float)
        # Index of the point into the output vtkPoints 
        barycentreIds = np.zeros(nCaps, dtype=int)

        # This is locator we'll use to find the points of the boundaries 
        # the points list
        locator = vtkOctreePointLocator()
        locator.SetDataSet(input)
        locator.BuildLocator()

        # First count will be at index 0
        stripLenIdx = 0
        for iCap in xrange(nCaps):
            # The polyline is a number of points, followed by a list
            # of the indices of the points that make it up. The 
            # final index will be the same as the first one.
            stripLen = int(lineIndices.GetTuple1(stripLenIdx))

            nNewTris[iCap] = stripLen - 1
            # First point id after the length of the strip
            pointIndexRanges[iCap] = (stripLenIdx + 1, stripLenIdx + 1 + nNewTris[iCap])

            # Now compute barycentre of face
            # R = (Sum_i m_i r_i) / (Sum_i m_i)
            sum_miri = np.zeros(3, dtype=float)
            sum_mi = 0.

            # For efficiency, we will also construct the triangles
            # making up the face now. The barycentre will of course be
            # in the wrong place, but we'll fix that at the end
            barycentreIds[iCap] = outPoints.InsertNextPoint(0., 0., 0.)

            for iLine in xrange(*pointIndexRanges[iCap]):
                start = np.array(linePoints.GetPoint(int(lineIndices.GetTuple1(iLine))))
                end = np.array(linePoints.GetPoint(int(lineIndices.GetTuple1(iLine + 1))))
                ri = 0.5 * (start + end)
                # Use length as a proxy of mass
                mi = np.sqrt(np.dot(end - start, end - start))

                sum_mi += mi
                sum_miri += mi * ri

                # Create the triangle now
                outCells.InsertNextCell(3)
                startPointId = locator.FindClosestPoint(start)
                outCells.InsertCellPoint(startPointId)
                outCells.InsertCellPoint(barycentreIds[iCap])
                endPointId = locator.FindClosestPoint(end)
                outCells.InsertCellPoint(endPointId)

                outCellData.InsertNextTuple1(self.Value)

                continue

            barycentres[iCap] = sum_miri / sum_mi
            outPoints.SetPoint(barycentreIds[iCap], barycentres[iCap])

            stripLenIdx += stripLen + 1 # +1 for the count
            continue

        return

    pass
class Face(object):
    def __init__(self, EdgePolyData):
        self.EdgePolyData = EdgePolyData
        self.nPoints = EdgePolyData.GetNumberOfPoints()
        
        self._Barycentre = None
        self._Perimeter = None
        self._MeanRadius = None
        self._FacePoints = None
        self._FacePolyData = None
        self._TriangulatedFace = None
        return
    
    @property
    def BaryCentre(self):
        if self._Barycentre is None:
            self._ComputeBarycentreAndPerimeter()
        return self._Barycentre
    
    @property
    def Perimeter(self):
        if self._Perimeter is None:
            self._ComputeBarycentreAndPerimeter()
        return self._Perimeter
    
    @property
    def MeanRadius(self):
        if self._MeanRadius is None:
            self._ComputeMeanRadius()
        return self._MeanRadius
    
    @property
    def TargetSideLength(self):
        # We have a target side length of twice the average side length.
        # This is because there will have been on average one extra point per 
        # pre-clipping triangle.
        return 2 * (self.Perimeter / self.nPoints)
    
    @property
    def NumberOfBands(self):
        return int(self.MeanRadius / (np.sqrt(3.) / 2. * self.TargetSideLength))
    
    @property
    def FacePoints(self):
        if self._FacePoints is None:
            self._ComputeFacePoints()
        return self._FacePoints
    
    @property
    def TriangulatedFace(self):
        if self._TriangulatedFace is None:
            self._TriangulateFace()
        return self._TriangulatedFace
    
    def _ComputeBarycentreAndPerimeter(self):
        # Now compute barycentre of face
        # R = (Sum_i m_i r_i) / (Sum_i m_i)
        # Sum_i m_i == perimeter
        self._Barycentre = np.zeros(3, dtype=float)
        self._Perimeter = 0.
        
        linePoints = self.EdgePolyData.GetPoints()
        lineIndices = self.EdgePolyData.GetLines().GetData()
        
        assert int(lineIndices.GetTuple1(0)) == self.nPoints + 1
        for iPoint in xrange(self.nPoints):
            # +1 since zeroth entry is number of points along the line
            iStart = int(lineIndices.GetTuple1(iPoint + 1))
            start = np.array(linePoints.GetPoint(iStart))
            iEnd = int(lineIndices.GetTuple1(iPoint + 2))
            end = np.array(linePoints.GetPoint(iEnd))
            # CofM of segment
            ri = 0.5 * (start + end)
            # Use length as a proxy of mass
            mi = np.sqrt(np.dot(end - start, end - start))
            
            self._Perimeter += mi
            self._Barycentre += mi * ri
            continue
        self._Barycentre /= self._Perimeter
        return
    
    def _ComputeMeanRadius(self):
        self._MeanRadius = 0.
        linePoints = self.EdgePolyData.GetPoints()
        lineIndices = self.EdgePolyData.GetLines().GetData()
        
        # Approximate the average radius of the gap
        for iPoint in xrange(self.nPoints):
            iStart = int(lineIndices.GetTuple1(iPoint + 1))
            pt = np.array(linePoints.GetPoint(iStart))
            self._MeanRadius += np.sqrt(np.sum((pt - self._Barycentre) ** 2))
            continue
        self._MeanRadius /= self.nPoints
        return
    
    def _ComputeFacePoints(self):
        self._FacePoints = vtkPoints()
        self._FacePoints.DeepCopy(self.EdgePolyData.GetPoints())
        linePoints = self.EdgePolyData.GetPoints()
        lineIndices = self.EdgePolyData.GetLines().GetData()

        rFrac = np.linspace(0., 1., self.NumberOfBands + 1)[1:-1]
        # Go round the cap again, adding new points
        # We want to add points in rings, spaced approximately evenly,
        # such that we get near equilateral triangles (except at the edge, 
        # where we know there are extra points due to the clipping).
        
        # There are nSides actual points, but we will work as if there are
        # (nSides/2). The perimeter is proportional to the radius, hence 
        # the number of new points is too. If rFrac is the reduced radius,
        # then:
        # nNew = (nSides/2) * rFrac
        nNew = ((self.nPoints / 2.) * rFrac).astype(int)
        
        # Easiest way to place these is by putting them rFrac of the way
        # from barycentre to the edge, every (nSides/nNew) edge sites.
        # This won't be perfect, but we don't care as we'll throw these 
        # points at the Delaunay triangulator in a moment 
        every = (2. / rFrac).astype(int)
        # Count how many we have actually added.
        nAdded = np.zeros(nNew.shape, dtype=int)
        bary = self.BaryCentre
        for iPoint in xrange(self.nPoints):
            which = np.where(iPoint % every == 0)[0]
            iStart = int(lineIndices.GetTuple1(iPoint + 1))
            pt = np.array(linePoints.GetPoint(iStart))
            
            for iBand in which:
                # If we've done enough, skip
                if nAdded[iBand] == nNew[iBand]: continue
                nAdded[iBand] += 1
                self._FacePoints.InsertNextPoint(rFrac[iBand] * (pt - bary) + bary)
                continue
                
            continue
        
        return
    
    def _TriangulateFace(self):
        # Now triangulate the face.
        # TODO: make this transform to the correct plane. This will work
        # for now as we're clipping in the XY plane
        triangulator = vtkDelaunay2D()
        inPD = vtkPolyData()
        inPD.SetPoints(self.FacePoints)
        triangulator.SetInput(inPD)
        triangulator.SetSource(self.EdgePolyData)
#        triangulator.BoundingTriangulationOn()
        triangulator.Update()
        self._TriangulatedFace = triangulator.GetOutput()
        return 
    
    pass

class PolyDataClipCapAndLabeller(vtkProgrammableFilter):
    """vtkFilter for capping a vtkPolyData surface, and labeling the
    cap with an integer cell data value.
    """
    def __init__(self, Value=None, ClipFunction=None):
        self.SetExecuteMethod(self._Execute)
        self.Value = Value
        self.ClipFunction = ClipFunction
        return
    
    def SetValue(self, val):
        self.Value = val
        return
    def GetValue(self):
        return self.Value
    
    def SetClipFunction(self, val):
        self.ClipFunction = val
        return
    def GetClipFunction(self):
        return self.ClipFunction
    
    def _MakeClipper(self, pdSource):
        clipper = vtkClipPolyData()
        clipper.SetInputConnection(pd)
        clipper.SetClipFunction(self.ClipFunction)
        # TODO: insert filter to get part closest to seed point
        return clipper.GetOutputPoint()
        
    def _GetOpenFaces(self, pd):
        # Gets any edges of the mesh
        edger = vtkFeatureEdges()
        edger.BoundaryEdgesOn()
        edger.FeatureEdgesOff()
        edger.NonManifoldEdgesOff()
        edger.ManifoldEdgesOff()
        edger.SetInput(pd)

        # Converts the edges to a polyline
        stripper = vtkStripper()
        stripper.SetInputConnection(edger.GetOutputPort())
        stripper.Update()
        
        lines = stripper.GetOutput()
        nLines = lines.GetNumberOfLines()
        
        # For easy cases, skip the copying
        if nLines == 0:
            return []
        elif nLines == 1:
            return [Face(lines)]
        
        ans = []
        allLinePoints = lines.GetPoints()
        allLineIndices = lines.GetLines().GetData()
        # First count will be at index 0
        lineLenIdx = 0
        for iLine in xrange(nLines):
            # The polyline is a number of points, followed by a list
            # of the indices of the points that make it up. The 
            # final index will be the same as the first one.
            newLineLen = int(allLineIndices.GetTuple1(lineLenIdx))
            
            # create points and connections for the line
            newLinePoints = vtkPoints()
            assert allLineIndices.GetTuple1(lineLenIdx + 1) == allLineIndices.GetTuple1(lineLenIdx + newLineLen), "Loop must be closed!"
            # -1 since it is a closed loop
            newLinePoints.SetNumberOfPoints(newLineLen - 1)
            newLineConnections = vtkCellArray()
            
            # This vtkIdTypeArray has to have the extra one to close up and first store the length
            newLineData = newLineConnections.GetData()
            newLineData.SetNumberOfTuples(newLineLen + 1)
            newLineData.SetTuple1(0, newLineLen)
            # Copy in the points and connections, going round one polyline
            for iNew in xrange(newLineLen - 1):
                iAll = int(allLineIndices.GetTuple1(lineLenIdx + 1 + iNew))
                newLinePoints.SetPoint(iNew, allLinePoints.GetPoint(iAll))
                newLineData.SetTuple1(iNew + 1, iNew)
                continue
            # Close the polyline
            newLineData.SetTuple1(newLineLen, 0)
            
            # Create polydata with the new points+connections
            newLinePolyData = vtkPolyData()
            newLinePolyData.SetPoints(newLinePoints)
            newLinePolyData.SetLines(newLineConnections)
            # create the Face and add to the output list
            ans.append(Face(newLinePolyData))
            
            lineLenIdx += newLineLen + 1 # +1 for the count
            continue
            
        return ans
    
    def _Execute(self, *args):
        assert isinstance(self.Value, int)
        assert isinstance(self.ClipFunction, vtkImplicitFunction)
        
        clippedSource = self._MakeClipper(self.GetInputConnection(0,0))
        clippedSource.Update()
        clipped = clippedSource.GetOutput()
        
        faces = self._GetOpenFaces(clipped)
        pdb.set_trace()
        
        output = self.GetPolyDataOutput()
        output.PrepareForNewData()
        
        for face in faces:
            
            
            continue

        return

    pass

class PolyLinesToCapConverter(vtkProgrammableFilter):
    """Given an input vtkPolyData object, containing the output of
    a vtkStripper (i.e. PolyLines), convert this to polygons
    covering the surface enclosed by the PolyLines.    
    """
    def __init__(self):
        self.SetExecuteMethod(self._Execute)

    def _Execute(self, *args):
        # Get the strips
        strips = self.GetPolyDataInput()

        # Get the out PD
        out = self.GetPolyDataOutput()
        out.PrepareForNewData()
        # Do a trick to set the points and poly of the cap
        out.SetPoints(strips.GetPoints())
        out.SetPolys(strips.GetLines())

        return
    pass


class XmlWriter(object):
    def __init__(self, profile):
        self.profile = profile
        return

    @staticmethod
    def indent(elem, level=0):
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                XmlWriter.indent(elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i
        return

    def Write(self):
        self.root = Element('hemelbsettings')
        self.DoSimulation()
        self.DoGeometry()
        self.DoIolets()
        self.DoVisualisation()

        self.indent(self.root)

        xmlFile = file(self.profile.OutputXmlFile, 'wb')
        xmlFile.write('<?xml version="1.0" ?>\n')
        ElementTree(self.root).write(xmlFile)
        return

    def DoSimulation(self):
        sim = SubElement(self.root, 'simulation')
        sim.set('cycles', str(3))
        sim.set('cyclesteps', str(1000))
        return

    def DoGeometry(self):
        geom = SubElement(self.root, 'geometry')
        geom.set('voxelsize', str(self.profile.VoxelSize))

        data = SubElement(geom, 'datafile')
        data.set('path', os.path.relpath(self.profile.OutputConfigFile,
                                         os.path.split(self.profile.OutputXmlFile)[0]))
        return

    def DoIolets(self):
        inlets = SubElement(self.root, 'inlets')
        outlets = SubElement(self.root, 'outlets')

        for io in self.profile.Iolets:
            if isinstance(io, Inlet):
                iolet = SubElement(inlets, 'inlet')
            elif isinstance(io, Outlet):
                iolet = SubElement(outlets, 'outlet')
            else:
                continue
            pressure = SubElement(iolet, 'pressure')
            pressure.set('mean', str(io.Pressure.x))
            pressure.set('amplitude', str(io.Pressure.y))
            pressure.set('phase', str(io.Pressure.z))

            normal = SubElement(iolet, 'normal')
            normal.set('x', str(io.Normal.x))
            normal.set('y', str(io.Normal.y))
            normal.set('z', str(io.Normal.z))

            position = SubElement(iolet, 'position')
            position.set('x', str(io.Centre.x))
            position.set('y', str(io.Centre.y))
            position.set('z', str(io.Centre.z))
            continue
        return

    def DoVisualisation(self):
        vis = SubElement(self.root, 'visualisation')

        centre = SubElement(vis, 'centre')
        centre.set('x', '0.0')
        centre.set('y', '0.0')
        centre.set('z', '0.0')

        orientation = SubElement(vis, 'orientation')
        orientation.set('longitude', '45.0')
        orientation.set('latitude', '45.0')

        display = SubElement(vis, 'display')
        display.set('zoom', '1.0')
        display.set('brightness', '0.03')

        range = SubElement(vis, 'range')
        range.set('maxvelocity', str(0.1))
        range.set('maxstress', str(0.1))

        return

    pass
