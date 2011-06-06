import numpy as np
import os.path
from xml.etree.ElementTree import Element, SubElement, ElementTree

from vtk import vtkClipPolyData, vtkAppendPolyData, vtkPlane, \
    vtkStripper, vtkFeatureEdges, vtkPolyDataConnectivityFilter, \
    vtkPolyDataNormals, vtkProgrammableFilter, \
    vtkOBBTree, vtkTriangleFilter, vtkCleanPolyData, vtkIntArray

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
        
#        self.generator.SetFullSurfaceSource(profile.SurfaceSource)
        
        clippedSrc = Clipper(profile).ConstructClipPipeline()
        self.generator.SetClippedSurfaceSource(clippedSrc)
        
        # Create the locator
        clippedSrc.Update()
        clippedSurface = clippedSrc.GetOutput()
        self.Locator = vtkOBBTree()
        self.Locator.SetDataSet(clippedSurface)
        self.Locator.BuildLocator()
        # Again, we need to make sure this is kept alive and not GC'ed
        self.generator.SetLocator(self.Locator)
        return
    
    def Execute(self):
        """Forward this to the C++ implementation.
        """
#        self.generator.Execute()
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

        self.ClippedSurface = self.ConstructClipPipeline()
        return

    def ConstructClipPipeline(self):
        """This constructs a VTK pipeline to clip the vtkPolyData read
        from the STL file against the Iolets. It also adds a scalar
        value to each polygon indicating which Iolet it represents, or
        -1 if it is just a wall, and computes polygon normals.
        
        Note that this does NOT EXECUTE the pipeline.
        """
        # Add the Iolet id -1 to all cells first
        adder = IntegerAdder(Value= -1)
        adder.SetInputConnection(self.SurfaceSource.GetOutputPort())

        # Have the name pdSource first point to the input, then loop
        # over IOlets, clipping and capping.
        pdSource = adder
        for i, iolet in enumerate(self.Iolets):
#            pdb.set_trace()

            # TODO: switch from this simple ImplicitPlane (i.e.
            # infinite plane) clipping to excising a thin cuboid
            # from the Iolet and keeping the piece next to the
            # SeedPoint.

            # Create a plane to clip against.
            plane = vtkPlane()
            plane.SetOrigin(iolet.Centre.x, iolet.Centre.y, iolet.Centre.z)
            plane.SetNormal(iolet.Normal.x, iolet.Normal.y, iolet.Normal.z)

            # Clips to give the surface we want
            clipper = vtkClipPolyData()
            clipper.SetClipFunction(plane)
            clipper.SetInputConnection(pdSource.GetOutputPort())

            # Gets only the parts that are connected to the seed
            # point.
            filter = vtkPolyDataConnectivityFilter()
            filter.SetInputConnection(clipper.GetOutputPort())
            filter.SetClosestPoint(self.SeedPoint.x, self.SeedPoint.y, self.SeedPoint.z)
            filter.SetExtractionModeToClosestPointRegion()

            # Gets any edges of the mesh
            edger = vtkFeatureEdges()
            edger.BoundaryEdgesOn()
            edger.FeatureEdgesOff()
            edger.NonManifoldEdgesOff()
            edger.ManifoldEdgesOff()
            edger.SetInputConnection(filter.GetOutputPort())
            # Converts the edges to a polyline
            stripper = vtkStripper()
            stripper.SetInputConnection(edger.GetOutputPort())
            stripper.Update()
            # Turns the poly line to a polygon
            faceMaker = PolyLinesToCapConverter()
            faceMaker.SetInputConnection(stripper.GetOutputPort())
            # Triangulates it
            triangulator = vtkTriangleFilter()
            triangulator.SetInputConnection(faceMaker.GetOutputPort())

            # Adds the index of the iolet to the triangles
            idAdder = IntegerAdder(Value=i)
            idAdder.SetInputConnection(triangulator.GetOutputPort())
            
            # Joins the two PolyData together
            joiner = vtkAppendPolyData()
            joiner.AddInputConnection(filter.GetOutputPort())
            joiner.AddInputConnection(idAdder.GetOutputPort())
            
            # Set the source of the next iteraction to the capped
            # surface producer.
            pdSource = joiner
            continue

        # This will tidy up the final PolyData, merging points etc to
        # make sure it's a closed surface.
        cleaner = vtkCleanPolyData()
        cleaner.SetTolerance(0.)
        cleaner.SetInputConnection(pdSource.GetOutputPort())

        # Adds cell normals to the PolyData
        normer = vtkPolyDataNormals()
        normer.SetInputConnection(cleaner.GetOutputPort())
        normer.ComputeCellNormalsOn()
        normer.ComputePointNormalsOff()
        
        return normer

    pass


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
