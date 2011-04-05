from contextlib import contextmanager
import xdrlib
import numpy as np
from math import sqrt
import bisect
import os.path
from xml.etree.ElementTree import Element, SubElement, ElementTree

from vtk import vtkClipPolyData, vtkAppendPolyData, vtkPlane, \
    vtkStripper, vtkFeatureEdges, vtkPolyDataConnectivityFilter, \
    vtkPoints, vtkIdList, vtkPolyDataNormals, vtkProgrammableFilter, \
    vtkOBBTree, vtkTriangleFilter, vtkCleanPolyData, vtkIntArray

from .Iolets import Inlet, Outlet, Iolet
from .OutputGenerationHelpers import Domain, DomainSmartBlockIterator, ClassifySite

import pdb

np.seterr(divide='ignore')

class ConfigGenerator(object):
    """This object is in charge of creating the input for HemeLB from
    the supplied Model.Profile object. The process is coordinated by
    Execute.
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
        self.Locator = vtkOBBTree()
        self.Locator.SetTolerance(0.)
        
        self.IsFirstSite = True
        return

    def Execute(self):
        """Create the output based on our configuration.
        
        1) Add an index to the Iolets within their type.
        
        2) Clip and cap the STL model; this is performed by a VTK
        pipeline (see ConstructClipPipeline).

        3) Create the Writer and Domain objects. The Domain initially
        has no MacroBlocks allocated, only a 3D array of Nones of the
        appropriate size.
        
        4) Iterate over the set of MacroBlocks (and their contained
        LatticeSites), constructing them in a lazy fashion; the blocks
        are deleted when they can no longer be used to save memory,
        see Domain.SmartIterBlocks for details.

        5) Each site is fully classified (see ClassifySite) when its
        turn comes to be written. The classification takes care to
        examine each link between sites only once for efficiency. It
        therefore partially updates the other site's attributes.
        
        The Writer object deals with writing the preamble and header
        of the output.
        """
        # Work out the indices of the IOlets 
        nIn = 0
        nOut = 0
        for io in self.Iolets:
            if isinstance(io, Inlet):
                io.Index = nIn
                nIn += 1
            elif isinstance(io, Outlet):
                io.Index = nOut
                nOut += 1
                pass
            continue

        # Run the pipeline as we need it to build the locator
        self.ClippedSurface.Update()
        surface = self.ClippedSurface.GetOutput()

        # Create the locator
        self.Locator.SetDataSet(surface)
        self.Locator.BuildLocator()


        domain = Domain(self.VoxelSize, surface.GetBounds())

        writer = Writer(OutputConfigFile=self.OutputConfigFile,
                        StressType=self.StressType,
                        BlockSize=domain.BlockSize,
                        BlockCounts=domain.BlockCounts,
                        VoxelSize=domain.VoxelSize,
                        Origin=domain.Origin)

        for block in DomainSmartBlockIterator(domain):
            # Open the BlockStarted context of the writer; this will
            # deal with flushing the state to the file (or not, in the
            # case where there are no fluid sites).
            with writer.BlockStarted() as blockWriter:
                for site in block.IterSites():
                    ClassifySite(self, site)
                    # cache the type cos it's probably slow to compute
                    type = site.Type
                    cfg = site.Config
                    blockWriter.pack_uint(cfg)
                    
                    if type == SOLID_TYPE:
                        # Solid sites, we don't do anything
                        continue

                    # Increase count of fluid sites
                    blockWriter.IncrementFluidSitesCount()

                    if cfg == FLUID_TYPE:
                        # Pure fluid sites don't need any more data
                        continue

                    # It must be INLET/OUTLET/EDGE

                    if type == INLET_TYPE or type == OUTLET_TYPE:
                        for n in site.BoundaryNormal: blockWriter.pack_double(n)
                        blockWriter.pack_double(site.BoundaryDistance)
                        pass

                    if site.IsEdge:
                        for n in site.WallNormal: blockWriter.pack_double(n)
                        blockWriter.pack_double(site.WallDistance)

                    for cd in site.CutDistances: blockWriter.pack_double(cd)

                    continue

                pass
            continue # for block in domain...

        writer.Close()

        # Write the XML file
        XmlWriter(self).Write()
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

class Writer(object):
    """Generates the HemeLB input file.
    
    
    """
    def __init__(self, OutputConfigFile=None, StressType=None,
                 BlockSize=8, BlockCounts=None, VoxelSize=None, Origin=None):
        """This opens the file and writes the preamble and a dummy
        header. All keyword arguments are required.
        
        OutputConfigFile - path to output file
        
        StressType - an integer indicating what stress calculation
        method to use

        BlockSize - number of sites along one dimension of a block

        BlockCounts - number of blocks along x-, y- and z-directions

        VoxelSize - length of a voxel, in metres

        Origin - location of the (0,0,0) site in physical coordinates
        
        """
        self.OutputConfigFile = OutputConfigFile
        self.StressType = StressType
        self.BlockSize = BlockSize
        self.BlockCounts = BlockCounts
        self.VoxelSize = VoxelSize
        # Make sure this in metres
        self.Origin = Origin

        # Truncate, and open for read & write in binary mode
        self.file = file(OutputConfigFile, 'w+b')

        encoder = xdrlib.Packer()
        # Write the preamble, starting with the stress type
        encoder.pack_int(StressType)

        # Blocks in each dimension
        for count in BlockCounts:
            encoder.pack_uint(int(count))
            continue
        # Sites along 1 dimension of a block
        encoder.pack_uint(BlockSize)

        # Voxel Size, in metres
        encoder.pack_double(VoxelSize)
        # Position of site index (0,0,0) in block index (0,0,0), in
        # metres in the STL file's coordinate system
        for ori in Origin:
            encoder.pack_double(ori)
            continue

        # Write this to the file
        self.file.write(encoder.get_buffer())

        # Reset, we're going to write a dummy header now
        encoder.reset()
        # For each block
        for ignored in xrange(np.prod(BlockCounts)):
            encoder.pack_uint(0) # n fluid sites
            encoder.pack_uint(0) # n bytes
            continue
        # Note the start of the header
        self.headerStart = self.file.tell()
        # Write it
        self.file.write(encoder.get_buffer())
        # Note the start of the body
        self.bodyStart = self.file.tell()

        self.HeaderEncoder = xdrlib.Packer()
        return

    def Close(self):
        """Rewrites the header map and closes the file.
        """
        head = self.HeaderEncoder.get_buffer()
        assert len(head) == (self.bodyStart - self.headerStart)

        self.file.seek(self.headerStart)
        self.file.write(head)
        self.file.close()
        return

    @contextmanager
    def BlockStarted(self):
        """A context manager (see PEP 343) for use in a 'with'
        statement.
        
        This context manager will return an enhanced xdrlib.Packer
        instance which should be used to write the data for a single
        MacroBlock within the body of the 'with' statement. Every time
        a fluid site is written, the IncrementFluidSitesCount method
        of the Packer should be called. This is used at the end of the
        'with' block to write a record in the header.
         
        """
        # The encode we will yiled
        encoder = xdrlib.Packer()

        # This function will be set to the IncrementFluidSitesCount
        # attribute of the Packer
        def incrementor():
            """Increase the count of fluid sites that have been
            written by the active BlockStarted context.
            """
            incrementor.count += 1
            return
        incrementor.count = 0
        
        encoder.IncrementFluidSitesCount = incrementor
        # Give the altered XDR encoder back to our caller
        yield encoder
        
        # Write our record into the header buffer
        self.HeaderEncoder.pack_uint(incrementor.count)

        if incrementor.count == 0:
            # We mustn't write anything to the file, so note that it
            # takes zero bytes
            self.HeaderEncoder.pack_uint(0)
        else:
            # Write the block to the main file 

            blockStart = self.file.tell()
            self.file.write(encoder.get_buffer())
            blockEnd = self.file.tell()
            self.HeaderEncoder.pack_uint(blockEnd - blockStart)
            pass

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



# Site types
SOLID_TYPE = 0b00
FLUID_TYPE = 0b01
INLET_TYPE = 0b10
OUTLET_TYPE = 0b11

BOUNDARIES = 3
INLET_BOUNDARY = 0
OUTLET_BOUNDARY = 1
WALL_BOUNDARY = 2

SITE_TYPE_BITS = 2
BOUNDARY_CONFIG_BITS = 14
BOUNDARY_DIR_BITS = 4
BOUNDARY_ID_BITS = 10

BOUNDARY_CONFIG_SHIFT = SITE_TYPE_BITS;
BOUNDARY_DIR_SHIFT = BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
BOUNDARY_ID_SHIFT = BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

#===============================================================================
# Horrifying bit-fiddling masks courtesy of Marco.
# Comments show the bit patterns.
#===============================================================================

SITE_TYPE_MASK = ((1 << SITE_TYPE_BITS) - 1)
# 0000 0000  0000 0000  0000 0000  0000 0011
# These give the *_TYPE given above

BOUNDARY_CONFIG_MASK = ((1 << BOUNDARY_CONFIG_BITS) - 1) << BOUNDARY_CONFIG_SHIFT
# 0000 0000  0000 0000  1111 1111  1111 1100
# These bits are set if the lattice vector they correspond to takes one to a solid site
# The following hex digits give the index into LatticeSite.neighbours
# ---- ----  ---- ----  DCBA 9876  5432 10--

BOUNDARY_DIR_MASK = ((1 << BOUNDARY_DIR_BITS) - 1) << BOUNDARY_DIR_SHIFT
# 0000 0000  0000 1111  0000 0000  0000 0000
# No idea what these represent. As far as I can tell, they're unused.

BOUNDARY_ID_MASK = ((1 << BOUNDARY_ID_BITS) - 1) << BOUNDARY_ID_SHIFT
# 0011 1111  1111 0000  0000 0000  0000 0000
# These bits together give the index of the inlet/outlet/wall in the output XML file

PRESSURE_EDGE_MASK = 1 << 31
# 1000 0000  0000 0000  0000 0000  0000 0000
