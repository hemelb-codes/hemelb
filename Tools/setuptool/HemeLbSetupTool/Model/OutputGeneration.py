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

import pdb

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
        # Pull in the StlReader too
        self.StlReader = profile.StlReader

        self.ClippedSurface = self.ConstructClipPipeline()
        self.Locator = vtkOBBTree()
        self.Locator.SetTolerance(0.)

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

        i = -1
        j = -1
        for block in domain.SmartIterBlocks():
            i += 1
            # Open the BlockStarted context of the writer; this will
            # deal with flushing the state to the file (or not, in the
            # case where there are no fluid sites).
            with writer.BlockStarted() as blockWriter:
                for site in block.IterSites():
                    j += 1

                    self.ClassifySite(site)
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

    def IterHitsForLink(self, start, end):
        """Given a pair of sites, yield the itersections one at a
        time.
        """
        # Intersect against the STL surface
        stlHits = StlHitList(self.Locator, start, end)

        # Now start yielding intersections
        topStl = stlHits.next()

        while topStl is not None:
            yield topStl
            topStl = stlHits.next()
            continue

        return

    def ClassifySite(self, site):
        """Perform classification of the supplied sites. Note that
        this will alter the connected sites that have yet to be
        classified, as we wish to examine each link only once.

        Each site must have its IsFluid and IsEdge flags set, along
        with the CutDistances array (at appropriate indices), {Wall,
        Boundary}x{Distance, Normal}.
        """
        if site.IsFluid is None:
            # We're the first site; if we're not, the IsFluid flag
            # will have been set to True/False below
            if self.Locator.InsideOrOutside(site.Position) < 0:
                # vtkOBBTree.InsideOrOutside returns -1 for inside
                site.IsFluid = True
            else:
                site.IsFluid = False
                pass
            pass

        for i, neigh in site.EnumerateLaterNeighbours():
            # Check our neighbours, who are further on in what would
            # be a conventional array of all sites.

            nHits = 0

            for hitPoint, hitObj in self.IterHitsForLink(site.Position, neigh.Position):
                nHits += 1
                if nHits == 1:
                    # First hit, assign the distance (in STL units for
                    # now) to first intersection (i.e. CutDistance)
                    # and the ID of the vtkPolygon which we
                    # intersected (CutCellId). We also now know we're
                    # an edge site.
                    
                    site.CutDistances[i] = sqrt(np.sum((hitPoint -
                                                        site.Position) ** 2))
                    site.IsEdge = True
                    site.CutCellIds[i] = hitObj

                    pass

                continue

            # If we had any hits, we need to set the last one for the
            # reverse link; hitPoint and hitObj assigned above will
            # handily have the final values from the loop.
            if nHits > 0:
                neigh.IsEdge = True
                neigh.CutDistances[i] = sqrt(np.sum((hitPoint - neigh.Position) ** 2))
                neigh.CutCellIds[i] = hitObj
                pass

            # Set the neighbour's fluid flag
            if nHits % 2 == 0:
                # Even nHits => we crossed the surface and then back
                # to where we were.
                neigh.IsFluid = site.IsFluid
            else:
                # Odd nHits => we're now on the other side of the
                # surface.
                neigh.IsFluid = not site.IsFluid
            continue

        if not site.IsFluid or not site.IsEdge:
            # Nothing more to do for solid sites or simple fluid sites
            return
        
        # These CutDistances need to be fractions of the corresponding
        # lattice vector, rather than distances in lattice units or
        # physical units. HOWEVER, we need to know the distances to
        # the walls and Iolets in plain lattice units below, so only
        # scale by the VoxelSize for now and scale to vector fractions
        # once we're done.
        site.CutDistances /= self.VoxelSize

        # Get the normals and scalars associated with the surface
        # polygons. The scalars hold the index of the Iolet which they
        # represent, -1 meaning they aren't an Iolet.
        celldata = self.ClippedSurface.GetOutput().GetCellData()
        normals = celldata.GetNormals()
        ioletIds = celldata.GetScalars()

        site.IsEdge = False
        for i, hitCellId in enumerate(site.CutCellIds):
            # The site.CutCellsIds array is initialised to -1 and then
            # updated with the index of the first vtkPolygon they
            # intersect.
            if hitCellId == -1:
                # We didn't hit in this direction.
                continue
            
            ioletId = ioletIds.GetValue(hitCellId)
            if ioletId >= 0:
                # It's an iolet
                if site.CutDistances[i] < site.BoundaryDistance:
                    # It is the closest yet
                    io = site.Iolet = self.Iolets[ioletId]
                    # TODO: confirm this should be in lattice units 
                    site.BoundaryDistance = site.CutDistances[i]
                    site.BoundaryNormal[:] = io.Normal.x, io.Normal.y, io.Normal.z
                    site.BoundaryId = io.Index
                    pass
                
            else:
                # It's wall
                site.IsEdge = True
                if site.CutDistances[i] < site.WallDistance:
                    # If it's the closest point yet, store it
                    site.WallDistance = site.CutDistances[i]
                    site.WallNormal[:] = normals.GetTuple3(hitCellId)
                    pass
                pass

            continue

        # Scale to be fractions of lattice vectors instead of
        # distances in lattice units; see above for more.
        site.CutDistances /= LatticeSite.latticeVectorNorms
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
        adder.SetInputConnection(self.StlReader.GetOutputPort())

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

class StlHitList(object):
    """Object that finds the intersections given by a vtkLocator the
    line joining a start and end LatticeSite. It then gives a simple
    iterator-like method (next) for getting these points, which when
    it finishes, returns None rather than raise StopIteration.
    """
    def __init__(self, locator, start, end):
        """locator - vtkLocator to use to test for intersections.
        start - starting LatticeSite
        end - ending LatticeSite
        """
        self.HitPoints = vtkPoints()
        self.HitCellIds = vtkIdList()
        # Find intersections
        locator.IntersectWithLine(start, end, self.HitPoints, self.HitCellIds)
        # Cache the number of them
        self.len = self.HitPoints.GetNumberOfPoints()
        # Track where we are in the iteration
        self.i = 0
        return

    def next(self):
        """Move our internal internal counter on one and return the
        new value, or None in the case that it doesn't exist.
        """
        i = self.i
        if i < self.len:
            self.i += 1
            return self[i]
        else:
            return None

    def __len__(self):
        return self.len

    def __getitem__(self, i):
        if i < 0: raise IndexError
        if i >= self.len: raise IndexError

        hit = self.HitPoints.GetPoint(i)

        return hit, self.HitCellIds.GetId(i)
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
        self.SetExecuteMethod(self.Execute)

    def Execute(self, *args):
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

class PolyLinesToRefinedCapConverter(vtkProgrammableFilter):
    """Given an input vtkPolyData object, containing the output of
    a vtkStripper (i.e. PolyLines), convert this to polygons
    covering the surface enclosed by the PolyLines.    
    """
    def __init__(self):
        self.SetExecuteMethod(self.Execute)

    def Execute(self, *args):
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
        self.OutputConfigFile = OutputConfigFile
        self.StressType = StressType
        self.BlockSize = BlockSize
        self.BlockCounts = BlockCounts
        # Make sure this in in metres
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
        # Position of site index (0,0,0) in block index (0,0,0), in metres in the STL file's coordinate system
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
        head = self.HeaderEncoder.get_buffer()
        assert len(head) == (self.bodyStart - self.headerStart)

        self.file.seek(self.headerStart)
        self.file.write(head)
        self.file.close()
        return

    @contextmanager
    def BlockStarted(self):
        encoder = xdrlib.Packer()


        def incrementor():
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


class Domain(object):
    def __init__(self, VoxelSize, SurfaceBounds, BlockSize=8):
        self.VoxelSize = VoxelSize
        self.BlockSize = BlockSize
        # VTK standard order of (x_min, x_max, y_min, y_max, z_min, z_max)
        bb = np.array(SurfaceBounds)
        bb.shape = (3, 2)

        origin = []
        blocks = []
        for min, max in bb:
            size = max - min
            # int() truncates, we add 2 to make sure there's enough
            # room for the sites just outside.
            nSites = int(size / VoxelSize) + 2

            # The extra space
            extra = nSites * VoxelSize - size
            # We want to balance this equally with the placement of
            # the first site.
            siteZero = min - 0.5 * extra

            nBlocks = nSites / BlockSize
            remainder = nSites % BlockSize
            if remainder:
                nBlocks += 1
                pass

            origin.append(siteZero)
            blocks.append(nBlocks)
            continue

        self.Origin = np.array(origin)
        self.BlockCounts = np.array(blocks)

#        for ijk, val in np.ndenumerate(self.blocks):
#            self.blocks[ijk] = MacroBlock(ijk=ijk, size=BlockSize)
#            continue

        return

    def CalcPositionFromIndex(self, index):
        return self.Origin + self.VoxelSize * np.array(index)

    def GetBlock(self, *blockIjk):
        val = self.blocks[blockIjk]
        if val is None:
            val = self.blocks[blockIjk] = MacroBlock(domain=self, ijk=blockIjk, size=self.BlockSize)
            pass
        return val

    def GetSite(self, *globalSiteIjk):
        blockIjk = [i / self.BlockSize for i in globalSiteIjk]
        block = self.GetBlock(*blockIjk)
        return block.GetSite(*globalSiteIjk)

    def SmartIterBlocks(self):
        # Fill the blocks with Nones
        # TODO: make this choose the best ordering of indices for memory efficiency
        self.blocks = np.empty(self.BlockCounts, dtype=object)
        maxInds = [l - 1 for l in self.BlockCounts]

        for ijk, val in np.ndenumerate(self.blocks):
            if val is None:
                # If the block hasn't been created, do so
                val = self.blocks[ijk] = MacroBlock(ijk=ijk, domain=self, size=self.BlockSize)
                pass

            yield val

            # Delete any unnecessary blocks now
            for i in range(ijk[0] - 1, ijk[0] + 1):
                if i < 0: continue
                if i == ijk[0] and i != maxInds[0]: continue
                for j in range(ijk[1] - 1, ijk[1] + 1):
                    if j < 0: continue
                    if j == ijk[1] and j != maxInds[1]: continue
                    for k in range(ijk[2] - 1, ijk[2] + 1):
                        if k < 0: continue
                        if k == ijk[2] and k != maxInds[2]: continue
                        self.blocks[i, j, k] = None
                        continue
                    continue
                continue


            continue # nd.enumerate(self.blocks)

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

class MacroBlock(object):
    def __init__(self, domain=None, ijk=(0, 0, 0), size=8):
        self.domain = domain
        self.size = size
        self.ijk = (ijk)
        self.sites = np.empty((size, size, size), dtype=object)
        for localSiteIjk, ignored in np.ndenumerate(self.sites):
            globalSiteIjk = tuple(self.ijk[i] * self.size + localSiteIjk[i]
                                  for i in xrange(3))

            self.sites[localSiteIjk] = LatticeSite(block=self, ijk=globalSiteIjk)
        return

    def IterSites(self):
        for site in self.sites.flat:
            yield site
            continue
        return

    def NdEnumerateSites(self):
        for item in np.ndenumerate(self.sites):
            yield item
            continue
        return

    def GetLocalSite(self, *localSiteIjk):
        return self.sites[localSiteIjk]

    def GetSite(self, *globalSiteIjk):
        localSiteIjk = tuple(globalSiteIjk[i] - self.ijk[i] * self.size
                             for i in xrange(3))

        # Check if the coords belong to another block, i.e. any of
        # the local ones outside the range [0, self.size)
        if any(map(lambda x: x < 0 or x >= self.size, localSiteIjk)):
            return self.domain.GetSite(*globalSiteIjk)

        return self.GetLocalSite(*localSiteIjk)

    def CalcPositionFromIndex(self, index):
        return self.domain.CalcPositionFromIndex(index)
    pass

class LatticeSite(object):
    # This ordering of the lattice vectors must be the same as in the HemeLB source.
    neighbours = np.array([[ 1, 0, 0],
                           [-1, 0, 0],
                           [ 0, 1, 0],
                           [ 0, -1, 0],
                           [ 0, 0, 1],
                           [ 0, 0, -1],
                           [ 1, 1, 1],
                           [-1, -1, -1],
                           [ 1, 1, -1],
                           [-1, -1, 1],
                           [ 1, -1, 1],
                           [-1, 1, -1],
                           [ 1, -1, -1],
                           [-1, 1, 1]])
    latticeVectorNorms = np.sqrt(np.sum(neighbours ** 2, axis= -1))
    laterNeighbourInds = [0, 2, 4, 6, 8, 10, 12]
    def __init__(self, block=None, ijk=(0, 0, 0)):
        self.block = block
        self.ijk = ijk
        self.Position = block.CalcPositionFromIndex(ijk)

        self.IsFluid = None
        self.IsEdge = None
        self.Iolet = None
        self.BoundaryId = None

        # Attributes that will be updated by Profile.ClassifySite
#        self.Type = None
        self.BoundaryNormal = np.zeros(3)
        self.BoundaryDistance = float('inf')

        self.WallNormal = np.zeros(3)
        self.WallDistance = float('inf')

        self.CutDistances = 1. / np.zeros(len(self.neighbours), dtype=float)
        self.CutCellIds = -np.ones(len(self.neighbours), dtype=int)
        return

    @property
    def Type(self):
        if self.IsFluid:
            if isinstance(self.Iolet, Inlet):
                return INLET_TYPE
            elif isinstance(self.Iolet, Outlet):
                return OUTLET_TYPE
            else:
                return FLUID_TYPE
        else:
            return SOLID_TYPE

    @property
    def Config(self):
        cfg = self.Type
        if not self.IsFluid:
            return self.Type

        # Fluid sites now
        if self.Type == FLUID_TYPE and not self.IsEdge:
            # Simple fluid sites
            return self.Type

        # A complex one

        if self.IsEdge:
            # Bit fiddle the boundary config. See comment below
            # by BOUNDARY_CONFIG_MASK for definition.
            boundary = 0
            for i, neigh in self.EnumerateNeighbours():
                if not neigh.IsFluid:
                    # If the lattice vector is cut, set the flag
                    boundary |= 1 << i
                    pass
                continue
            # Shift the boundary bit field to the appropriate
            # place and set these bits in the type
            cfg |= boundary << BOUNDARY_CONFIG_SHIFT

            if boundary:
                # Set this bit if we've hit any solid sites
                cfg |= PRESSURE_EDGE_MASK
                pass
            pass

        if self.Type != FLUID_TYPE:
            # It must be an inlet or outlet
            # Shift the index left and set the bits
            cfg |= self.BoundaryId << BOUNDARY_ID_SHIFT
            pass

        return cfg

    def EnumerateNeighbours(self):
        """Create an iterator that enumerates our neighbours.
        """
        domain = self.block.domain
        shape = domain.BlockCounts * domain.BlockSize
        for i, neigh in enumerate(self.neighbours):
            ind = self.ijk + neigh
            # Skip any out of range
            if np.any(ind < 0) or np.any(ind >= shape):
                continue

            yield i, self.block.GetSite(ind[0], ind[1], ind[2])
            continue
        return

    def EnumerateLaterNeighbours(self):
        """Create an iterator that enumerates our neighbours who're later in the array.
        """
        domain = self.block.domain
        shape = domain.BlockCounts * domain.BlockSize
        for i in self.laterNeighbourInds:
            neigh = self.neighbours[i]
            ind = self.ijk + neigh
            # Skip any out of range
            if np.any(ind < 0) or np.any(ind >= shape):
                continue

            yield i, self.block.GetSite(ind[0], ind[1], ind[2])
            continue
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
