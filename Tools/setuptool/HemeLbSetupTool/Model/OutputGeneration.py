from contextlib import contextmanager
import xdrlib
import numpy as np

from vtk import vtkClipPolyData, vtkAppendPolyData, vtkPlane, vtkStripper, vtkPolyData, \
    vtkFeatureEdges, vtkPolyDataConnectivityFilter, vtkSelectEnclosedPoints, vtkCellLocator, \
    vtkPoints, vtkIdList, vtkPolyDataNormals, vtkProgrammableFilter
import pdb

class ConfigGenerator(object):
    def __init__(self, profile):
        # Copy the profile's _Args
        for k in profile._Args:
            setattr(self, k, getattr(profile, k))
            continue
        # Pull in the StlReader too
        self.SurfaceSource = profile.StlReader
        return
    
    def Execute(self):
        """Create the output based on our configuration.
        """
        clipper = self.ConstructClipPipeline()
        clipper.Update()
        # Run the pipeline as we need it to build the locator
        surface = clipper.GetOutput()
        
        self.Locator = locator = vtkCellLocator()
        locator.SetDataSet(surface)
        locator.BuildLocator()
        
        checker = self.IsEnclosedChecker = vtkSelectEnclosedPoints()
        checker.Initialize(surface)
        
        domain = Domain(self.VoxelSize, surface.GetBounds())
        self.ExternalPoint = domain.Origin
        
        writer = Writer(OutputConfigFile=self.OutputConfigFile,
                        StressType=self.StressType,
                        BlockSize=domain.BlockSize,
                        BlockCounts=domain.BlockCounts)
        
        
        for block in domain.SmartIterBlocks():
            # Open the BlockStarted context of the writer; this will
            # deal with flushing the state to the file (or not, in
            # the case where there are no fluid sites).
            with writer.BlockStarted() as blockWriter:
                for site in block.IterSites():
                    self.ClassifySite(site)
                    
                    blockWriter.pack_uint(site.Type)
                    if site.Type == SOLID:
                        # Solid sites, we don't do anything
                        continue
                    # Increase count of fluid sites
                    blockWriter.IncrementFluidSitesCount()
                    
                    if site.Type == FLUID:
                        # Pure fluid sites don't need any more data
                        continue
                    
                    # It must be INLET/OUTLET/EDGE
                    for n in site.Normal: blockWriter.pack_double(n)
                    blockWriter.pack_double(site.BoundaryDistance)
                    
                    for cd in site.CutDistances: blockWriter.pack_double(cd)
                    
                    continue
                
                pass
            continue # for block in domain...
        
        writer.RewriteHeader()
        return

    def IsSiteEnclosed(self, site):
        if site.IsFluid is not None:
            return site.IsFluid
        
        ans = site.IsFluid = self.IsEnclosedChecker.IsInsideSurface(site.Position)
        return ans
    
    def ClassifySite(self, site):
        
        if not self.IsSiteEnclosed(site):
            return
        
        for i, neigh in enumerate(site.IterNeighbours()):
            if not self.IsSiteEnclosed(neigh):
                # We're an edge
                site.IsEdge = True
                site.CutDistances[i] = 1.0
                pass
            continue
        
        return
    
    def ConstructClipPipeline(self):
        """Clip the PolyData read from the STL file against the planes.
        """
        # Have the name pdSource first point to the input, then loop
        # over IOlets, clipping and capping.
        pdSource = self.SurfaceSource
        
        for iolet in self.Iolets:
#            pdb.set_trace()
            plane = vtkPlane()
            plane.SetOrigin(iolet.Centre.x, iolet.Centre.y, iolet.Centre.z)
            plane.SetNormal(iolet.Normal.x, iolet.Normal.y, iolet.Normal.z)
            # TODO: switch from this simple ImplicitPlane (i.e.
            # infinite plane) clipping to excising a thin cuboid
            # from the Iolet and keeping the piece next to the
            # SeedPoint.
             
            # Clip: gives the surface we want
            clipper = vtkClipPolyData()
            clipper.SetClipFunction(plane)
            clipper.SetInputConnection(pdSource.GetOutputPort())
            
            filter = vtkPolyDataConnectivityFilter()
            filter.SetInputConnection(clipper.GetOutputPort())
            filter.SetClosestPoint(self.SeedPoint.x, self.SeedPoint.y, self.SeedPoint.z)
            filter.SetExtractionModeToClosestPointRegion()
            
            # Get any edges of the mesh, will give us a PolyLine
            edger = vtkFeatureEdges()
            edger.BoundaryEdgesOn()
            edger.FeatureEdgesOff()
            edger.NonManifoldEdgesOff()
            edger.ManifoldEdgesOff()
            edger.SetInputConnection(filter.GetOutputPort())
            
            stripper = vtkStripper()
            stripper.SetInputConnection(edger.GetOutputPort())
            stripper.Update()
            
            # TODO: Want to convert this to an algorithm to make the whole process a pipeline
            faceMaker = PolyLinesToCapConverter()
            faceMaker.SetInputConnection(stripper.GetOutputPort())
#            face = vtkPolyData()
#            face.SetPoints(stripper.GetOutput().GetPoints())
#            face.SetPolys(stripper.GetOutput().GetLines())
            
            # Join the two together
            joiner = vtkAppendPolyData()
            joiner.AddInputConnection(filter.GetOutputPort())
            joiner.AddInputConnection(faceMaker.GetOutputPort())
#            filter.Update()
#            joiner.AddInput(filter.GetOutput())
#            joiner.AddInput(face)
            
            # Set the source to the capped surface
            pdSource = joiner
            continue
        
        normer = vtkPolyDataNormals()
        normer.SetInputConnection(pdSource.GetOutputPort())
        normer.ComputeCellNormalsOn()
        normer.ComputePointNormalsOff()

        return normer

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

class Writer(object):
    """Generates the HemeLB input file.
    """
    def __init__(self, OutputConfigFile=None, StressType=None, BlockSize=8, BlockCounts=None):
        self.OutputConfigFile = OutputConfigFile
        self.StressType = StressType
        self.BlockSize = BlockSize
        self.BlockCounts = BlockCounts

        # Truncate, and open for read & write in binary mode
        self.file = file(OutputConfigFile,'w+b')
                
        encoder = xdrlib.Packer()
        # Write the preamble, starting with the stress type
        encoder.pack_int(StressType)
        
        # Blocks in each dimension
        for count in BlockCounts:
            encoder.pack_uint(int(count))
            continue
        # Sites along 1 dimension of a block
        encoder.pack_uint(BlockSize)
        
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
    
    def RewriteHeader(self):
        head = self.HeaderEncoder.get_buffer()
        assert len(head) == (self.bodyStart - self.headerStart)
        
        self.file.seek(self.headerStart)
        self.file.write(head)
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
        bb.shape = (3,2)
        
        origin = []
        blocks = []
        for min, max in bb:
            size = max-min
            # int() truncates, we add 2 to make sure there's enough
            # room for the sites just outside.
            nSites = int(size/VoxelSize) + 2
            
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
        self.BlockCounts= np.array(blocks)
        
#        for ijk, val in np.ndenumerate(self.blocks):
#            self.blocks[ijk] = MacroBlock(ijk=ijk, size=BlockSize)
#            continue
        
        return
    
    def CalcPositionFromIndex(self, index):
        return self.Origin + self.VoxelSize * np.array(index)
    
    def GetBlock(self, *blockIjk):
        val = self.blocks[blockIjk]
        if val is None:
            val = self.blocks[blockIjk] = MacroBlock(ijk=blockIjk, size=self.BlockSize)
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
            for i in range(ijk[0]-1, ijk[0]+1):
                if i < 0: continue
                if i == ijk[0] and i != maxInds[0]: continue
                for j in range(ijk[1]-1, ijk[1]+1):
                    if j < 0: continue
                    if j == ijk[1] and j != maxInds[1]: continue
                    for k in range(ijk[2]-1, ijk[2]+1):
                        if k < 0: continue
                        if k == ijk[2] and k != maxInds[2]: continue
                        self.blocks[i,j,k] = None
                        continue
                    continue
                continue
            
                        
            continue # nd.enumerate(self.blocks)
        
        return
    pass

class MacroBlock(object):
    def __init__(self, domain=None, ijk=(0,0,0), size=8):
        self.domain = domain
        self.size = size
        self.ijk = (ijk)
        self.sites = np.empty((size,size,size), dtype=object)
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
        if any(map(lambda x: x<0 or x>=self.size, localSiteIjk)):
            return self.domain.GetSite(*globalSiteIjk)
        
        return self.GetSite(*localSiteIjk)
    
    def CalcPositionFromIndex(self, index):
        return self.domain.CalcPositionFromIndex(index)
    pass

class LatticeSite(object):
    # This ordering of the lattice vectors must be the same as in the HemeLB source.
    neighbours = np.array([[ 1, 0, 0],
                           [-1, 0, 0],
                           [ 0, 1, 0],
                           [ 0,-1, 0],
                           [ 0, 0, 1],
                           [ 0, 0,-1],
                           [ 1, 1, 1],
                           [-1,-1,-1],
                           [ 1, 1,-1],
                           [-1,-1, 1],
                           [ 1,-1, 1],
                           [-1, 1,-1],
                           [ 1,-1,-1],
                           [-1, 1, 1]])
    
    def __init__(self, block=None, ijk=(0,0,0)):
        self.block = block
        self.ijk = ijk
        self.Position = block.CalcPositionFromIndex(ijk)
        
        self.IsFluid = None
        self.IsEdge = None
        self.IsInlet = None
        self.IsOutlet = None
        
        # Attributes that will be updated by Profile.ClassifySite
#        self.Type = None
        self.Normal = None
        self.BoundaryDistance = None
        self.CutDistances = -1. * np.ones(len(self.neighbours))
        
        return
    @property
    def Type(self):
        if self.IsFluid:
            if self.IsEdge:
                return EDGE
            
            return FLUID
        
        if self.IsFluid is None:
            return None
        
        return SOLID
    def IterNeighbours(self):
        """Create an iterator for our neighbours.
        """
        for ijk in self.IterNeighbourIndices:
            yield self.block.GetSite(ijk)
            continue
        return
    
    def IterNeighbourIndices(self):
        """Create an iterator for the indices of our neighbours.
        """
        domain = self.block.domain
        shape = domain.BlockCounts * domain.BlockSize
        for neigh in self.neighbours:
            ind = self.ijk + neigh
            # Skip any out of range
            if np.any(ind<0) or np.any(ind>=shape):
                continue
            
            yield tuple(ind)
            continue
        return
    
    pass

# Site types
SOLID = 0
FLUID = 1
