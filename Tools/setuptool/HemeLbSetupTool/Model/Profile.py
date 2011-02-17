import os.path
import pickle
from copy import copy
import numpy as np
import xdrlib
from contextlib import contextmanager

from vtk import vtkSTLReader
from vtk import vtkClipPolyData, vtkAppendPolyData, vtkPlane, vtkStripper, vtkPolyData, \
    vtkFeatureEdges, vtkPolyDataConnectivityFilter, vtkSelectEnclosedPoints

from HemeLbSetupTool.Util.Observer import Observable, ObservableList
from HemeLbSetupTool.Model.SideLengthCalculator import AverageSideLengthCalculator
from HemeLbSetupTool.Model.Vector import Vector
import pdb
class Profile(Observable):
    """This class represents the parameters necessary to perform a
    setup for HemeLb and supplies the functionality to do it.

    The required parameters are below with defaults and all must be
    specified to actually create the setup files.
    
    """
    # Required parameters and defaults.
    _Args = {'StlFile': None,
             'Iolets': ObservableList(),
             'VoxelSize': None,
             'SeedPoint': Vector(),
             'OutputConfigFile': None,
             'OutputXmlFile': None,
             'StressType': 1}
    
    def __init__(self, **kwargs):
        """Required arguments may be set here through keyword arguments.
        """
        # Set attributes on this instance according to the keyword
        # args given here or the default dict if they aren't present.
        for a, default in Profile._Args.iteritems():
            setattr(self, a,
                    kwargs.pop(a, copy(default)))
            continue
        # Raise an error on a kwarg we don't understand
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)

        # We need a reader to get the polydata
        self.StlReader = vtkSTLReader()
        # And a way to estimate the voxel size
        self.sider = AverageSideLengthCalculator()
        self.sider.SetInputConnection(self.StlReader.GetOutputPort())

        # When the STL changes, we should reset the voxel size and
        # update the vtkSTLReader.
        self.AddObserver('StlFile', self.OnStlFileChanged)
        
        # Dependencies for properties
        self.AddDependency('HaveValidStlFile', 'StlFile')
        self.AddDependency('HaveValidOutputXmlFile', 'OutputXmlFile')
        self.AddDependency('HaveValidOutputConfigFile', 'OutputConfigFile')
        self.AddDependency('HaveValidSeedPoint', 'SeedPoint.x')
        self.AddDependency('HaveValidSeedPoint', 'SeedPoint.y')
        self.AddDependency('HaveValidSeedPoint', 'SeedPoint.z')
        self.AddDependency('IsReadyToGenerate', 'HaveValidStlFile')
        self.AddDependency('IsReadyToGenerate', 'HaveValidOutputXmlFile')
        self.AddDependency('IsReadyToGenerate', 'HaveValidOutputConfigFile')
        self.AddDependency('IsReadyToGenerate', 'HaveValidSeedPoint')
        return
    
    def OnStlFileChanged(self, change):
        self.StlReader.SetFileName(self.StlFile)
        self.VoxelSize = self.sider.GetOutputValue()
        return
    
    @property
    def HaveValidStlFile(self):
        """Read only property indicating if our STL file is valid.
        """
        return IsFileValid(self.StlFile, ext='.stl', exists=True)

    @property
    def HaveValidSeedPoint(self):
        if np.isfinite(self.SeedPoint.x) and np.isfinite(self.SeedPoint.y) and np.isfinite(self.SeedPoint.z):
            return True
        return False
    
    @property
    def HaveValidOutputXmlFile(self):
        return IsFileValid(self.OutputXmlFile, ext='.xml')
    @property
    def HaveValidOutputConfigFile(self):
        return IsFileValid(self.OutputConfigFile, ext='.dat')
    
    @property
    def IsReadyToGenerate(self):
        """Read only property indicating if we have enough information
        to do the setup.
        """
        if not self.HaveValidSeedPoint:
            return False
        if not self.HaveValidOutputXmlFile:
            return False
        if not self.HaveValidOutputConfigFile:
            return False
        if not self.HaveValidStlFile:
            return False
        return True
    
    def LoadFromFile(self, filename):
        restored = pickle.Unpickler(file(filename)).load()
        self.CloneFrom(restored)
        return
    
    def Save(self, filename):
        outfile = file(filename, 'w')
        pickler = pickle.Pickler(outfile)
        pickler.dump(self)
        return
    
    def Generate(self):
        """Create the output based on our configuration.
        """
        capper = self._Capper = self.ClipAndCapSurface()
        capper.Update()
        surface = capper.GetOutput()
        
        domain = Domain(self.VoxelSize, surface.GetBounds())
        
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
        try:
            checker = self.IsEnclosedChecker
        except AttributeError:
            checker = self.IsEnclosedChecker = vtkSelectEnclosedPoints()
            checker.Initialize(self._Capper.GetOutput())
            pass
        
        return checker.IsInsideSurface(site.Position)
    
    def ClassifySite(self, site):
        
        if self.IsSiteEnclosed(site):
            site.Type = FLUID
        else:
            site.Type = SOLID
            pass
        
        return
    
    def ClipAndCapSurface(self):
        """Clip the PolyData read from the STL file against the planes.
        """
        # Have the name pdSource first point to the input, then loop
        # over IOlets, clipping and capping.
        pdSource = self.StlReader
        
        for iolet in self.Iolets:
#            pdb.set_trace()
            plane = vtkPlane()
            plane.SetOrigin(iolet.Centre.x, iolet.Centre.y, iolet.Centre.z)
            plane.SetNormal(iolet.Normal.x, iolet.Normal.y, iolet.Normal.z)
            # Clip: gives the surface we want
            clipper = vtkClipPolyData()
            clipper.SetClipFunction(plane)
            clipper.SetInputConnection(pdSource.GetOutputPort())
            
            filter = vtkPolyDataConnectivityFilter()
            filter.SetInputConnection(clipper.GetOutputPort())
            filter.SetClosestPoint(self.SeedPoint.x, self.SeedPoint.y, self.SeedPoint.z)
            filter.SetExtractionModeToClosestPointRegion()
            
            # Cut: generates the cap
#            cutter = vtkCutter()
#            cutter.SetCutFunction(plane)
#            cutter.SetInputConnection(pdSource.GetOutputPort())
            
            edger = vtkFeatureEdges()
            edger.SetInputConnection(filter.GetOutputPort())
            edger.BoundaryEdgesOn()
            edger.FeatureEdgesOff()
            edger.NonManifoldEdgesOff()
            edger.ManifoldEdgesOff()
            
            stripper = vtkStripper()
            stripper.SetInputConnection(edger.GetOutputPort())
            stripper.Update()
            
            # Want to convert this to an algorithm to make the whole process a pipeline
            face = vtkPolyData()
            face.SetPoints(stripper.GetOutput().GetPoints())
            face.SetPolys(stripper.GetOutput().GetLines())
            
            # Join the two together
            joiner = vtkAppendPolyData()
#            joiner.AddInputConnection(clipper.GetOutputPort())
#            joiner.AddInputConnection( cutter.GetOutputPort())
            filter.Update()
            joiner.AddInput(filter.GetOutput())
            joiner.AddInput(face)
            
            # Set the source to the capped surface
            pdSource = joiner
            continue
        
        return pdSource
    
    def ResetVoxelSize(self, ignored=None):
        """Action to reset the voxel size to its default value.
        """
        self.VoxelSize = self.sider.GetOutputValue()
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
    def __init__(self, block=None, ijk=(0,0,0)):
        self.block = block
        self.ijk = ijk
        self.Position = block.CalcPositionFromIndex(ijk)
        
        # Attributes that will be updated by Profile.ClassifySite
        self.Type = None
        self.Normal = None
        self.BoundaryDistance = None
        self.CutDistances = None
        
        return
    
    pass

# Site types
SOLID = 0
FLUID = 1

def IsFileValid(path, ext=None, exists=None):
    if not isinstance(path, (str, unicode)):
        return False
    if path == '':
        return False
    
    if exists is not None:
        if os.path.exists(path) != exists:
            return False
        pass
    
    if ext is not None:
        ending = os.path.splitext(path)[1]
        if ending != ext:
            return False
        pass
    return True
