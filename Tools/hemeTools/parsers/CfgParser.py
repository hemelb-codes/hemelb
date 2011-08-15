#!/usr/bin/env python

import xdrlib
import numpy as np
import pdb

class Domain(object):
    def __init__(self, cfgFileName):
        data = file(cfgFileName).read()
        headerStart = 5*4 + 4*8
        preambleLoader = xdrlib.Unpacker(data[0:headerStart])

        self.StressType = preambleLoader.unpack_uint()
        self.BlockCounts = np.array([preambleLoader.unpack_uint() for i in xrange(3)])
        self.BlockSize = preambleLoader.unpack_uint()

        self.VoxelSize = preambleLoader.unpack_double()
        self.Origin = np.array([preambleLoader.unpack_double() for i in xrange(3)])

        self.TotalBlocks = np.prod(self.BlockCounts)
        
        headerEnd = headerStart + self.TotalBlocks * 2 * 4

        headerLoader = xdrlib.Unpacker(data[headerStart:headerEnd])
        self.BlockFluidSiteCounts = np.zeros(self.TotalBlocks, dtype=np.uint)
        self.BlockDataLength = np.zeros(self.TotalBlocks, dtype=np.uint)
        blockStarts = np.zeros(self.TotalBlocks, dtype=np.uint)
        blockEnds = np.zeros(self.TotalBlocks, dtype=np.uint)
        
        i = [0,0,0]
        ijk = 0
        
        blockEnds[-1] = headerEnd
        
        for i[0] in xrange(self.BlockCounts[0]):
            for i[1] in xrange(self.BlockCounts[1]):
                for i[2] in xrange(self.BlockCounts[2]):
                    self.BlockFluidSiteCounts[ijk] = headerLoader.unpack_uint()
                    self.BlockDataLength[ijk] = headerLoader.unpack_uint()
                    blockStarts[ijk] = blockEnds[ijk-1]
                    blockEnds[ijk] = blockStarts[ijk] + self.BlockDataLength[ijk]
                    ijk += 1
        
        self.Blocks = np.zeros(self.TotalBlocks, dtype=object)

        i = [0,0,0]
        ijk = 0
        for i[0] in xrange(self.BlockCounts[0]):
            for i[1] in xrange(self.BlockCounts[1]):
                for i[2] in xrange(self.BlockCounts[2]):
                    blockLoader = xdrlib.Unpacker(data[blockStarts[ijk]:blockEnds[ijk]])
                    
                    self.Blocks[ijk] = Block(self, np.array(i), self.BlockFluidSiteCounts[ijk], blockLoader)
                    ijk += 1
        
        self.Blocks.shape = self.BlockCounts
        self.Sites = np.zeros(self.BlockCounts * self.BlockSize, dtype=object)
        i = np.zeros(3, dtype=int)
        ijk = 0
        for i[0] in xrange(self.BlockCounts[0]):
            for i[1] in xrange(self.BlockCounts[1]):
                for i[2] in xrange(self.BlockCounts[2]):
                    minInd = i*self.BlockSize
                    maxInd = (i+1)*self.BlockSize
                    self.Sites[minInd[0]:maxInd[0],
                               minInd[1]:maxInd[1],
                               minInd[2]:maxInd[2]] = self.Blocks[i[0], i[1], i[2]].Sites
                    ijk += 1
                    
        return
    pass

class Block(object):
    def __init__(self, domain, bIndex, nFluid, loader):
        self.Domain = domain
        self.Index = bIndex
        self.nFluidSites = nFluid
        self.Sites = np.zeros(domain.BlockSize**3, dtype=object)
        
        allSolid = nFluid == 0
        
        i = [0,0,0]
        ijk = 0
        for i[0] in xrange(domain.BlockSize):
            for i[1] in xrange(domain.BlockSize):
                for i[2] in xrange(domain.BlockSize):
                    self.Sites[ijk] = Site(self, domain.BlockSize * bIndex + i, loader, allSolid)
                    ijk += 1
                    
        self.Sites.shape = (domain.BlockSize, domain.BlockSize, domain.BlockSize)
        return
    pass

class Site(object):
    def __init__(self, block, sIndex, loader, allSolid=False):
        self.Block = block
        self.Index = sIndex
        if allSolid:
            self.Config = SOLID_TYPE
            return
        
        self.Config = loader.unpack_uint()
        self.Position = block.Domain.Origin + block.Domain.VoxelSize * sIndex
        if self.Config == SOLID_TYPE:
            return
        if self.Config == FLUID_TYPE:
            return
        
        if self.Type == INLET_TYPE or self.Type == OUTLET_TYPE:
            self.BoundaryNormal = np.array([loader.unpack_double() for i in xrange(3)])
            self.BoundaryDistance = loader.unpack_double()

        if self.Edge:
            self.WallNormal = np.array([loader.unpack_double() for i in xrange(3)])
            self.WallDistance = loader.unpack_double()
        self.CutDistances = np.array([loader.unpack_double() for i in xrange(14)])
        
        return

    @property
    def Type(self):
        return GetType(self.Config)
    @property
    def BoundaryConfig(self):
        return GetType(self.Config)
    @property
    def BoundaryId(self):
        return GetBoundaryId(self.Config)
    @property
    def Edge(self):
        return bool(GetPressureEdge(self.Config))
    
    pass

def GetType(cfg):
    return cfg & SITE_TYPE_MASK
def GetBoundaryConfig(cfg):
    return (cfg & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT
def GetBoundaryDirMask(cfg):
    return (cfg & BOUNDARY_DIR_MASK) >> BOUNDARY_DIR_SHIFT
def GetBoundaryIdMask(cfg):
    return (cfg & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT
def GetPressureEdge(cfg):
    return cfg & PRESSURE_EDGE_MASK

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

BOUNDARY_CONFIG_SHIFT = SITE_TYPE_BITS
BOUNDARY_DIR_SHIFT = BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS
BOUNDARY_ID_SHIFT = BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS

# ===============================================================================
# Horrifying bit-fiddling masks courtesy of Marco.
# Comments show the bit patterns.
# ===============================================================================
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

PRESSURE_EDGE_MASK = 1
PRESSURE_EDGE_MASK <<= 31
# 1000 0000  0000 0000  0000 0000  0000 0000

def GetSiteAttr(domain, attr, otype=int):
    getter = np.vectorize(lambda s: getattr(s, attr), otypes=[otype])
    return getter(domain.Sites)

if __name__ == "__main__":
    import sys
    from enthought.tvtk.api import tvtk as vtk

    dom = Domain(sys.argv[1])
    #types = GetSiteAttr(dom, 'Type').astype(int)
    # pos = GetSiteAttr(dom, 'Position')

    #points = vtk.StructuredPoints(dimensions=types.shape,
    #                              origin=dom.Origin,
    #                              spacing=3*[dom.VoxelSize,])
    #points.point_data.scalars = types.transpose().flatten()

    #writer = vtk.XMLDataSetWriter(file_name='types.vti')
    #writer.input = points
    #writer.write()

    i=0
    j=0

    nSites=0
    nFluid=0
    nSolid=0
    nInlet=0
    nOutlet=0

    for i in xrange(dom.TotalBlocks):
        nSites += dom.BlockFluidSiteCounts[i]
        
    i = [0,0,0]
    for i[0] in xrange(dom.BlockCounts[0] * dom.BlockSize):
        for i[1] in xrange(dom.BlockCounts[1] * dom.BlockSize):
            for i[2] in xrange(dom.BlockCounts[2] * dom.BlockSize):
                if dom.Sites[i[0], i[1], i[2]].Type == FLUID_TYPE:
                    nFluid += 1
                elif dom.Sites[i[0], i[1], i[2]].Type == SOLID_TYPE:
                    nSolid += 1
                elif dom.Sites[i[0], i[1], i[2]].Type == INLET_TYPE:
                    nInlet += 1
                elif dom.Sites[i[0], i[1], i[2]].Type == OUTLET_TYPE:
                    nOutlet += 1

    print "Total sites: " + str(nSites) 
    print "Fluid: " + str(nFluid)
    print "Solid: " + str(nSolid)
    print "Inlet: " + str(nInlet)
    print "Outlet: " + str(nOutlet) 
    
