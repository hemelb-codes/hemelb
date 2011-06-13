#!/usr/bin/env python

import xdrlib
import numpy as np
import pdb

def assert_equal(var1, var2, msg):
    if var1 != var2:
        print 'ERROR: ' + msg + ' (var1 had ' + str(var1) + ', var2 had ' + str(var2) + ')'

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

        print '-----\nInfo: stress type = ' + str(self.StressType) + ', sites per block = ' + str(self.BlockSize ** 3) + \
          ', blocks: ' + str(self.BlockCounts[0]) + 'x' + str(self.BlockCounts[1]) + 'x' + str(self.BlockCounts[2]) +  \
          ', vox size = ' + str(self.VoxelSize) + ', origin = ' + str(self.Origin) + '\n-----\n'

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

        print '-----\nFluid Site Count = ' + str(np.sum(self.BlockFluidSiteCounts)) + '\n-----\0n'
    
        # CONSTRAINT 1: the length of the file must be equal to the value we calculate from the headers.
        assert_equal(headerEnd + np.sum(self.BlockDataLength), len(data), 'File length appears incorrect.')
        
        self.Blocks = np.zeros(self.TotalBlocks, dtype=object)

        i = [0,0,0]
        ijk = 0
        for i[0] in xrange(self.BlockCounts[0]):
            for i[1] in xrange(self.BlockCounts[1]):
                for i[2] in xrange(self.BlockCounts[2]):
                    if ijk % BLOCK_REPORT_FREQ == 0:
                        print 'Reading block ' + str(ijk) +  ' of ' + str(self.TotalBlocks)
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
        elif not self.Edge:
            assert_equal(false, true, 'Site doesn\'t appear to have any fitting type.')    

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

BLOCK_REPORT_FREQ = 200

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

    dom = Domain(sys.argv[1])

    # Neighbour deltas are the vectors which when added to a fluid site position vector give the location of a putative neighbour.
    neighs = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1], \
              [1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,1,-1],[-1,-1,1],[-1,-1,-1]]

    domainExtremeSite = dom.BlockCounts * dom.BlockSize - 1

    # Now we iterate over every site on every block, and check the following constraints.
    # CONSTRAINT 2: the number of non-solid sites on each block must equal that declared in the header.
    # CONSTRAINT 3: every site which has a fluid type (and isn't an edge) can have no invalid or solid-type neighbours.
    # CONSTRAINT 4: every site which is not a fluid type must have at least one invalid or solid-type neighbour.
    # CONSTRAINT 5: if the site is not of fluid or solid type, at least one of the cut distances must be between 0 and 1
    # CONSTRAINT 6: if the site is not of fluid or solid type, all cut distances must be between 0 and 1 or be infinite
    # CONSTRAINT 7: if the site is of inlet or outlet type, the magnitude of the boundary normal must be 1
    # CONSTRAINT 8: if the site is of inlet or outlet type, the boundary distance must be between 0 and 1
    # CONSTRAINT 9: if the site is an edge, the magnitude of the wall normal must be 1
    # CONSTRAINT 10: if the site is an edge, the distance to the wall must be between 0 and root(3)

    i = np.zeros(3, dtype=int)
    ijk = 0
    for i[0] in xrange(dom.BlockCounts[0]):
        for i[1] in xrange(dom.BlockCounts[1]):
            for i[2] in xrange(dom.BlockCounts[2]):

                if ijk % BLOCK_REPORT_FREQ == 0:
                    print 'Checking block: ' + str(ijk) + ' of ' + str(dom.TotalBlocks)

                numFluid = 0
                block = dom.Blocks[i[0],i[1],i[2]]

                j = np.zeros(3, dtype=int)

                minInd = i*dom.BlockSize
                maxInd = (i+1)*dom.BlockSize

                for j[0] in xrange(minInd[0], maxInd[0], 1):
                    for j[1] in xrange(minInd[1], maxInd[1], 1):
                        for j[2] in xrange(minInd[2], maxInd[2], 1):

                            site = dom.Sites[j[0],j[1],j[2]]
                            type = site.Type
                            edge = site.Edge

                            if type == SOLID_TYPE:
                               continue
                            else:
                                numFluid += 1

                            nonFluidNeighs = 0
                            offender = [0,0,0]
                            invalid = False
                            for delta in neighs:
                                neigh = j + delta
                                if min(neigh) < 0 or min(domainExtremeSite - neigh) < 0:
                                    nonFluidNeighs += 1
                                    offender = neigh
                                    invalid = True
                                    break

                                elif dom.Sites[neigh[0],neigh[1],neigh[2]].Type == SOLID_TYPE:
                                    nonFluidNeighs += 1
                                    offender = neigh
                                    break

                            # CON 3
                            if type == FLUID_TYPE and not edge:
                                assert_equal(nonFluidNeighs, 0, 'Fluid site at ' + str(j) + ' has type, edge of ' + str(type) + ', ' + str(edge) + ' and has some neighbours that aren\'t valid or are solids, e.g. site at ' + str(offender) + (' which is out-of-bounds' if invalid else ('with type ' + str(dom.Sites[offender[0], offender[1], offender[2]].Type))))

                            # CON 4
                            if type != FLUID_TYPE or edge:
                                assert_equal(nonFluidNeighs > 0, True, 'Site is non-normal fluid but has all valid fluid neighbours.')
                                
                            if type != FLUID_TYPE:
                                # CON 5
                                assert_equal(min(site.CutDistances) < 1., True, 'Non-fluid site had no cut distances between 0 and 1')
                                for cutDist in site.CutDistances:
                                    if cutDist < 0. or (cutDist > 1 and cutDist != float('inf')):
                                        # CON 6
                                        assert_equal(False, True, 'Invalid cut distance of ' + str(cutDist))

                            if type == INLET_TYPE or type == OUTLET_TYPE:
                                normMag = np.sum(site.BoundaryNormal ** 2) ** 0.5
                                # CON 7
                                assert_equal( abs(normMag - 1) < 0.001, True, 'Boundary normal had non-unity magnitude: ' + str(site.BoundaryNormal))
                                # CON 8
                                assert_equal(site.BoundaryDistance > 0. and site.BoundaryDistance < 1., True, 'Boundary distance was not in [0,1]: ' \
                                    + str(site.BoundaryDistance))

                            if edge:
                                normMag = np.sum(site.WallNormal ** 2) ** 0.5
                                # CON 9
                                assert_equal( abs(normMag - 1) < 0.001, True, 'Wall normal had non-unity magnitude: ' + str(site.WallNormal))
                                # CON 10
                                assert_equal(site.WallDistance > 0. and site.WallDistance < (3.0**0.5), True, 'Wall distance was not in [0,root(3)]: ' \
                                    + str(site.WallDistance))

                # CON 2
                assert_equal(numFluid, block.nFluidSites, 'Fluid site counts not equal.')
                ijk += 1
