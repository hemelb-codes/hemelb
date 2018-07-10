# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import numpy as np

from hemeTools.parsers.geometry.simple import ConfigLoader
from hemeTools.parsers.geometry.generic import Domain

from hemeTools.parsers.octree import SectionTree

class HeaderEndException(Exception):
    pass

class CountingLoader(ConfigLoader):
    def __init__(self, filename):
        self.GmyFileName = filename
        self.Domain = Domain()
        self.Domain.VoxelSize = 1.0
        self.Domain.Origin = np.zeros(3, dtype=float)
        self.File = file(self.GmyFileName)

    def OnEndHeader(self):
        raise HeaderEndException()
    pass

def CountFluidSitesGmy(filename, verbosity=1):
    summary = 'FileName = "{filename}"\n'
    for var in 'BlockCounts BlockSize TotalFluidSites BlocksWithFluidSites'.split():
        summary += '\t%s = {%s}\n' % (var,var)
        continue
    
    ldr = CountingLoader(filename)
    try:
        ldr.Load()
    except HeaderEndException:
        # OK
        pass
    else:
        # Error
        raise RunTimeError("The expected exception was not raised!")
    
    BlockCounts = ldr.Domain.BlockCounts
    BlockSize = ldr.Domain.BlockSize
    TotalBlocks = np.prod(BlockCounts)
    
    FluidSitesPerBlock = ldr.Domain.BlockFluidSiteCounts
    BytesPerBlock = ldr.BlockDataLength
    UncompressedBytesPerBlock = ldr.BlockUncompressedDataLength

    TotalFluidSites = FluidSitesPerBlock.sum()
    BlocksWithFluidSites = np.sum(FluidSitesPerBlock>0)
    if verbosity == 0:
        return TotalFluidSites
    elif verbosity == 1:
        return summary.format(**locals())

    # verbosity >= 2
    summary = summary.format(**locals())
    summary += '\t\tIndex\tSites\tZipped\tUnzipped\n'
    for i in xrange(TotalBlocks):
        summary += '\t\t%d\t%d\t%d\t%d\n' % (i, FluidSitesPerBlock[i], BytesPerBlock[i], UncompressedBytesPerBlock[i])
        continue
    
    return summary

def CountFluidSitesOct(filename, verbosity=1):
    tree = SectionTree(filename)
    
    def valid_child_count(level):
        ds = tree.indices[level]
        last_8 = tree.indices[level][-8:]
        valid = (last_8 != SectionTree.NA)
        last_i = int(last_8[valid].max())
        if level > 1:
            last_i /= 8
        return last_i + 1
    
    tot_fluid = valid_child_count(1)
    
    if verbosity == 0:
        return tot_fluid
    
    nlevels = int(tree.levels)
    boxsize = 2**nlevels
    
    summary = '''FileName = "{filename}"
\tLevels = {nlevels}
\tBoxSize = {boxsize}
\tTotalFluidSites = {tot_fluid}
'''.format(**locals())

    if verbosity == 1:
        return summary
    
    # v >= 2
    summary += '\t\tLevel\tNodes\tFrac\tCumulative fraction\n'
    t = '\t\t{: 2d}\t{: %dd}\t{:.2f}\t{:.2e}\n' % (int(np.ceil(np.log10(tot_fluid)))+1)
    
    # Root node always there
    summary += t.format(nlevels, 1, 1.0, 1.0)
    last = 1
    cum_frac = 1.0
    for lvl in xrange(nlevels-1, 0, -1):
        x = valid_child_count(lvl+1)
        frac = x/(last*8.0)
        cum_frac *= frac
        summary += t.format(lvl, x, frac, cum_frac)
        last = x
        
    frac = tot_fluid/(8.0*last)
    cum_frac *= frac
    summary += t.format(0, tot_fluid, frac, cum_frac)
    return summary

def get_magic(fn):
    '''Gets the magic number from the file.'''
    with open(fn) as f:
        return f.read(8)

magic = {
    'hlb!gmy\x04': CountFluidSitesGmy,
    '\x89HDF\r\n\x1a\n': CountFluidSitesOct
}

def CountFluidSites(filename, verbosity=1):
    return magic[get_magic(filename)](filename, verbosity)

def main():
    import sys
    verbosity = 0
    inputs = sys.argv[1:]
    while inputs[0] == '-v':
        inputs.pop(0)
        verbosity += 1
    
    for fn in inputs:
        print CountFluidSites(fn, verbosity)
            
