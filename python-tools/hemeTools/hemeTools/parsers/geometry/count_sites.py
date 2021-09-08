# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from .simple import ConfigLoader
from .generic import Domain
import numpy as np


class HeaderEndException(Exception):
    pass


class CountingLoader(ConfigLoader):
    def __init__(self, filename):
        self.GmyFileName = filename
        self.Domain = Domain()
        self.Domain.VoxelSize = 1.0
        self.Domain.Origin = np.zeros(3, dtype=float)
        self.File = open(self.GmyFileName, "rb")

    def OnEndHeader(self):
        raise HeaderEndException()

    pass


def CountFluidSites(filename, verbosity=1):
    summary = 'FileName = "{filename}"\n'
    for var in "BlockCounts BlockSize TotalFluidSites BlocksWithFluidSites".split():
        summary += "\t%s = {%s}\n" % (var, var)
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
    BlocksWithFluidSites = np.sum(FluidSitesPerBlock > 0)
    if verbosity == 0:
        return TotalFluidSites
    elif verbosity == 1:
        return summary.format(**locals())

    # verbosity >= 2
    summary = summary.format(**locals())
    summary += "\t\tIndex\tSites\tZipped\tUnzipped\n"
    for i in range(TotalBlocks):
        summary += "\t\t%d\t%d\t%d\t%d\n" % (
            i,
            FluidSitesPerBlock[i],
            BytesPerBlock[i],
            UncompressedBytesPerBlock[i],
        )
        continue

    return summary


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--verbose", "-v", action="count", default=0)
    p.add_argument("inputs", nargs="+")
    args = p.parse_args()

    for fn in args.inputs:
        print(CountFluidSites(fn, args.verbose))


if __name__ == "__main__":
    main()
