#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
"""Tools to decompress and recompress HemeLB geometry files.

These can be useful for debugging.
"""

import argparse
import xdrlib
import zlib

import numpy as np
from six.moves import range

from .. import HemeLbMagicNumber
from . import GeometryMagicNumber
from .simple import ConfigLoader


class CompressionBase(ConfigLoader):
    """Most of the functionality goes here.

    Subclasses need to override _LoadBlock to actually do the
    (de)compression."""

    def __init__(self, filename, outfilename):
        is_gmy = sniff_gmy(filename)
        ConfigLoader.__init__(self, filename, gmy=is_gmy)
        self.OutputFileName = outfilename

    def OnEndPreamble(self):
        # Copy the preamble
        pos = self.File.tell()
        self.File.seek(0)

        self.OutFile = open(self.OutputFileName, "wb")
        self.OutFile.write(self.File.read(self.PreambleBytes))
        self.File.seek(pos)

    def OnEndHeader(self):
        # Write dummy values
        nBlocks = self.Domain.TotalBlocks
        self.OutFile.write((4 * 3 * nBlocks) * b"\0")

    def OnEndBody(self):
        # Write a real header
        packer = xdrlib.Packer()
        for i in range(self.Domain.TotalBlocks):
            packer.pack_uint(self.Domain.BlockFluidSiteCounts[i])
            packer.pack_uint(self.BlockDataLength[i])
            packer.pack_uint(self.BlockUncompressedDataLength[i])

        self.OutFile.seek(self.PreambleBytes)
        self.OutFile.write(packer.get_buffer())
        self.OutFile.close()


class Decompressor(CompressionBase):
    def _LoadBlock(self, domain, bIdx, bIjk):
        if domain.BlockFluidSiteCounts[bIjk] == 0:
            return
        compressed = self.File.read(self.BlockDataLength[bIjk])
        uncompressed = zlib.decompress(compressed)
        self.BlockDataLength[bIjk] = self.BlockUncompressedDataLength[bIjk]
        self.OutFile.write(uncompressed)
        return


class Compressor(CompressionBase):
    def _LoadBlock(self, domain, bIdx, bIjk):
        if domain.BlockFluidSiteCounts[bIjk] == 0:
            return
        uncompressed = self.File.read(self.BlockDataLength[bIjk])
        compressed = zlib.compress(uncompressed)
        self.BlockUncompressedDataLength[bIjk] = len(uncompressed)
        self.BlockDataLength[bIjk] = len(compressed)
        self.OutFile.write(compressed)
        return


def sniff_gmy(filename):
    """Look for HemeLB Geometry magic numbers at the start of the file."""
    with open(filename, "rb") as f:
        data = f.read(8)
    up = xdrlib.Unpacker(data)
    hlbNumber = up.unpack_uint()
    gmyNumber = up.unpack_uint()
    return hlbNumber == HemeLbMagicNumber and gmyNumber == GeometryMagicNumber


def mk_argparser():
    argp = argparse.ArgumentParser()
    argp.add_argument("input")
    argp.add_argument("output")
    return argp


def compress_main():
    args = mk_argparser().parse_args()
    comp = Compressor(args.input, args.output)
    comp.Load()


def decompress_main():
    args = mk_argparser().parse_args()
    decomp = Decompressor(args.input, args.output)
    decomp.Load()
