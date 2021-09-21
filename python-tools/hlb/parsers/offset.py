# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from pathlib import Path
import xdrlib
import numpy as np

from . import HemeLbMagicNumber

OffsetMagicNumber = 0x6F666604
OffsetVersionNumber = 1

# Header has (as uint32):
# - HemeLbMagicNumber
# - OffsetMagicNumber
# - OffsetVersionNumber
# - Number of ranks
OffsetHeaderLength = 4 * 4


class OffsetFile:
    def __init__(self, path):
        if isinstance(path, str):
            path = Path(path)
        self.path = path
        data = self.path.read_bytes()
        self.ReadHeader(data[:OffsetHeaderLength])
        self.ReadData(data[OffsetHeaderLength:])

    def ReadHeader(self, header):
        unp = xdrlib.Unpacker(header)
        assert unp.unpack_uint() == HemeLbMagicNumber
        assert unp.unpack_uint() == OffsetMagicNumber
        assert unp.unpack_uint() == OffsetVersionNumber
        self.NumberOfRanks = unp.unpack_uint()

    def ReadData(self, data):
        unp = xdrlib.Unpacker(data)
        self.Data = np.zeros(self.NumberOfRanks + 1, dtype=np.uint64)
        last = 0
        for i in range(self.NumberOfRanks + 1):
            self.Data[i] = unp.unpack_uhyper()
            assert self.Data[i] >= last
            last = self.Data[i]
