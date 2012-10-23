# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import pdb
import xdrlib
import numpy as N

class Image(object):
    
    def __init__(self, filename):
        self.filename = filename
        
        f = file(filename)
        reader = xdrlib.Unpacker(f.read())
        self.mode = reader.unpack_uint()
        self.pressure_threshold = (reader.unpack_float(), reader.unpack_float())
        self.velocity_max = reader.unpack_float()
        self.stress_max = reader.unpack_float()
        self.screen = (reader.unpack_uint(), reader.unpack_uint())
        self.nPixels = reader.unpack_uint()

        tempPixels = N.zeros((self.nPixels,4), dtype=N.uint32)

        for i in xrange(self.nPixels):
            for j in xrange(4):
                tempPixels[i,j] = reader.unpack_uint()
                continue
            continue
        sortedIndices = N.argsort(tempPixels[:,0])
        self.pixels = tempPixels[sortedIndices].view(dtype=[('index', N.uint16, 2),
                                                            ('r', N.uint8, 4),
                                                            ('g', N.uint8, 4),
                                                            ('b', N.uint8, 4)]).view(N.recarray)

        # The above statement puts the pixel indices in the wrong order (i.e. j,i rather than i,j).
        # This is corrected here.
        self.pixels.index[:,0] = self.pixels.index[:, 0, [1,0]]

        return
    
    def old__init__(self, filename):
        self.filename = filename
        
        f = file(filename)
        reader = xdrlib.Unpacker(f.read())
        self.mode = reader.unpack_uint()
        self.pressure_threshold = (reader.unpack_float(), reader.unpack_float())
        self.velocity_max = reader.unpack_float()
        self.stress_max = reader.unpack_float()
        self.screen = (reader.unpack_uint(), reader.unpack_uint())
        self.nPixels = reader.unpack_uint()

        tempPixels = N.zeros(self.nPixels, dtype=type(self).RowType)

        # Masks and shifts to get the bits for each vis mode
        masks = [(2**32-1) ^ (2**24-1),
                (2**24-1) ^ (2**16-1),
                (2**16-1) ^ (2**8-1),
                (2**8-1)]
        shifts = [24, 16, 8, 0]

        conversions = zip(masks, shifts)
        for i in xrange(self.nPixels):
            ind = reader.unpack_uint()
            tempPixels[i]['index'][0] = (ind & ((2**32-1) ^ (2**16-1))) >> 16
            tempPixels[i]['index'][1] = ind & (2**16-1)
            r = reader.unpack_uint()
            g = reader.unpack_uint()
            b = reader.unpack_uint()

            for j, (mask, shift) in enumerate(conversions):
                tempPixels[i]['r'][j] = (r & mask) >> shift
                tempPixels[i]['g'][j] = (g & mask) >> shift
                tempPixels[i]['b'][j] = (b & mask) >> shift
                continue

            # assert r == (
            #     (int(tempPixels[i]['r'][0]) << 24) +
            #     (int(tempPixels[i]['r'][1]) << 16) +
            #     (int(tempPixels[i]['r'][2]) <<  8) +
            #     (int(tempPixels[i]['r'][3]) <<  0)
            #     )
            # assert g == (
            #     (int(tempPixels[i]['g'][0]) << 24) +
            #     (int(tempPixels[i]['g'][1]) << 16) +
            #     (int(tempPixels[i]['g'][2]) <<  8) +
            #     (int(tempPixels[i]['g'][3]) <<  0)
            #     )
            # assert b == (
            #     (int(tempPixels[i]['b'][0]) << 24) +
            #     (int(tempPixels[i]['b'][1]) << 16) +
            #     (int(tempPixels[i]['b'][2]) <<  8) +
            #     (int(tempPixels[i]['b'][3]) <<  0)
            #     )
            
            continue
        sortedIndices = N.argsort(tempPixels['index'][:, 0]<<16 + tempPixels['index'][:, 1])
        self.pixels = tempPixels[sortedIndices]
        
        return

    def __eq__(self, other):
        """Return True if the Images are identical.
        """
        attrs = ('mode', 'pressure_threshold', 'velocity_max', 'stress_max', 'screen')
        for at in attrs:
            if not getattr(self, at) == getattr(other, at):
                return False
            continue
        
        if not N.alltrue(self.pixels.view(dtype=N.uint32) == other.pixels.view(dtype=N.uint32)):
            return False
        
        return True

    def almost_eq(self, other, tol=1):
        """Same as __eq__, except allow data array to differ by up to
        tol (default 1).
        """
        
        attrs = ('mode', 'pressure_threshold', 'velocity_max', 'stress_max', 'screen')
        for at in attrs:
            if not getattr(self, at) == getattr(other, at):
                print 'differed on ' + at + ': ' + str(getattr(self, at)) + ' and ' + str(getattr(other, at))
                return False
            continue

        if not N.alltrue(self.pixels.index == other.pixels.index):
            for i in range(self.pixels.index.shape[0]):
                if not N.alltrue(self.pixels.index[i] == other.pixels.index[i]):
                    print 'Images differed on indices, first at index: ' + str(i)
                    print 'left array had ' + str(self.pixels.index[i])
                    print 'right array had ' + str(other.pixels.index[i])
                    return False

        for attr in ('r', 'g', 'b'):
            d = self.pixels[attr] - other.pixels[attr] + 128
            if N.min(d) < (128-tol) or N.max(d) > (128+tol):
                minArg = N.argmin(d, 0)[0][0]
                maxArg = N.argmax(d, 0)[0][0]
                print 'Differed on ' + attr + ' channel with delta range:'
                print str(N.min(d) - 128) + ' between ' + str(self.pixels[minArg]) + ' and ' + str(other.pixels[minArg])
                print 'and'
                print str(N.max(d) - 128) + ' between ' + str(self.pixels[maxArg]) + ' and ' + str(other.pixels[maxArg])
                return False
            continue
        
        return True
    
    pass

if __name__ == "__main__":
    import sys
    images = [Image(arg) for arg in sys.argv[1:]]
    
