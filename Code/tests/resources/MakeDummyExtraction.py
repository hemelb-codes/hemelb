
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import xdrlib
import numpy as np

# Set things up
HemeLbMagicNumber = 0x686c6221
ExtractionMagicNumber = 0x78747204
Version = 0x2

voxelSizeMetres = 0.3e-3
originMetres = (0.034, 0.001, 0.074)

# Name, number of floats
fields = (
    ('Pressure', 1),
    ('Velocity', 3),
    )
fieldCount = len(fields)


siteCount = 64
# Positions as for a FourCube
gridPositions = np.mgrid[:4,:4,:4].reshape(3, siteCount).transpose()

############################################################
# Main Header

mainHeaderEncoder = xdrlib.Packer()

mainHeaderEncoder.pack_uint(HemeLbMagicNumber)
mainHeaderEncoder.pack_uint(ExtractionMagicNumber)
mainHeaderEncoder.pack_uint(Version)

mainHeaderEncoder.pack_double(voxelSizeMetres)

for o in originMetres:
    mainHeaderEncoder.pack_double(o)

mainHeaderEncoder.pack_uhyper(siteCount)
mainHeaderEncoder.pack_uint(fieldCount)

# Not quite done with this encoder; need to add the field header's length, so create that first.

############################################################
# Field Header

fieldHeaderEncoder = xdrlib.Packer()

for fieldName, fieldLength in fields:
    fieldHeaderEncoder.pack_string(fieldName)
    fieldHeaderEncoder.pack_uint(fieldLength)

fieldHeaderLength = len(fieldHeaderEncoder.get_buffer())

mainHeaderEncoder.pack_uint(fieldHeaderLength)

############################################################
# Body

bodyEncoder = xdrlib.Packer()

# Need a source of data; use a linear congruential PRNG as no requirement
# for it to be "good". Implentation from Numerical Recipes
class PRNG(object):
    def __init__(self, seed):
        """m = modulus
        a = multiplier
        c = increment
        
        """
        self.m = 2**32
        self.a = np.uint32(1664525)
        self.c = np.uint32(1013904223)

        self.state = np.uint32(seed)
        
        self._norm = np.float32(np.uint32(-1))
        return
    
    def rand(self):
        """Get a random int by
        State[n+1] = (a State[n] + c) mod m
        Since we have m == 2**32 the overflow does the modulus operation for us.
        """
        self.state = (self.a * self.state + self.c)
        return self.state
    
    def uniform(self):
        return self.rand() / self._norm

    pass

# Random number generator with seed of time of writing to produce body data.
prng = PRNG(1358)
    
for timestep in (0, 100):
    # Timestep
    bodyEncoder.pack_uhyper(timestep)
    
    for iSite in xrange(siteCount):
        # Grid pos
        for g in gridPositions[iSite]:
            bodyEncoder.pack_uint(g)
            continue
        
        # Pressure
        pressureMmHg = 80. + 2 * prng.uniform()
        bodyEncoder.pack_float(pressureMmHg)

        # Velocity
        for iIgnored in xrange(3):
            velocityMetresPerSecond = 0.01 * prng.uniform()
            bodyEncoder.pack_float(velocityMetresPerSecond)
            continue
        
        continue
    continue

# Write the file
xtrFile = file('dummy.xtr', 'w')
xtrFile.write(mainHeaderEncoder.get_buffer())
xtrFile.write(fieldHeaderEncoder.get_buffer())
xtrFile.write(bodyEncoder.get_buffer())
xtrFile.close()
