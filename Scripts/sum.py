#! /usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# Create a new checkpoint file from two existing files by performing an
# element by element sum.

import collections
import xdrlib
import sys
import os

# ******************************************************************************
# Define some constants for later use.
# ******************************************************************************

HemeLbMagicNumber = 0x686c6221
extractionMagicNumber = 0x78747204
extractionVersionNumber = 4

sizeOfUint32 = 4
sizeOfUint64 = 8
sizeOfFloat = 4
sizeOfDouble = 8

# ******************************************************************************
# Define a function to read the preamble, header and timestep from a checkpoint
# file.
# A checkpoint file is an extraction file..
# ******************************************************************************

FileInfo = collections.namedtuple('FileInfo', ['voxelSize', 'origin', 'numSites', 'numberOfComponents', 'timestep'])

def readCheckpointFilePreambleEtc(checkpointFile, firstOffset, lastOffset):
    # Read the preamble of the checkpoint file, i.e.
    # read two magic numbers, a version number, the voxel size,
    # the origin of the grid and the number of sites.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    preamble = checkpointFile.read(preambleBytes)
    up1 = xdrlib.Unpacker(preamble)
    firstMagicNumber = up1.unpack_uint()
    assert firstMagicNumber == HemeLbMagicNumber, "firstMagicNumber (%r) is not equal to the HemeLbMagicNumber" % firstMagicNumber
    secondMagicNumber = up1.unpack_uint()
    assert secondMagicNumber == extractionMagicNumber, "secondMagicNumber (%r) is not equal to the extractionMagicNumber (%r)" % (secondMagicNumber, extractionMagicNumber)
    versionNumber = up1.unpack_uint()
    assert versionNumber == extractionVersionNumber, "versionNumber (%r) is not the correct version number (%r)" % (versionNumber, extractionVersionNumber)
    voxelSize = up1.unpack_double()
    x0 = up1.unpack_double()
    y0 = up1.unpack_double()
    z0 = up1.unpack_double()
    origin = (x0, y0, z0)
    numSites = up1.unpack_uhyper()

    # Read the header of the checkpoint file, i.e. the number of fields
    # and length of the field header
    infoBytes = 2*sizeOfUint32
    info = checkpointFile.read(infoBytes)
    up2 = xdrlib.Unpacker(info)
    numberOfFields = up2.unpack_uint()
    assert numberOfFields == 1, "numberOfFields (%r) is not equal to 1" % numberOfFields
    lengthOfFieldHeader = up2.unpack_uint()
    assert lengthOfFieldHeader == 32, "lengthOfFieldHeader (%r) is not equal to 32" % lengthOfFieldHeader

    # Read the field header of the checkpoint file.
    fieldHeader = checkpointFile.read(lengthOfFieldHeader)
    up3 = xdrlib.Unpacker(fieldHeader)

    fieldName = up3.unpack_string()
    assert fieldName == "distributions", "fieldName (%r) is not equal to 'distributions'" % fieldName
    numberOfComponents = up3.unpack_uint()
    fieldOffset = up3.unpack_double()
    assert fieldOffset == 0.0, "fieldOffset (%r) is not equal to 0.0" % fieldOffset

    # Read the time step.
    timestepBytes = checkpointFile.read(sizeOfUint64)
    up3 = xdrlib.Unpacker(timestepBytes)
    timestep = up3.unpack_uhyper()
    
    filePosition = checkpointFile.tell()
    assert firstOffset == filePosition, "filePosition (%r) is not equal to firstOffset (%r)" % (filePosition, firstOffset)

    return FileInfo(voxelSize, origin, numSites, numberOfComponents, timestep)

# ******************************************************************************
# Open three checkpoint files for reading. The files must have a common structure.
# Also open an associated file which contains information on where in the
# checkpoint files the data relating to each MPI rank is located.
# ******************************************************************************
def createSummedDistributionFile(checkpointFileName1,
                                    checkpointFileName2,
                                    summedCheckpointFileName):

    sys.stderr.write("createSummedDistributionFile:args:\n")
    sys.stderr.write("checkpointFileName1: " + checkpointFileName1 + "\n")
    sys.stderr.write("checkpointFileName2: " + checkpointFileName2 + "\n")
    sys.stderr.write("summedCheckpointFileName: " + summedCheckpointFileName + "\n")

    # Open a checkpoint file for reading.
    print "Checkpoint file 1:", checkpointFileName1
    checkpointFile1 = open(checkpointFileName1, 'r')

    # Open a checkpoint file for reading.
    print "Checkpoint file 2:", checkpointFileName2
    checkpointFile2 = open(checkpointFileName2, 'r')

    # Open a (checkpoint) offset file for reading.
    offsetFileName = checkpointFileName1[:-3] + "off"
    print "Offset file", offsetFileName
    offsetFile = open(offsetFileName, 'r')

# ******************************************************************************
#  Read all of the data from the offset file.
# ******************************************************************************

    # Read the number of MPI processes on which the simulation was run from
    #   the offset file (a uint64 at the start of the file).
    rawNumProcs = offsetFile.read(sizeOfUint64)

    upO1 = xdrlib.Unpacker(rawNumProcs)
    numProcs = upO1.unpack_uhyper()

    # Read all of the offsets ((numProcs+1)*uint64) from the offset file.
    numberOfOffsets = numProcs+1
    numberOfOffsetBytes = numberOfOffsets*sizeOfUint64
    rawOffsets = offsetFile.read(numberOfOffsetBytes)

    offsetFile.close()

    upO2 = xdrlib.Unpacker(rawOffsets)
    upO2.unpack_uhyper() # Discard.
    distributionOffsets = []
    for r in range(numProcs):
        distributionOffsets.append(upO2.unpack_uhyper())
    
    firstOffset = distributionOffsets[0]

    lastOffset = distributionOffsets[numProcs-1]

# ******************************************************************************
# Read the preambles, headers and timesteps from the checkpoint files.
# ******************************************************************************

    fileInfo1 = readCheckpointFilePreambleEtc(checkpointFile1, firstOffset, lastOffset)
    fileInfo2 = readCheckpointFilePreambleEtc(checkpointFile2, firstOffset, lastOffset)

# ******************************************************************************
# Read all of the position and distribution data from the checkpoint files.
# ******************************************************************************

    readLength = lastOffset - firstOffset

#checkpointFile1.seek(firstOffset)
    data1 = checkpointFile1.read(readLength)
    up03 = xdrlib.Unpacker(data1)

    checkpointFile1.close()

    data2 = checkpointFile2.read(readLength)
    up04 = xdrlib.Unpacker(data2)

    checkpointFile2.close()

# ******************************************************************************
# Open a checkpoint file for writing. This will hold the summed distributions.
# Also open for writing the associated file which contains information on where
# in the preceding file the data relating to each MPI rank is located.
# ******************************************************************************

    # Open a new checkpoint file for writing.
    print "Creating checkpoint file", summedCheckpointFileName
    summedCheckpointFile = open(summedCheckpointFileName, 'w')

    # Open a (checkpoint) offset file for writing.
    checkpointOffsetFileName = summedCheckpointFileName[:-3] + "off"
    print "offset file", checkpointOffsetFileName
    checkpointOffsetFile = open(checkpointOffsetFileName, 'w')

    # Pack the preamble of the checkpoint file, i.e.
    # two magic numbers, a version number, the voxel size,
    # the origin of the grid and the number of sites.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    p01 = xdrlib.Packer()
    p01.pack_uint(HemeLbMagicNumber)
    p01.pack_uint(extractionMagicNumber)
    p01.pack_uint(extractionVersionNumber)
    p01.pack_double(fileInfo1.voxelSize)
    p01.pack_double(fileInfo1.origin[0])
    p01.pack_double(fileInfo1.origin[1])
    p01.pack_double(fileInfo1.origin[2])
    p01.pack_uhyper(fileInfo1.numSites)

    # Pack the header of the checkpoint file, i.e. the number of fields
    # and length of the field header
    numberOfFields = 1
    p01.pack_uint(numberOfFields)
    lengthOfFieldHeader = 32
    p01.pack_uint(lengthOfFieldHeader)

    # Pack the field header of the checkpoint file.
    fieldName = "distributions"
    p01.pack_string(fieldName)
    p01.pack_uint(fileInfo1.numberOfComponents)
    p01.pack_double(0.0)

    summedCheckpointFile.write(p01.get_buffer())

    # Start packing the contents of the offset file.
    p01Off = xdrlib.Packer()
    p01Off.pack_uhyper(numProcs)
    p01Off.pack_uhyper(summedCheckpointFile.tell())

    # Pack the time step.
    p02 = xdrlib.Packer()
    p02.pack_uhyper(fileInfo1.timestep)
    summedCheckpointFile.write(p02.get_buffer())

# ******************************************************************************
# Decode the distributions that have been read and compute the summed values.
# ******************************************************************************

    # A segment consists of a grid point (with three coordinates) and
    # the distributions at that point.
    lengthOfSegment = 3*sizeOfUint32 + fileInfo1.numberOfComponents*sizeOfFloat

    numberOfSegments = readLength / lengthOfSegment
    assert numberOfSegments == fileInfo1.numSites, "numberOfSegments (%r) is not equal to numSites (%r)" % (numberOfSegments, fileInfo1.numSites)

    for r in range(numProcs-1):
        p01Off.pack_uhyper(summedCheckpointFile.tell())
        p03 = xdrlib.Packer()
        startOffset = distributionOffsets[r]
        endOffset = distributionOffsets[r+1]
        nData = (endOffset - startOffset) / lengthOfSegment
        for n in range(nData):
            x = up03.unpack_uint()
            y = up03.unpack_uint()
            z = up03.unpack_uint()
            p03.pack_uint(x)
            p03.pack_uint(y)
            p03.pack_uint(z)
            up04.unpack_uint() # Discard 2nd x value.
            up04.unpack_uint() # Discard 2nd y value.
            up04.unpack_uint() # Discard 2nd z value.
            for c in range(fileInfo1.numberOfComponents):
                distribution1 = up03.unpack_float()
                distribution2 = up04.unpack_float()
                distribution = distribution1 + distribution2
                p03.pack_float(distribution)
        summedCheckpointFile.write(p03.get_buffer())
    p01Off.pack_uhyper(summedCheckpointFile.tell())
    checkpointOffsetFile.write(p01Off.get_buffer())
         
    summedCheckpointFile.close()
    checkpointOffsetFile.close()
