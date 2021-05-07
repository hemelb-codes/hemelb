# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import xdrlib
import sys
import os

def createCoarseCheckpointFile(fineCheckpointFileName,
                               coarseRankFileName,
                               coarseCheckpointFileName):

    sys.stderr.write("createCoarseCheckpointFile:args:\n")
    sys.stderr.write("fineCheckpointFileName: " + fineCheckpointFileName + "\n")
    sys.stderr.write("coarseRankFileName: " + coarseRankFileName + "\n")
    sys.stderr.write("coarseCheckpointFileName: " + coarseCheckpointFileName + "\n")

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
# Open a fine checkpoint file.
# Also open the associated file which contains information on where in the
# preceding file the data relating to each MPI rank is located.
# ******************************************************************************

    # Open a fine checkpoint file for reading.
    print "Fine checkpoint file", fineCheckpointFileName
    fineCheckpointFile = open(fineCheckpointFileName, 'r')

    # Open a fine (checkpoint) offset file for reading.
    fineOffsetFileName = fineCheckpointFileName[:-3] + "off"
    print "Fine offset file", fineOffsetFileName
    fineOffsetFile = open(fineOffsetFileName, 'r')

# ******************************************************************************
#  Read all of the data from the fine offset file.
# ******************************************************************************

    # Read the number of MPI processes on which the fine simulation was run from
    #   the fine offset file (a uint64 at the start of the file).
    rawNumFineProcs = fineOffsetFile.read(sizeOfUint64)

    upFCO1 = xdrlib.Unpacker(rawNumFineProcs)
    numFineProcs = upFCO1.unpack_uhyper()

    # Read all of the fine offsets ((numFineProcs+1)*uint64) from the offset file.
    numberOfOffsets = numFineProcs+1
    numberOfOffsetBytes = numberOfOffsets*sizeOfUint64
    rawFineOffsets = fineOffsetFile.read(numberOfOffsetBytes)

    fineOffsetFile.close()

    upFCO2 = xdrlib.Unpacker(rawFineOffsets)
    upFCO2.unpack_uhyper() # Discard.
    firstOffset = upFCO2.unpack_uhyper()

    upFCO2.set_position(numFineProcs*sizeOfUint64)
    lastOffset = upFCO2.unpack_uhyper()

# ******************************************************************************
# Read all of the data from the fine checkpoint file.
# A checkpoint file is an extraction file..
# ******************************************************************************

    # Read the preamble of the checkpoint file, i.e.
    # read two magic numbers, a version number, the voxel size,
    # the origin of the grid and the number of sites.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    preamble = fineCheckpointFile.read(preambleBytes)
    upFC1 = xdrlib.Unpacker(preamble)
    firstMagicNumber = upFC1.unpack_uint()
    assert firstMagicNumber == HemeLbMagicNumber, "firstMagicNumber (%r) is not equal to the HemeLbMagicNumber" % firstMagicNumber
    secondMagicNumber = upFC1.unpack_uint()
    assert secondMagicNumber == extractionMagicNumber, "secondMagicNumber (%r) is not equal to the extractionMagicNumber (%r)" % (secondMagicNumber, extractionMagicNumber)
    versionNumber = upFC1.unpack_uint()
    assert versionNumber == extractionVersionNumber, "versionNumber (%r) is not the correct version number (%r)" % (versionNumber, extractionVersionNumber)
    fineVoxelSize = upFC1.unpack_double()
    x0F = upFC1.unpack_double()
    y0F = upFC1.unpack_double()
    z0F = upFC1.unpack_double()
    fineOrigin = (x0F, y0F, z0F)
    numFineSites = upFC1.unpack_uhyper()

    # Read the header of the checkpoint file, i.e. the number of fields
    # and length of the field header
    infoBytes = 2*sizeOfUint32
    info = fineCheckpointFile.read(infoBytes)
    upFC2 = xdrlib.Unpacker(info)
    numberOfFields = upFC2.unpack_uint()
    assert numberOfFields == 1, "numberOfFields (%r) is not equal to 1" % numberOfFields
    lengthOfFieldHeader = upFC2.unpack_uint()
    assert lengthOfFieldHeader == 32, "lengthOfFieldHeader (%r) is not equal to 32" % lengthOfFieldHeader

    # Read the field header of the checkpoint file.
    fieldHeader = fineCheckpointFile.read(lengthOfFieldHeader)
    upFC3 = xdrlib.Unpacker(fieldHeader)

    fieldName = upFC3.unpack_string()
    assert fieldName == "distributions", "fieldName (%r) is not equal to 'distributions'" % fieldName
    numberOfComponents = upFC3.unpack_uint()
    fieldOffset = upFC3.unpack_double()
    assert fieldOffset == 0.0, "fieldOffset (%r) is not equal to 0.0" % fieldOffset

    # Read the time step.
    timestepBytes = fineCheckpointFile.read(sizeOfUint64)
    upFC3 = xdrlib.Unpacker(timestepBytes)
    timestep = upFC3.unpack_uhyper()
    
    filePosition = fineCheckpointFile.tell()
    assert firstOffset == filePosition, "filePosition (%r) is not equal to firstOffset (%r)" % (filePosition, firstOffset)

    # Read all of the fine position and distribution data from the checkpoint file.
    readLength = lastOffset - firstOffset
    fineCheckpointFile.seek(firstOffset)
    data = fineCheckpointFile.read(readLength)
    upFC4 = xdrlib.Unpacker(data)

    fineCheckpointFile.close()

# ******************************************************************************
# Decode the distributions that have been read and insert those corresponding
# to coarse grid points into a dictionary.
# ******************************************************************************

    # A segment consists of a grid point (with three coordinates) and
    # the distributions at that point.
    lengthOfSegment = 3*sizeOfUint32 + numberOfComponents*sizeOfFloat

    numberOfSegments = readLength / lengthOfSegment
    assert numberOfSegments == numFineSites, "numberOfSegments (%r) is not equal to numFineSites (%r)" % (numberOfSegments, numFineSites)

    coarseFlowField = {} # Coarse flow field.
    for i in range(numberOfSegments):
        x = upFC4.unpack_uint()
        y = upFC4.unpack_uint()
        z = upFC4.unpack_uint()
        distributions = []
        for j in range(numberOfComponents):
            distribution = upFC4.unpack_float()
            distributions.append(distribution)
        numberOfOddIndices = x%2 + y%2 + z%2
        if numberOfOddIndices == 0:
            point = (x/2, y/2, z/2) # Create a coarse grid point.
            # Copy the fine distributions to the coarse distributions.
            coarseFlowField[point] = distributions

# ******************************************************************************
# Open the file containing information about the MPI rank on which each coarse
# grid point is located.
# Also open the associated file which contains information on where in the
# preceding file the data relating to each MPI rank is located.
# ******************************************************************************

    # Open a coarse rank file for reading.
    print "Coarse rank file", coarseRankFileName
    coarseRankFile = open(coarseRankFileName, 'r')

    # Open a coarse (rank) offset file for reading.
    coarseRankOffsetFileName = coarseRankFileName[:-3] + "off"
    print "Coarse offset file", coarseRankOffsetFileName
    coarseRankOffsetFile = open(coarseRankOffsetFileName, 'r')

# ******************************************************************************
#  Read all of the data from the coarse offset file.
# ******************************************************************************

    # Read the number of MPI processes on which the coarse simulation was run from
    #   the coarse offset file (a uint64 at the start of the file).
    rawNumCoarseProcs = coarseRankOffsetFile.read(sizeOfUint64)

    upCRO1 = xdrlib.Unpacker(rawNumCoarseProcs)
    numCoarseProcs = upCRO1.unpack_uhyper()

    # Read all of the coarse offsets ((numCoarseProcs+1)*uint64) from the offset file.
    numberOfOffsets = numCoarseProcs+1
    numberOfOffsetBytes = numberOfOffsets*sizeOfUint64
    rawCoarseOffsets = coarseRankOffsetFile.read(numberOfOffsetBytes)
    upCRO2 = xdrlib.Unpacker(rawCoarseOffsets)

    coarseRankOffsetFile.close()

    upCRO2.unpack_uhyper() # Discard.
    coarseRankOffsets = []
    for r in range(numCoarseProcs):
        coarseRankOffsets.append(upCRO2.unpack_uhyper())

# ******************************************************************************
# Now process the coarse rank file.
# Once the the number of MPI processes on which the coarse simulation was run has
# been read the grid point data is processed rank by rank.
# ******************************************************************************

    # Read the preamble of the coarse rank file.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    preamble = coarseRankFile.read(preambleBytes)
    upCR1 = xdrlib.Unpacker(preamble)
    firstMagicNumber = upCR1.unpack_uint()
    assert firstMagicNumber == HemeLbMagicNumber, "firstMagicNumber (%r) is not equal to the HemeLbMagicNumber" % firstMagicNumber
    secondMagicNumber = upCR1.unpack_uint()
    assert secondMagicNumber == extractionMagicNumber, "secondMagicNumber (%r) is not the correct version number (%r)" % (secondMagicNumber, extractionMagicNumber)
    versionNumber = upCR1.unpack_uint()
    assert versionNumber == extractionVersionNumber, "versionNumber (%r) is not the correct version number (%r)" % (versionNumber, extractionVersionNumber)
    coarseVoxelSize = upCR1.unpack_double()
    assert coarseVoxelSize == 2*fineVoxelSize, "coarseVoxelSize (%r) is not twice the fineVoxelSize (%r)" % (fineVoxelSize, coarseVoxelSize)
    x0C = upCR1.unpack_double()
    y0C = upCR1.unpack_double()
    z0C = upCR1.unpack_double()
    coarseOrigin = (x0C, y0C, z0C)
    numCoarseSites = upCR1.unpack_uhyper()

    # Read the header of the coarse rank file.
    infoBytes = 2*sizeOfUint32
    info = coarseRankFile.read(infoBytes)
    upCR2 = xdrlib.Unpacker(info)
    numberOfFields = upCR2.unpack_uint()
    assert numberOfFields == 1, "numberOfFields (%r) is not equal to 1" % numberOfFields
    lengthOfFieldHeader = upCR2.unpack_uint()
    assert lengthOfFieldHeader == 24, "lengthOfFieldHeader (%r) is not equal to 24" % lengthOfFieldHeader

    # Read the field header of the coarse rank file.
    fieldHeader = coarseRankFile.read(lengthOfFieldHeader)
    upCR3 = xdrlib.Unpacker(fieldHeader)

    fieldName = upCR3.unpack_string()
    assert fieldName == "mpirank", "fieldName (%r) is not equal to 'mpirank'" % fieldName
    numberOfFloats = upCR3.unpack_uint()
    assert numberOfFloats == 1, "numberOfFloats (%r) is not equal to 1" % numberOfFloats
    offset = upCR3.unpack_double()
    assert offset == 0.0, "offset (%r) is not equal to 0.0" % offset

    timestepBytes = coarseRankFile.read(sizeOfUint64)
    upCR4 = xdrlib.Unpacker(timestepBytes)
    upCR4.unpack_uhyper() # Discard the timestep.

    filePosition = coarseRankFile.tell()
    firstCoarseRankOffset = coarseRankOffsets[0]
    lastCoarseRankOffset = coarseRankOffsets[numCoarseProcs-1]
    assert firstCoarseRankOffset == filePosition, "filePosition (%r) is not equal to firstCoarseRankOffset (%r)" % (filePosition, firstCoarseRankOffset)

    # Read all of the coarse position and rank data from the rank file.
    readLength = lastCoarseRankOffset - firstCoarseRankOffset
    coarseRankFile.seek(firstCoarseRankOffset)
    data = coarseRankFile.read(readLength)
    upCR5 = xdrlib.Unpacker(data)

    coarseRankFile.close()

# ******************************************************************************
# Open a coarse checkpoint file for writing.
# Also open for writing the associated file which contains information on where
# in the preceding file the data relating to each MPI rank is located.
# ******************************************************************************

    # Open a coarse checkpoint file for writing.
    print "Coarse checkpoint file", coarseCheckpointFileName
    coarseCheckpointFile = open(coarseCheckpointFileName, 'w')

    # Open a coarse (checkpoint) offset file for writing.
    coarseCheckpointOffsetFileName = coarseCheckpointFileName[:-3] + "off"
    print "Coarse offset file", coarseCheckpointOffsetFileName
    coarseCheckpointOffsetFile = open(coarseCheckpointOffsetFileName, 'w')

    # Pack the preamble of the checkpoint file, i.e.
    # two magic numbers, a version number, the voxel size,
    # the origin of the grid and the number of sites.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    pCC1 = xdrlib.Packer()
    pCC1.pack_uint(HemeLbMagicNumber)
    pCC1.pack_uint(extractionMagicNumber)
    pCC1.pack_uint(versionNumber)
    pCC1.pack_double(coarseVoxelSize)
    pCC1.pack_double(x0C)
    pCC1.pack_double(y0C)
    pCC1.pack_double(z0C)
    pCC1.pack_uhyper(numCoarseSites)

    # Pack the header of the checkpoint file, i.e. the number of fields
    # and length of the field header
    numberOfFields = 1
    pCC1.pack_uint(numberOfFields)
    lengthOfFieldHeader = 32
    pCC1.pack_uint(lengthOfFieldHeader)

    # Pack the field header of the checkpoint file.
    fieldName = "distributions"
    pCC1.pack_string(fieldName)
    pCC1.pack_uint(numberOfComponents)
    pCC1.pack_double(fieldOffset)

    coarseCheckpointFile.write(pCC1.get_buffer())

    pCCO1 = xdrlib.Packer()
    pCCO1.pack_uhyper(numCoarseProcs)
    pCCO1.pack_uhyper(coarseCheckpointFile.tell())

    # Pack the time step.
    pCC2 = xdrlib.Packer()
    pCC2.pack_uhyper(timestep) # Using the time step from the fine checkpoint file.
    coarseCheckpointFile.write(pCC2.get_buffer())

# ******************************************************************************
# Decode the coarse rank data that have been read and pack the corresponding
# distributions in preparation for writing them to a file.
# ******************************************************************************

    # A segment consists of a grid point (with three coordinates) and
    # the rank on which that point is located.
    lengthOfSegment = 3*sizeOfUint32 + sizeOfFloat

    for r in range(numCoarseProcs-1):
        pCCO1.pack_uhyper(coarseCheckpointFile.tell())
        pCC3 = xdrlib.Packer()
        startOffset = coarseRankOffsets[r]
        endOffset = coarseRankOffsets[r+1]
        nData = (endOffset - startOffset) / lengthOfSegment
        for n in range(nData):
            x = upCR5.unpack_uint()
            y = upCR5.unpack_uint()
            z = upCR5.unpack_uint()
            upCR5.unpack_float() # Discard the rank.
            pCC3.pack_uint(x)
            pCC3.pack_uint(y)
            pCC3.pack_uint(z)
            point = (x, y, z) # Create a coarse grid point.
            distribution = coarseFlowField[point]
            for c in distribution:
                pCC3.pack_float(c)
        coarseCheckpointFile.write(pCC3.get_buffer())
    pCCO1.pack_uhyper(coarseCheckpointFile.tell())
    coarseCheckpointOffsetFile.write(pCCO1.get_buffer())
         
    coarseCheckpointFile.close()
    coarseCheckpointOffsetFile.close()
