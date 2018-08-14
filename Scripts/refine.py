# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import xdrlib
import sys
import os

# ******************************************************************************
# Define a function for later use.
# ******************************************************************************

def computeFlow(point, coarseFlowField, neighbouringPoints, numberOfComponents):
    numberOfNeighbours = len(neighbouringPoints)
    fineFlow = []
    for i in range(numberOfComponents):
        fineFlow.append(0.0)
    for neighbour in neighbouringPoints:
        coarseFlow = coarseFlowField[neighbour]
        for j in range(numberOfComponents):
            fineFlow[j] = fineFlow[j] + coarseFlow[j]
    for k in range(numberOfComponents):
        fineFlow[k] = fineFlow[k]/numberOfNeighbours
    return fineFlow
            
def createFineCheckpointFile(coarseCheckpointFileName,
                             fineRankFileName,
                             fineCheckpointFileName):

    sys.stderr.write("createFineCheckpointFile:args:\n")
    sys.stderr.write("coarseCheckpointFileName: " + coarseCheckpointFileName + "\n")
    sys.stderr.write("fineRankFileName: " + fineRankFileName + "\n")
    sys.stderr.write("fineCheckpointFileName: " + fineCheckpointFileName + "\n")

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
# Open a coarse checkpoint file.
# Also open the associated file which contains information on where in the
# preceding file the data relating to each MPI rank is located.
# ******************************************************************************

    # Open a coarse checkpoint file for reading.
    print "Coarse checkpoint file", coarseCheckpointFileName
    coarseCheckpointFile = open(coarseCheckpointFileName, 'r')

    # Open a coarse (checkpoint) offset file for reading.
    coarseOffsetFileName = coarseCheckpointFileName[:-3] + "off"
    print "Coarse offset file", coarseOffsetFileName
    coarseOffsetFile = open(coarseOffsetFileName, 'r')

# ******************************************************************************
#  Read all of the data from the coarse offset file.
# ******************************************************************************

    # Read the number of MPI processes on which the coarse simulation was run from
    #   the coarse offset file (a uint64 at the start of the file).
    rawNumCoarseProcs = coarseOffsetFile.read(sizeOfUint64)

    upCCO1 = xdrlib.Unpacker(rawNumCoarseProcs)
    numCoarseProcs = upCCO1.unpack_uhyper()

    # Read all of the coarse offsets ((numCoarseProcs+1)*uint64) from the offset file.
    numberOfOffsets = numCoarseProcs+1
    numberOfOffsetBytes = numberOfOffsets*sizeOfUint64
    rawCoarseOffsets = coarseOffsetFile.read(numberOfOffsetBytes)

    coarseOffsetFile.close()

    upCCO2 = xdrlib.Unpacker(rawCoarseOffsets)
    upCCO2.unpack_uhyper() # Discard.
    firstOffset = upCCO2.unpack_uhyper()

    upCCO2.set_position(numCoarseProcs*sizeOfUint64)
    lastOffset = upCCO2.unpack_uhyper()

# ******************************************************************************
# Read all of the data from the coarse checkpoint file.
# A checkpoint file is an extraction file..
# ******************************************************************************

    # Read the preamble of the checkpoint file, i.e.
    # read two magic numbers, a version number, the voxel size,
    # the origin of the grid and the number of sites.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    preamble = coarseCheckpointFile.read(preambleBytes)
    upCC1 = xdrlib.Unpacker(preamble)
    firstMagicNumber = upCC1.unpack_uint()
    assert firstMagicNumber == HemeLbMagicNumber, "firstMagicNumber (%r) is not equal to the HemeLbMagicNumber" % firstMagicNumber
    secondMagicNumber = upCC1.unpack_uint()
    assert secondMagicNumber == extractionMagicNumber, "secondMagicNumber (%r) is not equal to the extractionMagicNumber" % (secondMagicNumber, extractionMagicNumber)
    versionNumber = upCC1.unpack_uint()
    assert versionNumber == extractionVersionNumber, "versionNumber (%r) is not the correct version number (%r)" % (versionNumber, extractionVersionNumber)
    coarseVoxelSize = upCC1.unpack_double()
    x0C = upCC1.unpack_double()
    y0C = upCC1.unpack_double()
    z0C = upCC1.unpack_double()
    coarseOrigin = (x0C, y0C, z0C)
    numCoarseSites = upCC1.unpack_uhyper()

    # Read the header of the checkpoint file, i.e. the number of fields
    # and length of the field header
    infoBytes = 2*sizeOfUint32
    info = coarseCheckpointFile.read(infoBytes)
    upCC2 = xdrlib.Unpacker(info)
    numberOfFields = upCC2.unpack_uint()
    assert numberOfFields == 1, "numberOfFields (%r) is not equal to 1" % numberOfFields
    lengthOfFieldHeader = upCC2.unpack_uint()
    assert lengthOfFieldHeader == 32, "lengthOfFieldHeader (%r) is not equal to 32" % lengthOfFieldHeader

    # Read the field header of the checkpoint file.
    fieldHeader = coarseCheckpointFile.read(lengthOfFieldHeader)
    upCC3 = xdrlib.Unpacker(fieldHeader)

    fieldName = upCC3.unpack_string()
    assert fieldName == "distributions", "fieldName (%r) is not equal to 'distributions'" % fieldName
    numberOfComponents = upCC3.unpack_uint()
    fieldOffset = upCC3.unpack_double()
    assert fieldOffset == 0.0, "fieldOffset (%r) is not equal to 0.0" % fieldOffset

    # Read the time step.
    timestepBytes = coarseCheckpointFile.read(sizeOfUint64)
    upCC4 = xdrlib.Unpacker(timestepBytes)
    timestep = upCC4.unpack_uhyper()
    
    filePosition = coarseCheckpointFile.tell()
    assert firstOffset == filePosition, "filePosition (%r) is not equal to firstOffset (%r)" % (filePosition, firstOffset)

    # Read all of the coarse position and distribution data from the checkpoint file.
    readLength = lastOffset - firstOffset
    coarseCheckpointFile.seek(firstOffset)
    coarseCheckpointData = coarseCheckpointFile.read(readLength)
    upCC5 = xdrlib.Unpacker(coarseCheckpointData)

    coarseCheckpointFile.close()

# ******************************************************************************
# Decode the distributions that have been read and insert them into a dictionary.
# ******************************************************************************

    # A segment consists of a grid point (with three coordinates) and
    # the distributions at that point.
    lengthOfSegment = 3*sizeOfUint32 + numberOfComponents*sizeOfFloat

    numberOfSegments = readLength / lengthOfSegment
    assert numberOfSegments == numCoarseSites, "numberOfSegments (%r) is not equal to numCoarseSites (%r)" % (numberOfSegments, numCoarseSites)

    coarseFlowField = {} # Coarse flow field on fine grid.
    for i in range(numberOfSegments):
        x = upCC5.unpack_uint()
        y = upCC5.unpack_uint()
        z = upCC5.unpack_uint()
        point = (2*x, 2*y, 2*z) # Create a fine grid point.
        # Copy the coarse distributions to the fine distributions.
        distributions = []
        for j in range(numberOfComponents):
            distribution = upCC5.unpack_float()
            distributions.append(distribution)
        coarseFlowField[point] = distributions

# ******************************************************************************
# Open the file containing information about the MPI rank on which each fine
# grid point is located.
# Also open the associated file which contains information on where in the
# preceding file the data relating to each MPI rank is located.
# ******************************************************************************

    # Open a fine rank file for reading.
    print "Fine rank file", fineRankFileName
    fineRankFile = open(fineRankFileName, 'r')

    # Open a fine (rank) offset file for reading.
    fineOffsetFileName = fineRankFileName[:-3] + "off"
    print "Fine offset file", fineOffsetFileName
    fineOffsetFile = open(fineOffsetFileName, 'r')

# ******************************************************************************
#  Read all of the data from the fine offset file.
# ******************************************************************************

# Read the number of MPI processes on which the fine simulation was run from
#   the fine offset file (a uint64 at the start of the file).
    rawNumFineProcs = fineOffsetFile.read(sizeOfUint64)

    upFRO1 = xdrlib.Unpacker(rawNumFineProcs)
    numFineProcs = upFRO1.unpack_uhyper()

    # Read all of the fine offsets ((numFineProcs+1)*uint64) from the offset file.
    numberOfOffsets = numFineProcs+1
    numberOfOffsetBytes = numberOfOffsets*sizeOfUint64
    rawFineOffsets = fineOffsetFile.read(numberOfOffsetBytes)
    upFRO2 = xdrlib.Unpacker(rawFineOffsets)

    fineOffsetFile.close()

    upFRO2.unpack_uhyper() # Discard.

    fineRankOffsets = []
    for r in range(numFineProcs):
        fineRankOffsets.append(upFRO2.unpack_uhyper())

    upCCO2.set_position(numFineProcs*sizeOfUint64)

# ******************************************************************************
# Now process the fine rank file.
# Once the the number of MPI processes on which the fine simulation was run has
# been read the grid point data is processed rank by rank.
# ******************************************************************************

    # Read the preamble of the fine rank file.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    preamble = fineRankFile.read(preambleBytes)
    upFR1 = xdrlib.Unpacker(preamble)
    firstMagicNumber = upFR1.unpack_uint()
    assert firstMagicNumber == HemeLbMagicNumber, "firstMagicNumber (%r) is not equal to the HemeLbMagicNumber" % firstMagicNumber
    secondMagicNumber = upFR1.unpack_uint()
    assert secondMagicNumber == extractionMagicNumber, "secondMagicNumber (%r) is not the correct version number (%r)" % (secondMagicNumber, extractionMagicNumber)
    versionNumber = upFR1.unpack_uint()
    assert versionNumber == extractionVersionNumber, "versionNumber (%r) is not the correct version number (%r)" % (versionNumber, extractionVersionNumber)
    fineVoxelSize = upFR1.unpack_double()
    assert coarseVoxelSize == 2*fineVoxelSize, "coarseVoxelSize (%r) is not twice the fineeVoxelSize (%r)" % (coarseVoxelSize, fineVoxelSize)
    x0F = upFR1.unpack_double()
    y0F = upFR1.unpack_double()
    z0F = upFR1.unpack_double()
    fineOrigin = (x0F, y0F, z0F)
    numFineSites = upFR1.unpack_uhyper()

    # Read the header of the fine rank file.
    infoBytes = 2*sizeOfUint32
    info = fineRankFile.read(infoBytes)
    upFR2 = xdrlib.Unpacker(info)
    numberOfFields = upFR2.unpack_uint()
    assert numberOfFields == 1, "numberOfFields (%r) is not equal to 1" % numberOfFields
    lengthOfFieldHeader = upFR2.unpack_uint()
    assert lengthOfFieldHeader == 24, "lengthOfFieldHeader (%r) is not equal to 24" % lengthOfFieldHeader

    # Read the field header of the fine rank file.
    fieldHeader = fineRankFile.read(lengthOfFieldHeader)
    upFR3 = xdrlib.Unpacker(fieldHeader)

    fieldName = upFR3.unpack_string()
    assert fieldName == "mpirank", "fieldName (%r) is not equal to 'mpirank'" % fieldName
    numberOfFloats = upFR3.unpack_uint()
    assert numberOfFloats == 1, "numberOfFloats (%r) is not equal to 1" % numberOfFloats
    offset = upFR3.unpack_double()
    assert offset == 0.0, "offset (%r) is not equal to 0.0" % offset

    timestepBytes = fineRankFile.read(sizeOfUint64)
    upFR4 = xdrlib.Unpacker(timestepBytes)
    upFR4.unpack_uhyper() # Discard the timestep.

    filePosition = fineRankFile.tell()
    firstFineRankOffset = fineRankOffsets[0]
    lastFineRankOffset = fineRankOffsets[numFineProcs-1]
    assert firstFineRankOffset == filePosition, "filePosition (%r) is not equal to firstFineRankOffset (%r)" % (filePosition, firstFineRankOffset)

    # Read all of the fine position and rank data from the rank file.
    readLength = lastFineRankOffset - firstFineRankOffset
    fineRankFile.seek(firstFineRankOffset)
    fineRankData = fineRankFile.read(readLength)
    upFR5 = xdrlib.Unpacker(fineRankData)

    fineRankFile.close()

# ******************************************************************************
# Open a fine checkpoint file for writing.
# Also open for writing the associated file which contains information on where
# in the preceding file the data relating to each MPI rank is located.
# ******************************************************************************

    # Open a fine checkpoint file for writing.
    print "Fine checkpoint file", fineCheckpointFileName
    fineCheckpointFile = open(fineCheckpointFileName, 'w')

    # Open a fine (checkpoint) offset file for writing.
    fineCheckpointOffsetFileName = fineCheckpointFileName[:-3] + "off"
    print "Fine offset file", fineCheckpointOffsetFileName
    fineCheckpointOffsetFile = open(fineCheckpointOffsetFileName, 'w')

    # Pack the preamble of the checkpoint file, i.e.
    # two magic numbers, a version number, the voxel size,
    # the origin of the grid and the number of sites.
    preambleBytes = 3*sizeOfUint32 + 4*sizeOfDouble + sizeOfUint64
    pFC1 = xdrlib.Packer()
    pFC1.pack_uint(HemeLbMagicNumber)
    pFC1.pack_uint(extractionMagicNumber)
    pFC1.pack_uint(versionNumber)
    pFC1.pack_double(fineVoxelSize)
    pFC1.pack_double(x0C)
    pFC1.pack_double(y0C)
    pFC1.pack_double(z0C)
    pFC1.pack_uhyper(numFineSites)

    # Pack the header of the checkpoint file, i.e. the number of fields
    # and length of the field header
    numberOfFields = 1
    pFC1.pack_uint(numberOfFields)
    lengthOfFieldHeader = 32
    pFC1.pack_uint(lengthOfFieldHeader)

    # Pack the field header of the checkpoint file.
    fieldName = "distributions"
    pFC1.pack_string(fieldName)
    pFC1.pack_uint(numberOfComponents)
    pFC1.pack_double(fieldOffset)

    fineCheckpointFile.write(pFC1.get_buffer())

    pFCO1 = xdrlib.Packer()
    pFCO1.pack_uhyper(numFineProcs)
    pFCO1.pack_uhyper(fineCheckpointFile.tell())

    # Pack the time step.
    pFC2 = xdrlib.Packer()
    pFC2.pack_uhyper(timestep) # Using the time step from the fine checkpoint file.
    fineCheckpointFile.write(pFC2.get_buffer())

# ******************************************************************************
# Decode the rank data that have been read and insert entries into the dictionary.
# ******************************************************************************

    # A segment consists of a grid point (with three coordinates) and
    # the rank on which that point is located.
    lengthOfSegment = 3*sizeOfUint32 + sizeOfFloat

    maxNumNeighbours = [0, 0, 0, 0]
    minX = 1000000
    maxX = -1
    minY = 1000000
    maxY = -1
    minZ = 1000000
    maxZ = -1
    numBoundaryPoints = 0
    numAnomPoints = 0
    neighbouringCoarsePoints = {}
    for r in range(numFineProcs-1):
        pFCO1.pack_uhyper(fineCheckpointFile.tell())
        pFC3 = xdrlib.Packer()
        startOffset = fineRankOffsets[r]
        endOffset = fineRankOffsets[r+1]
        nData = (endOffset - startOffset) / lengthOfSegment
        for n in range(nData):
            x = upFR5.unpack_uint()
            y = upFR5.unpack_uint()
            z = upFR5.unpack_uint()
            upFR5.unpack_float() # Discard the rank.
            pFC3.pack_uint(x)
            pFC3.pack_uint(y)
            pFC3.pack_uint(z)
            point = (x, y, z) # Create a fine grid point.
            numberOfOddIndices = x%2 + y%2 + z%2
            minX = min(x, minX)
            maxX = max(x, maxX)
            minY = min(y, minY)
            maxY = max(y, maxY)
            minZ = min(z, minZ)
            maxZ = max(z, maxZ)
            neighbouringPoints = []
            for i in range(-1, 2):
                for j in range(-1, 2):
                    for k in range(-1, 2):
                        neighbouringPoint = (x+i, y+j, z+k)
                        isCoarsePoint = coarseFlowField.get(neighbouringPoint, False)
                        if isCoarsePoint is not False:
                            neighbouringPoints.append(neighbouringPoint)
            maxNumNeighbours[numberOfOddIndices] = max(maxNumNeighbours[numberOfOddIndices], len(neighbouringPoints))
            neighbouringCoarsePoints[point] = neighbouringPoints
            if len(neighbouringPoints) != pow(2, numberOfOddIndices):
                if not ((x==1) or (x==15) or (y==1) or (y==15) or (z==1) or (z==87)):
                    numAnomPoints = numAnomPoints + 1
                else:
                    numBoundaryPoints = numBoundaryPoints + 1
            fineFlow = computeFlow(point, coarseFlowField, neighbouringPoints, numberOfComponents)
            for c in fineFlow:
                pFC3.pack_float(c)
        fineCheckpointFile.write(pFC3.get_buffer())

    pFCO1.pack_uhyper(fineCheckpointFile.tell())
    fineCheckpointOffsetFile.write(pFCO1.get_buffer())
