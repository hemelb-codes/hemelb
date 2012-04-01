#! /usr/bin/env python
import xdrlib
import sys
from hemeTools.parsers import HemeLbMagicNumber

ExtractionMagicNumber = 0x78747278

def unpack(filename):
  file = open(filename, 'r')
  x = xdrlib.Unpacker(file.read())

  hlbMagicNumber = x.unpack_uint()
  extractionMagicNumber = x.unpack_uint()
  extractionVersionNumber = x.unpack_uint()

  if hlbMagicNumber != HemeLbMagicNumber:
    print "Expected HemeLb magic number %i, got %i" % (HemeLbMagicNumber, hlbMagicNumber)
    return
  
  if extractionMagicNumber != ExtractionMagicNumber:
    print "Expected Extraction magic number %i, got %i" % (ExtractionMagicNumber, extractionMagicNumber)
    return

  print "Extraction version number %i" % extractionVersionNumber

  fieldCount = x.unpack_uint()
  fieldNames = [x.unpack_string() for fieldNum in xrange(fieldCount)]
  fieldLengths = [x.unpack_uint() for fieldNum in xrange(fieldCount)]

  siteCount = x.unpack_uhyper()
  voxelSize = x.unpack_float()
  origin = [x.unpack_float() for fieldNum in xrange(3)]

  print "-----\nOpened file \"%s\"...\n-----" % filename
  print "File has %s fields: %s, for %s sites" % (fieldCount, fieldNames, siteCount)
  print "Each of these fields has %s floats" % fieldLengths
  print "Geometry origin: %s, voxel size: %s" % (origin, voxelSize)
 
  while x.get_position() < (len(x.get_buffer())):
    iterationNumber = x.unpack_hyper()

    print "\nIteration %s\n-----" % iterationNumber

    for siteNum in xrange(siteCount):
      position = [x.unpack_uhyper() for fieldNum in xrange(3)]
      field = [[x.unpack_float() for dimension in xrange(fieldLengths[fieldNum])] for fieldNum in xrange(fieldCount)]
      print "Site at %s: %s = %s" % (position, fieldNames, field)

if __name__ == "__main__":
  unpack(sys.argv[1])
