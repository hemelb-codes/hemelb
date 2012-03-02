#! /usr/bin/env python
import xdrlib
import sys

def unpack(filename):
  file = open(filename, 'r')
  x = xdrlib.Unpacker(file.read())

  fieldCount = x.unpack_uint()
  fieldNames = [x.unpack_string() for fieldNum in xrange(fieldCount)]
  fieldLengths = [x.unpack_uint() for fieldNum in xrange(fieldCount)]

  siteCount = x.unpack_uhyper()

  print "-----\nOpened file \"%s\"...\n-----" % filename
  print "File has %s fields: %s, for %s sites" % (fieldCount, fieldNames, siteCount)
  print "Each of these fields has %s floats" % fieldLengths
 
  while x.get_position() < (len(x.get_buffer())):
    iterationNumber = x.unpack_hyper()

    print "\nIteration %s\n-----" % iterationNumber

    for siteNum in xrange(siteCount):
      position = [x.unpack_uhyper() for fieldNum in xrange(3)]
      field = [[x.unpack_float() for dimension in xrange(fieldLengths[fieldNum])] for fieldNum in xrange(fieldCount)]
      print "Site at %s: %s = %s" % (position, fieldNames, field)

if __name__ == "__main__":
  unpack(sys.argv[1])
