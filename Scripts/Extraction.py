#! /usr/bin/env python
import xdrlib
import sys
from hemeTools.parsers.extraction import ExtractedProperty

def unpack(filename):
    propFile = ExtractedProperty(filename)

    print '# Dump of file "{}"'.format(filename)
    print '# File has {} sites.'.format(propFile.siteCount)
    print '# File has {} fields:'.format(propFile.fieldCount)
    for name, xdrType, memType, length, offset in propFile._fieldSpec:
        print '#     "{0}", length {1}'.format(name, length)
    print '# Geometry origin = {} m'.format(propFile.originMetres)
    print '# Voxel size = {} m'.format(propFile.voxelSizeMetres)

    header  = '# '+ ', '.join(name for name, xdrType, memType, length, offset in propFile._fieldSpec)
    print header
    
    for t in propFile.times:
        fields = propFile.GetByTimeStep(t)
        print "# Timestep {:d}".format(t)
        
        for row in fields:
            print ', '.join(str(getattr(row, name)) for name, xdrType, memType, length, offset in propFile._fieldSpec)

        print ""

if __name__ == "__main__":
  unpack(sys.argv[1])
