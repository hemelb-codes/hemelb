import xdrlib
import sys
from hemeTools.parsers.extraction import ExtractedProperty

header_fields=['siteCount','fieldCount','voxelSizeMetres','originMetres']

def extraction_loader(filename):
  return ExtractedProperty(filename)

def extraction_parser(content,pattern):
  if pattern in header_fields:
    return getattr(content,pattern)
  if type(pattern)==str:
    time='final'
    name=pattern
  else:
    time,name=pattern
  if time=='final':
    time=content.times[-1]
  fields=content.GetByTimeStep(time)
  return [getattr(row,name) for row in fields]
