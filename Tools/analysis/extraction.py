import xdrlib
import sys
from hemeTools.parsers.extraction import ExtractedProperty

header_fields=['siteCount','fieldCount','voxelSizeMetres','originMetres', 'times']

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
    if type(time) == str and time[0] == 'i':
      time = content.times[int(time[1:])]
  if time=='final':
    time=content.times[-1]
  fields=content.GetByTimeStep(time)
  return fields.field(name)
