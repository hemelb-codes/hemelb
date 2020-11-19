# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import xdrlib
import sys
from hemeTools.parsers.extraction import ExtractedProperty

header_fields = ['siteCount',
                 'fieldCount',
                 'voxelSizeMetres',
                 'originMetres',
                 'times']

def extraction_loader(filename):
    return ExtractedProperty(filename)

def extraction_parser(content, pattern):
    if pattern in header_fields:
        return getattr(content, pattern)
    if isinstance(pattern, str):
        time = 'final'
        name = pattern
    else:
        time, name = pattern
        if isinstance(time, str) and time[0] == 'i':
            time = content.times[int(time[1:])]
    if time == 'final':
        time = content.times[-1]
    fields = content.GetByTimeStep(time)
    return fields.field(name)
