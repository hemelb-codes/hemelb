# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 
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
