# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import os.path
import re
import numpy as np
from enthought import units as U

from .orderer import continuousOrder
from .series import SnapCollection

# The units used by HemeLB, which are canonical, in a way.
hlbUnits = {'position' : U.length.mm,
            'pressure' : U.pressure.torr, 
            'shear_strain_rate' : 1 / U.time.second, 
            'velocity' : U.speed.m_per_s, 
            'speed' : U.speed.m_per_s
            }

def CfxSnapshot(filename):
    """Factory function wrapping a CFX snapshot.
    
    Load the data with:
    >>> snap = CfxSnapshot(filename)

    Fields are constructed from the header line.
    
        
    """
    (__raw_row, fieldUnits) = parseHeader(filename)
    __raw_row = [('id', int),] + __raw_row
    fieldUnits['id'] = 1
    
                 # ('position', float, (3,)),
                 # ('strain_rate', float),
                 # ('speed', float),
                 # ('velocity', float, (3,)),
                 # ('wall_shear', float, (4,))]

    __readable_row = np.dtype(__raw_row[1:])
    row = np.dtype(__raw_row)
    
    noindex = np.genfromtxt(filename, skip_header=findStart(filename)+2,
                           delimiter=',',
                           dtype=__readable_row).view(np.recarray)
    index = np.recarray(shape=noindex.shape, dtype=row)
    index.id = np.arange(len(noindex))
    for el in __raw_row[1:]:
        key = el[0]
        index.__setattr__(key, U.convert(noindex.__getattribute__(key), fieldUnits[key], hlbUnits[key]))
        continue
    
    return continuousOrder(index)

def CfxCentreLineSnapshot(filename):
    """Factory function wrapping a CFX snapshot.
    
    Load the data with:
    >>> snap = CfxSnapshot(filename)

    Fields are constructed from the header line.
    
        
    """
    (__raw_row, fieldUnits) = parseHeader(filename, AllData=True)
    __raw_row = [('id', int),] + __raw_row
    fieldUnits['id'] = 1
    
                 # ('position', float, (3,)),
                 # ('strain_rate', float),
                 # ('speed', float),
                 # ('velocity', float, (3,)),
                 # ('wall_shear', float, (4,))]

    __readable_row = np.dtype(__raw_row[1:])
    row = np.dtype(__raw_row)
    
    noindex = np.genfromtxt(filename, skip_header=findStart(filename, AllData=True)+2,
                           delimiter=',',
                           dtype=__readable_row).view(np.recarray)
    index = np.recarray(shape=noindex.shape, dtype=row)
    index.id = np.arange(len(noindex))
    for el in __raw_row[1:]:
        key = el[0]
        index.__setattr__(key, U.convert(noindex.__getattribute__(key), fieldUnits[key], hlbUnits[key]))
        continue
    
    return index

fieldMap = {'velocity': 'speed',
            'velocity_u': ('velocity', 0, 3),
            'velocity_v': ('velocity', 1, 3),
            'velocity_w': ('velocity', 2, 3),
            'x' : ('position', 0, 3),
            'y' : ('position', 1, 3),
            'z' : ('position', 2, 3),
            }

def parseHeader(filename, AllData=False):
    '''Parses header to get data type for record array'''
    fieldRegEx = re.compile('(.+?)\s*\[\s*(.+?)\s*\]')
    
    f = file(filename)
    for i in range(findStart(filename, AllData=AllData)+2):
        header = f.readline()
        continue
    
    fields = []
    fieldUnits = dict()
    for field in header.split(','):
        match = fieldRegEx.match(field.strip())
        field, unit = match.group(1,2)
        field = field.lower().replace(' ', '_')
        fields.append(field)
        # Have to fiddle with the read-in units to be compatible with the enthought function - 
        # it supports "**" for powers and doesn't support spaces in the middle of units.
        # The latter are replaced with "*" characters.
        try:
            niceField = fieldMap[field]
            if isinstance(niceField, tuple): 
                niceField = niceField[0]
                pass
        except KeyError:
            niceField = field
            pass
        fieldUnits[niceField] = U.unit_parser.unit_parser.parse_unit(unit.replace('^', '**').replace(' ', '*'), suppress_unknown = False)
        continue
    
    dtype = []
    for i, field in enumerate(fields):
        try:
            val = fieldMap[field]
        except KeyError:
            val = field
            pass
        if isinstance(val, tuple):
            val, current, total = val
            for j in range(total):
                other = fieldMap[fields[i + j - current]]
                assert val == other[0]
                assert j == other[1]
                assert total == other[2]
                continue
            
            if current == 0:
                dtype.append((val, float, (total,)))
                pass
            
        else:
            dtype.append((val, float))
            pass
        
        continue
    
    return (dtype, fieldUnits)

def findStart(filename, AllData=False):
    """Find the start line."""

    if AllData == True:
        return -1

    for i, line in enumerate(file(filename)):
        if line.find('[Data]')>=0:
            return i
        continue
    
    raise ValueError("File '%s' does not appear to have a '[Data]' section." % self.filename)

class Cfx(SnapCollection):
    """Collection of CFX snapshots, as produced by Savvas."""
    
    loader = staticmethod(CfxSnapshot)
    snapPattern = '*.txt'

    @staticmethod
    def snapToTime(snap, delimiter):
        return float(
            os.path.splitext(os.path.basename(snap))[0].split(delimiter)[1][:-1]
            )

    def orderSnap(self, snap):
        return continuousOrder(snap)
    
    pass

class CfxCentreLine(SnapCollection):
    """Collection of CFX centre line snapshots, as produced by Hywel."""
    
    loader = staticmethod(CfxCentreLineSnapshot)
    snapPattern = '*.txt'

    @staticmethod
    def snapToTime(snap, delimiter):
        return float(
            os.path.splitext(os.path.basename(snap))[0].split(delimiter)[1][:-1]
            )

    pass

class SteadyCfx(Cfx):
    """A slight hack; this class returns the same snapshot for any
    requested time. Intended for use with steady state simulations
    that are time-independent."""
    
    def __init__(self, filename, tol=tol):
        base, self.snapPattern = os.path.split(filename)
        return Cfx.__init__(self, base, tol=tol)
    
    def __getitem__(self, time):
        return Cfx.__getitem__(self, 0.)
    
    @staticmethod
    def snapToTime(snap, delimiter):
        return 0.
    pass
