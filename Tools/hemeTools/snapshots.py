import numpy as N
import re
import xdrlib
from enthought import units as U

def HemeLbSnapshot(filename):
    """Guess which file format we were given and use the correct class
    to open it.
    """
    try:
        reader = xdrlib.Unpacker(file(filename).read(4))
        stable = reader.unpack_int()
        if stable == 0 or stable == 1:
            cls = XdrHemeLbSnapshot
        else:
            cls = TextHemeLbSnapshot
            pass
        
    except:
        cls = TextHemeLbSnapshot
        pass
    
    return cls(filename)

class BaseHemeLbSnapshot(N.recarray):
    """Base class wrapping a HemeLB snapshot.
    
    Snap is basically a numpy record array with the following fields:

        - id (int) -- an id number (basically the index of the point in the
            file

        - position (3x float) -- the position in input (STL) space
            (typically mm)

        - grid (3x int) -- the (x, y, z) coordinates in lattice units

        - pressure (float) -- the pressure in physical units (mmHg)

        - velocity (3x float) -- (x,y,z) components of the velocity field
            in physical units (m/s)

        - stress (float) -- the von Mises stress in physical units (Pa)

    Additionally it has the following properties:

    stable, voxel_size, bb_min, bb_max, bb_len, voxel_count 

    (see __readHeader for full details)

    """
    _raw_row = [('id', int),
                 ('position', float, (3,)),
                 ('grid', int, (3,)),
                 ('pressure', float),
                 ('velocity', float, (3,)),
                 ('stress', float)]
    _readable_row = N.dtype(_raw_row[2:])
    row = N.dtype(_raw_row)

    _attrs = ('stable', 'voxel_size', 'bb_min', 'bb_max', 'bb_len', 'voxel_count')
    header = len(_attrs)


    def __new__(cls, filename):
        """Create a new instance. Numpy array subclasses use this
        method instead of __init__ for initialization.

        """
        noindex = cls._load(filename)
        index = N.recarray(shape=noindex.shape, dtype=cls.row)
        index.id = N.arange(len(noindex))
        index.position = N.nan  # Don't know the coordinate file name
            # at this moment

        for el in cls._raw_row[2:]:
            key = el[0]
            index.__setattr__(key, noindex.__getattribute__(key))
            continue

        obj = index.view(cls)

        obj.stable, obj.voxel_size, obj.bb_min, \
            obj.bb_max, obj.bb_len, obj.voxel_count = cls._readHeader(filename)

        return obj

    def __array_finalize__(self, parent):
        """Numpy special method."""

        if parent is None:
            return
        for a in self._attrs:
            setattr(self, a, getattr(parent, a, None))
            continue

        return

    def computePosition(self, coordsFile):
        """Given the coordinate file from the segtool, calculate all
        the lattice positions' coordinates.

        """
        from .coordinates import Transformer
        trans = Transformer(coordsFile)

        self.position = trans.siteToStl(self.grid + self.bb_min)
        return

    pass

class TextHemeLbSnapshot(BaseHemeLbSnapshot):
    """Read a text snapshot.
    """
    @classmethod
    def _readHeader(cls, filename):
        """Read the header lines, according to:
        0- Flag for simulation stability, 0 or 1
        1- Voxel size in physical units (units of m)
        2- vertex coords of the minimum bounding box with minimum values (x, y and z values)
        3- vertex coords of the minimum bounding box with maximum values (x, y and z values)
        4- #voxels within the minimum bounding box along the x, y, z axes (3 values)
        5- total number of fluid voxels

        """

        f = file(filename)
        stable = int(f.readline())
        voxel_size = float(f.readline())
        bb_min = N.array([int(x) for x in f.readline().split()])
        bb_max = N.array([int(x) for x in f.readline().split()])
        bb_len = N.array([int(x) for x in f.readline().split()])
        voxel_count = int(f.readline())

        return stable, voxel_size, bb_min, bb_max, bb_len, voxel_count

    @classmethod
    def _load(cls, filename):
        return N.loadtxt(filename,
                         skiprows=cls.header,
                         dtype=cls._readable_row).view(N.recarray)

    pass

class XdrHemeLbSnapshot(BaseHemeLbSnapshot):
    """Read an XDR snapshot.
    """
    # int float 3x int 3x int 3x int int
    _headerLengthBytes = 4 + 8 + 3*4 + 3*4 + 3*4 + 4

    @classmethod
    def _readHeader(cls, filename):
        """Read the header lines, according to:
        0- Flag for simulation stability, 0 or 1
        1- Voxel size in physical units (units of m)
        2- vertex coords of the minimum bounding box with minimum values (x, y and z values)
        3- vertex coords of the minimum bounding box with maximum values (x, y and z values)
        4- #voxels within the minimum bounding box along the x, y, z axes (3 values)
        5- total number of fluid voxels

        """
        try:
            assert cls._headerCache[0] == filename
            header = cls._headerCache[1]
        except:
            reader = xdrlib.Unpacker(file(filename).read(_headerLengthBytes))
            header = cls._actuallyReadHeader(cls, reader)
            pass
        
        return (header['stable'], header['voxel_size'],
                header['bb_min'], header['bb_max'],
                header['bb_len'], header['voxel_count'])
    
    @classmethod
    def _actuallyReadHeader(cls, reader):
        header = {}
        header['stable'] = reader.unpack_int()
        header['voxel_size'] = reader.unpack_double()

        header['bb_min'] = N.array((reader.unpack_int(),
                                    reader.unpack_int(),
                                    reader.unpack_int()))
        header['bb_max'] = N.array((reader.unpack_int(),
                                    reader.unpack_int(),
                                    reader.unpack_int()))
        header['bb_len'] = N.array((reader.unpack_int(),
                                    reader.unpack_int(),
                                    reader.unpack_int()))
        header['voxel_count'] = reader.unpack_int();
        return header

    @classmethod
    def _load(cls, filename):
        reader = xdrlib.Unpacker(file(filename).read())
        header = cls._actuallyReadHeader(reader)

        ans = N.recarray((header['voxel_count'],), dtype=cls._readable_row)

        for i in xrange(header['voxel_count']):
            ans[i] = ((reader.unpack_int(),
                      reader.unpack_int(),
                      reader.unpack_int()),
                      reader.unpack_float(),
                      (reader.unpack_float(),
                      reader.unpack_float(),
                      reader.unpack_float()),
                      reader.unpack_float())

            continue
        reader.done()
        cls._headerCache = (filename, header)

        return ans
    pass

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

    __readable_row = N.dtype(__raw_row[1:])
    row = N.dtype(__raw_row)

    noindex = N.genfromtxt(filename, skip_header=findStart(filename)+2,
                           delimiter=',',
                           dtype=__readable_row).view(N.recarray)
    index = N.recarray(shape=noindex.shape, dtype=row)
    index.id = N.arange(len(noindex))
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

def parseHeader(filename):
    '''Parses header to get data type for record array'''
    fieldRegEx = re.compile('(.+?)\s*\[\s*(.+?)\s*\]')

    f = file(filename)
    for i in range(findStart(filename)+2):
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

def findStart(filename):
    """Find the start line."""
    for i, line in enumerate(file(filename)):
        if line.find('[Data]')>=0:
            return i
        continue

    raise ValueError("File '%s' does not appear to have a '[Data]' section." % self.filename)
