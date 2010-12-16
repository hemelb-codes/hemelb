"""A simplified interface into netCDF to write files suitable for
OpenDX to read.

"""

import netCDF4
import numpy as N
import os.path

def get_or_def_dim(ncf, name, length):
    try:
        dim = ncf.dimensions[name]
    except KeyError:
        dim = ncf.createDimension(name, length)
        dim = ncf.dimensions[name]
    return name

class Positions(object):
    """ABC for positions of data."""
    pass

class Data(object):
    """Wraps data to be written and coordinates its writing."""
    _default_rank=None
    def __init__(self, data, rank=None):
        """Create the data object. Unless you are using one of the
        specialized subclasses (ScalarData/VectorData/TensorData) you
        must specify the rank of the data.

        data -- the data you want to write

        rank -- optional if using one of the specialized subclasses,
        otherwise an integer in [1,2,3]
        
        """
        if rank is None:
            rank = self._default_rank
            pass
        assert rank is not None

        assert data.ndim > rank
        
        self.rank = rank
        self.rankstr = ['', ', vector', ', matrix'][rank]
        
        itemshape = []
        shape = list(data.shape)
        
        for i in range(rank):
            itemshape.insert(0, shape.pop())
            continue
        
        self.shape = tuple(shape)
        self.itemshape = tuple(itemshape)
        self.size = N.multiply.reduce(shape)
        self.ndim = len(shape)
        
        self.data = data
        return

    def define(self, ncf, name):
        """Define the necessary dimensions and variables."""
        ncdims = []
        for i, n in enumerate(self.shape):
            ncdims.append(get_or_def_dim(ncf, name+'_values_n%d' % i, n))
            continue
        
        for i,n in enumerate(self.itemshape):
            ncdims.append(get_or_def_dim(ncf, name+'_values_r%d' % i, n))
            continue
        
        return ncf.createVariable(name+'_values', 'f', ncdims)
    
    def nameString(self, name):
        return name + self.rankstr
    
    def writeto(self, var):
        """Write the data to the netCDF file."""
        var[:] = self.data.astype(N.float32)
        return
    
    pass

class ScalarData(Data):
    _default_rank = 0
    pass
class VectorData(Data):
    _default_rank = 1
    pass
class MatrixData(Data):
    _default_rank = 2
    pass

class IrregularPositions(Positions):
    """Deals with irregular/unconnected positions."""
    def __init__(self, dim, positions):
        """dim -- dimension of the space in which the positions lie

        positions -- the array of positions, should have shape ==
        (nPoints, dim)
        """
        assert dim in [1,2,3]
        self.dim = dim
        
        positions = N.atleast_2d(positions)
        assert positions.ndim == 2
        assert positions.shape[1] == dim

        self.size = positions.shape[0]
        self.shape = (self.size,)
        self.ndim = 1
        
        self.data = self.positions = positions
        return
    
    def define(self, ncf, basename):
        """Define the necessary dimensions and variables."""
        sizeDim = get_or_def_dim(ncf, basename+'_values_n0', self.size)
        dimDim = get_or_def_dim(ncf, basename+'_dim', self.dim)
        
        return ncf.createVariable(self.positionString(basename),
                                  'f', (sizeDim, dimDim))

    @staticmethod
    def positionString(basename):
        return basename+'_locations'

    def writeto(self, var):
        """Write the data to the netCDF file."""
        var[:] = self.data.astype(N.float32)
        return

    pass

class RegularPositions(Positions):
    """Deals with DX regular positions."""
    
    def __init__(self, axes):
        """Initialize the wrapping object.

        axes -- a list/tuple of Axis objects specifying the
        coordinates along each each.
        
        """
        dim = len(axes)
        assert dim in [1,2,3]
        self.dim = dim
        self.axes = axes
        for i, ax in enumerate(axes):
            ax.dim = dim
            ax.aID = i
            continue
        
        return
    
    def define(self, ncf, basename):
        """Define the necessary dimensions and variables."""
        return [a.define(ncf, basename)
                for i, a in enumerate(self.axes)]
        
    def positionString(self, basename):
        ans = ''
        for i in range(self.dim):
            ans += self.axes[i].mkvarname(basename)
            
            if self.dim>1:
                ans += ', product'
                pass
            
            ans += self.axes[i].regular
            if i < (self.dim - 1):
                ans += '; '
            continue
        return ans
    
    def writeto(self, vars):
        """Write the data to the netCDF file."""
        for i in range(self.dim):
            self.axes[i].writeto(vars[i])
        return
    
    pass


class Axis(object):
    """Represent an axis for regular positions."""
    
    def mkvarname(self, basename):
        return basename+'_axis_%d' % self.aID
    
    def define(self, ncf, basename):
        """Define the necessary dimensions and variables."""
        dimDim = get_or_def_dim(ncf, basename+'_naxes', self.dim)
        sizeDim = get_or_def_dim(ncf, self.mkdimname(basename), self.n)
        return ncf.defineVariable(self.mkvarname(basename), 'f', (sizeDim, dimDim))
    
    pass

class RegularAxis(Axis):
    """Regular axes have an origin and a delta. Positions are given
    by: x[i] = origin + i*delta

    """
    regular = ', compact'
    
    def __init__(self, origin=0., delta=1., dim=None, aID=None):
        """origin - starting point for coordinates

        delta - amount to move along axis per point

        The following 2 args probably shouldn't be used outwith this
        file.
        
        dim - (optional) which axis this is; if omitted will be calculated by the using field

        aID - (optional) a unique ID for the axis; if omitted will be calculated by the using field

        """
        self.origin = origin
        self.delta = delta
        self.n = 2
        self.dim = dim
        self.aID = aID
        return
    
    def mkdimname(self, basename):
        return basename+'compact_dim'
    
    
    def writeto(self, var):
        """Write the data to the netCDF file."""
        pos = N.zeros((self.n, self.dim), dtype=N.float32)
        pos[0, self.aID] = self.origin
        pos[1, self.aID] = self.delta
        var[:] = pos
        return
    pass

class IrregularAxis(Axis):
    """Irregular axes have a monotonically increasing set of values,
    passed to the constructor.
    
    """
    regular = ''
    def __init__(self, points, dim=None, aID=None):
        """points - the points along the axis
        
        The following 2 args probably shouldn't be used outwith this
        file.
        
        dim - (optional) which axis this is; if omitted will be
        calculated by the using field

        aID - (optional) a unique ID for the axis; if omitted will be
        calculated by the using field
        
        """
        # has to be floats not doubles
        points = N.array(points, dtype=N.float32).squeeze()
        points = N.atleast_1d(points)
        
        assert points.ndim == 2
        assert len(points) >= 1
        
        self.n = points.shape[0]
        self.dim = points.shape[1]
        self.points = points
        self.aID = aID
        return
    
    def mkdimname(self, basename):
        return self.mkvarname(basename)+'_len'
    
    def writeto(self, var):
        """Write the data to the netCDF file."""
        var[:] = self.points
        return
    
    pass

class Field(object):
    """Represents an OpenDX field.

    Coordinates writing of its dependent positions and dimensions.
    
    """
    
    def __init__(self, posns, **kwargs):
        """Initialize but do not write the field.

        posns - the positions for the data in the field.
        
        name - keyword arg giving the name of the field, default is
        'data'

        Further kwargs are assumed to be Data objects with the name
        given by the keyword's name.
        
        """
        self.name = kwargs.pop('name', 'data')
        self.positions = posns
        
        assert self.checkdict(kwargs, lambda key, val: isinstance(val, Data)),\
               "All keyword args must be Data istances"
        
        self.dataDict = kwargs
        return
    
    def writeto(self, ncf):
        """Write the field to the netCDF file, coordinating its
        children's definition and writing.

        """
        #ncf = CDF(filename, mode=CREATE)
        
        #ncf.definemode()
        pos = self.positions.define(ncf, self.name)

        datavars = {}
        for dname, d in self.dataDict.iteritems():
            fullvarname = "%s_%s" % (self.name, dname)
            dat = d.define(ncf, fullvarname)
            dat.field = d.nameString(fullvarname)
            dat.positions = self.positions.positionString(self.name)
            datavars[dname] = dat
            continue
        
        #ncf.enddef()

        self.positions.writeto(pos)
        
        [d.writeto(datavars[dname]) for dname, d in self.dataDict.iteritems()]
        
        return
    
    @classmethod
    def checkdict(cls, d, condition):
        """Checks if a condition is satisfied by all values of a
        dictionary. Returns True if so, False otherwise.

        d -- the dictionary

        condition -- a callable accepting 2 arguments: the key and the
        value
        
        """
        for key, val in d.iteritems():
            if condition(key, val):
                # we're ok here
                pass
            else:
                # test fails
                return False
            continue
        
        return True
    
    pass

class Unconnected(Field):
    def __init__(self, posns, **kwargs):
        assert isinstance(posns, IrregularPositions)
        
        Field.__init__(self, posns, **kwargs)
        
        assert self.checkdict(self.dataDict,
                              lambda key,val: val.ndim == 1), \
                              "data must be a list of positions"
        
        assert self.checkdict(self.dataDict,
                              lambda key,val: val.size == posns.size), \
                              "data & posns must have same size"
        
        return
    
    pass

class Connected(Field):
    def __init__(self, posns, **kwargs):
        assert isinstance(posns, RegularPositions)
        
        Field.__init__(self, posns, **kwargs)
        return


    pass

# class LoadedData(object):
#     def __init__(self, dname, vname, v):
#         self.vname = vname
#         self.dname = dname
#         self.v = v
        
# class LoadedField(object):
#     def __init__(self, name):
#         self.name = name
#         self.positions = None
#         self.dataDict = {}
#         return
    
#     def __getattribute__(self, name):
#         try:
#             return object.__getattribute__(self, name)
#         except AttributeError:
#             dd = object.__getattribute__(self,'dataDict')
#             if name in dd:
#                 return dd[name]
#             raise
#         return
    
#     pass

# class LoadedFields(object):
#     def __init__(self, ncfName):
#         self.ncf = CDF(ncfName)
#         vars = self.ncf.variables()
#         dims = self.ncf.dimensions()
#         self.fields = {}
#         for vname, vinfo in vars.iteritems():
#             v = self.ncf.var(vname)
#             attr = v.attributes()
#             if 'field' in attr:
#                 # it's a data field
#                 splitFieldAttr = attr['field'].split(',')
#                 try:
#                     rankStr = splitFieldAttr[1].strip()
#                 except IndexError:
#                     rankStr = 'scalar'
                    
#                 fname, dname  = splitFieldAttr[0].split('_')
#                 try:
#                     f = self.fields[fname]
#                 except KeyError:
#                     f = self.fields[fname] = LoadedField(fname)
#                     pass
                
#                 f.dataDict[dname] = Data(v.get(),
#                                          rank={'scalar': 0,
#                                                'vector': 1,
#                                                'matrix': 2}[rankStr])
#                 if f.positions is not None:
#                     continue
                
#                 posdesc = attr['positions']
#                 if 'product' in posdesc:
#                     # connected
#                     axes = []
#                     for dimdesc in posdesc.split(';'):
#                         parts = dimdesc.split(',')
#                         aID = parts[0].split('_')[2]
#                         dimvar = self.ncf.var(parts[0].strip())
#                         dimdat = dimvar.get()
#                         if 'compact' in parts:
#                             origin = dimdat[0,aID]
#                             delta = dimdat[1,aID]
#                             axes.append(RegularAxis(origin,delta))
#                         else:
#                             axes.append(IrregularAxis(dimdat))
#                             pass
#                         continue
#                     pos = RegularPositions(axes)
#                 else:
#                     # unconnected
#                     posvar = self.ncf.var(posdesc.strip())
#                     pos = IrregularPositions(dims[posvar.dimensions()[1]],
#                                              posvar.get())
#                     pass
#                 f.positions = pos
#                 pass
#             continue
#         return
    
#     def __getattribute__(self, name):
#         try:
#             return object.__getattribute__(self, name)
#         except AttributeError:
#             f = object.__getattribute__(self,'fields')
#             if name in f:
#                 return f[name]
#             raise
#         return

#     pass

class Series(object):
    """A series of fields connected through a series variable
    (e.g. time). You should ensure that you write the fields such that
    the series variable increases monotonically. The class will not
    however enforce this but attempting to read from an instance for
    which this has not been obeyed will have undefined results!
    
    """
    def __init__(self, dir, **kwargs):
        """dir - the directory containing the results

        mode - single character string for which mode to use for
        opening: 'r' = read; 'w' = write, clobbering any existing
        data; 'a' = append only
        
        seriesvar - name of the series variable (default is 'time');
        only meaningful in 'w' mode

        template - standard python string interpolation template for
        sub-file's names. Names will be created by (template % var) +
        '.ncdx' so these must have enough precision to be unique for
        the values you intend to use.

        """
        
        try:
            func = {
            'r': self.openForRead,
            'w': self.create,
            'a': self.openForWrite
            }[kwargs.pop('mode', 'r')]
        except KeyError:
            raise ValueError('Invalid mode specified')
        
        func(dir, **kwargs)
        return
    
    def openForRead(self, dir, **kwargs):
        """Open for reading. Use the constructor with mode='r'
        instead.

        """
        
        if not os.path.exists(dir):
            raise IOError('Series directory does not exist: %s' % dir)
        self.load(dir)
        return

    def openForWrite(self, dir, seriesvar='time', template='%.9d'):
        """Open a (possibly existing) series for writing. Use the
        constructor with mode='a' instead.

        """
        if not os.path.exists(dir):
            self.create(dir, seriesvar, template)
        else:
            self.load(dir)
        return
    
    def load(self, dir):
        """Open an existing series for writing. Use the
        constructor with mode='a' instead.

        """
        self.dir = dir
        self.indexfile = os.path.join(dir, 'index')
        if not os.path.exists(self.indexfile):
            raise IOError('Series index does not exist: %s' % self.indexfile)
        
        header = file(self.indexfile).readline()
        if header[0] == '#':
            _, self.seriesvar = header[1:].split()
        else:
            self.seriesvar = 'time'
            pass
        self.index = [[el[0], eval(el[1])] for el in N.loadtxt(self.indexfile, object)]
        return
    
    def create(self, dir, seriesvar='time', template='%.9d'):
        """Open for writing, replacing any existing data. Use the
        constructor with mode='w' instead.

        """
        
        if os.path.exists(dir):
            #clobber it
            import shutil
            shutil.rmtree(dir)
            pass
        
        self.dir = dir
        os.mkdir(dir)
        self.indexfile = os.path.join(dir, 'index')
        file(self.indexfile, 'w').write('# path %s\n' % seriesvar)
        self.seriesvar = seriesvar
        self.template = template
        self.index = []
        return
    
    def nextFile(self, var):
        """Open the next file for the series variable = var.

        """
        filebase = (self.template % var) + '.ncdx'
        filename = os.path.join(self.dir, filebase)

        self.index.append([filebase, var])
        file(self.indexfile, 'a').write('%s %s\n' % (filebase, repr(var)))
        ncf = netCDF4.Dataset(filename, 'w')
        return ncf
    
    def append(self, var, fields):
        """Append a new field.

        var - value of the series variable, should be greater than all
        existing

        fields - iterable of the fields to add the series at this time
        
        """
        ncf = self.nextFile(var)
        for f in fields:
            f.writeto(ncf)
            continue
        return

    def __len__(self):
        return len(self.index)
    
    def __getitem__(self, key):
        i = bisect_left(self.index, key)
        if self.index[i][1] == key:
            return LoadedFields(os.path.join(self.dir, self.index[i][0]))
        else:
            raise IndexError('Key "%s" not in index of "%s"' % (repr(key), self.dir))
        return
    
    def __iter__(self):
        for k in self.keys():
            yield self[k]
            
    def __contains__(self, key):
        i = bisect_left(self.index, key)
        if self.index[i][1] == key:
            return True
        else:
            return False
        return
    
    def keys(self):
        """Return a list of all the values of the series variable."""
        return [el[1] for el in self.index]
    pass

def bisect_left(a, x):
    lo = 0
    hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if a[mid][1] < x: lo = mid+1
        else: hi = mid
    return lo
def bisect_right(a, x):
    hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if x < a[mid][1]: hi = mid
        else: lo = mid+1
    return lo

if __name__ == "__main__":
    import pdb
    nPoints = 100
    pos = N.random.uniform(size=(nPoints,3))
    val = N.random.normal(size=nPoints)
    uField = Unconnected(IrregularPositions(3, pos), values=ScalarData(val))
    
    cdf = netCDF4.Dataset('test.nc', 'w', format='NETCDF4')
    uField.writeto(cdf)
    cdf.close()
