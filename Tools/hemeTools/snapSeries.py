import os.path
import glob
from .snapshots import HemeLbSnapshot, CfxSnapshot
import numpy as N

tol = 1e-8

class SnapCollection(object):
    """Abstract class for managing a collection of snapshots. Indexed
    by intra-cycle time."""
    
    period_s = 60.0/70.0
    
    def __init__(self, base, tol=tol, delimiter=" "):
        self.base = base
        self.tol = tol
        
        self.snapshots = glob.glob(os.path.join(base, self.snapPattern))
        self.snapmap = {}
        times = []
        for s in self.snapshots:
            t = self.snapToTime(s, delimiter) % self.period_s
            self.snapmap[t] = s
            times.append(t)
            continue
        times = N.array(times)
        times.sort()
        
        self.times = times
        return
 
    def __getitem__(self, time):
        i = self.times.searchsorted(time)
        n = len(self.times)
        
        if i < n and \
                self.times[i]-time <= self.tol:
            t = self.times[i]
        elif i > 0 and \
                self.times[i-1] - time <= self.tol:
            t = self.times[i-1]
        else:
            raise KeyError('No time close enough to "%s" tolerance %s' % (str(time), self.tol))
        
        snap = self.snapmap[t]
        
        return self.loader(snap)
    pass

        
class Heme(SnapCollection):
    """Collection of HemeLB snapshots in original format."""
    
    snapPattern = os.path.join('Snapshots', 'snapshot_*.asc')
    
    def __init__(self, base, timestep_s=None, tol=1e-8):
        if timestep_s is None:
            raise ValueError("Must supply keyword argument 'timestep_s', the length of a timestep in seconds.")
        self.timestep_s = timestep_s
        return SnapCollection.__init__(self, base, tol=tol)
    
    def snapToTime(self, snap, delimiter='_'):
        ts = int(
            os.path.splitext(os.path.basename(snap))[0].split('_')[1]
            )
        return ts *self.timestep_s
    
    def loader(self, snap):
        ans = HemeLbSnapshot(snap)
        ans.computePosition(os.path.join(self.base, 'coords.asc'))
        return ans
    
    pass

class Cfx(SnapCollection):
    """Collection of CFX snapshots, as produced by Savvas."""
    
    loader = staticmethod(CfxSnapshot)
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

class Diff(SnapCollection):
    def __init__(self, sc1, sc2, tolerance=1e-6):
        assert N.allclose(sc1.times, sc2.times, tolerance)
        self.times = sc1.times.copy()
        
        self.sc1 = sc1
        self.sc2 = sc2
        return
    
    def __getitem__(self, time):

        snap1 = self.sc1[time]
        snap2 = self.sc2[time]
        
        fields1 = set((k,v[0]) for k,v in snap1.dtype.fields.iteritems())
        fields2 = set((k,v[0]) for k,v in snap2.dtype.fields.iteritems())
        
        fieldsD = list(set.intersection(fields1, fields2))
        ans = N.recarray(shape=snap1.shape,
                         dtype=fieldsD)
        
        for f in fieldsD:
            name = f[0]
            
            if name in ['id', 'position']:
                assert N.allclose(snap1.__getattribute__(name),
                                  snap2.__getattribute__(name))
                
                ans.__setattr__(name, snap1.__getattribute__(name))
            else:
                ans.__setattr__(name, 
                                snap2.__getattribute__(name) - \
                                    snap1.__getattribute__(name))
                pass
            continue
        
        return ans
    pass

    
if __name__ == "__main__":
    
    pass
