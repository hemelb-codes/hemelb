# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import numpy as np
from .iohelp import *

class Interval(LoadableMixin, DumpableMixin):
    """Simple time interval representation.
    """
    def __init__(self, dt, start, n):
        '''
        dt - time step
        start - real start time
        n - number of steps to take
        '''

        self.dt = dt
        self.start = start
        self.n = n
        return

    def generate(self):
        return self.dt * (np.arange(self.n+1) + self.start)
    
    def subdivide(self, num_time_slices):
        if self.n % num_time_slices:
            raise ValueError("time interval cannot be subdivided exactly")
        
        slice_steps = self.n / num_time_slices
        return [Interval(self.dt, self.start + slice_steps*i, slice_steps) for i in range(num_time_slices)]

    def to_dict(self):
        return {
            'start': self.start,
            'n': self.n,
            'dt': self.dt
            }

    @classmethod
    def from_dict(cls, rawd):
        d = {}
        for k, v in rawd.items():
            if isinstance(v, str):
                ctx = rawd.copy()
                ctx.update(d)
                d[k] = eval(v, np.__dict__, ctx)
            else:
                d[k] = float(v)
        return cls(d['dt'], d['start'], d['n'])
    
    pass
