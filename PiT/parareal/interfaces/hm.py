# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# Solver for damped harmonic motion

import numpy as np
import scipy.integrate
import yaml

from ..interval import Interval
from ..iohelp import LoadableMixin

class Integrator(object):
    """Our system is:
    d^2x / dt^2 + b dx/dt + omega^2 x = 0

    Transform to 2 coupled equations:

    0: dx/dt = v
    1: dv/dt = -b v - omega^2 x

    Define y = (x, v)
    """
    
    def __init__(self, params):
        self.b = params['b']
        self.omega = params['omega']
        # The system coefficients
        self.mat = np.array([
            [0.0, 1.0],
            [-self.omega**2, -self.b]], dtype=float)

        return
        
    def rhs(self, y, t):
        """Evaluates the rhs of the system."""
        
        return np.matmul(self.mat, y)

    def jac(self, y, t):
        return self.mat

    pass

class CoarseInt(Integrator):
    def __call__(self, y0, iv):
        ts = iv.generate()
        vals, info = scipy.integrate.odeint(
            self.rhs,
            y0,
            ts,
            rtol=1e-2,
            atol=1e-3,
            full_output=True)

        return (ts, vals)
    pass
class FineInt(Integrator):
    def __call__(self, y0, iv):
        ts = iv.generate()
        vals, info = scipy.integrate.odeint(
            self.rhs,
            y0,
            ts,
            Dfun=self.jac,
            rtol=1e-4,
            atol=1e-6,
            h0=iv.dt,
            full_output=True)
        return (ts, vals)
    pass

class Problem(LoadableMixin):
    def __init__(self, params, ic, iv, solver):
        self.yi = np.array((ic['x'], ic['v']), dtype=float)
        self.iv = iv

        IntCls = {'fine': FineInt, 'coarse': CoarseInt}[solver]
        self.integrator = IntCls(params)
        return
        
    @classmethod
    def from_dict(cls, state):
        params = state['params']
        time = state['time']
        iv = Interval.from_dict(time)
        return cls(
            state['params'],
            state['initial'], iv,
            state['solver'])
    
    def run(self, outfile):
        ts, vals = self.integrator(self.yi, self.iv)
        data = np.concatenate((ts[:, np.newaxis], vals), axis=1)
        np.savetxt(outfile, data)
        
if __name__ == '__main__':
    import sys
    prob = Problem.from_file(sys.argv[1])
    prob.run(sys.argv[2])
    
    
