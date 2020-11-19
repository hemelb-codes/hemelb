# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from ..interval import Interval
from ..iohelp import LoadableMixin

class Integrator:
    def __init__(self, A):
        self.A = A
        self.Id = sp.identity(np.shape(A)[0])
    pass

class ImplicitEuler(Integrator):
    def __call__(self, y0, iv):
        
        for ignored in range(iv.n):
            y0 = spsolve(self.Id - iv.dt * self.A, y0)
        return y0

    pass

class Trapezoidal(Integrator):
    def __call__(self, y0, iv):
        for ignored in range(iv.n):
            b = np.copy(y0)
            b *= self.Id + 0.5 * iv.dt * self.A
            y0 = spsolve(self.Id - 0.5 * iv.dt * self.A, b)
        return y0
    

class Problem(LoadableMixin):
    '''Solves the heat equation:
    u_t = nu*u_xx
    on [0,1] with boundary values u(0) = u(1) = 0.0.
    
    '''

    def __init__(self, params, ic, iv, solver):
        '''
        params - the dict of parameters, keys being:
            * nx - number of total finite difference nodes 
            * nu - diffusivity

        ic - initial conditions - can be:
            * a callable f(x)
            * an array of appropriate shape and type
            * a scalar

        iv - time interval
        '''
        
        self.nx = params['nx']
        self.nu = params['nu']
        
        self.full_xs = np.linspace(0.0, 1.0, self.nx)
        self.dx = self.full_xs[1]
        self.xs = self.full_xs[1:-1]
        
        # Initial value is sin(pi*x)
        try:
            # First try calling it as a f(x) initial condition
            self.y0 = ic(self.xs)
        except TypeError:
            if instance(ic, np.ndarray):
                # Try assuming it's an appropriate ndarray-like
                assert ic.shape == self.xs.shape, \
                  'shape of initial condition data does not match x shape'
                assert y.dtype == float, 'wrong dtype'
                self.y0 = ic.copy()
            else:
                # maybe a scalar
                self.y0 = np.empty(self.xs.shape)
                self.y0[:] = ic
        
        self.iv = iv
        
        IntCls = {'trap': Trapezoidal, 'euler': ImplicitEuler}[solver.lower()]
        self.integrator = IntCls(self._fd_matrix())
        return
    
    def _fd_matrix(self):
        stencil = [1.0, -2.0, 1.0]
        N = self.nx - 2
        A = sp.diags(stencil, [-1, 0, 1],
                         shape=(N,N), format='csc')
        A *= self.nu / (self.dx**2)
        return A

    def run(self):
        return self.integrator(self.y0, self.iv)

    @classmethod
    def from_dict(cls, state):
        iv = Interval.from_dict(state['time'])
        ic = state['initial']
        if isinstance(ic, str):
            # Probably should be treated as a function
            # Assume it's in np context
            ic = eval('lambda x:' + ic, np.__dict__)

        return cls(
            state['params'],
            ic,
            iv,
            state['solver'])
    pass

