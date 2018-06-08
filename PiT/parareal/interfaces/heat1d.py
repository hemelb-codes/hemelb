import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp1d

from ..interval import Interval

import pdb

class Integrator:
    def __init__(self, A):
        self.A = A
        self.Id = sp.identity(np.shape(A)[0])
    pass

class ImplicitEuler(Integrator):
    def __call__(self, y0, iv):
        
        for ignored in range(iv.n):
            # I think I can skip this?
            b = np.copy(y0)
            y0 = spsolve(self.Id - iv.dt * self.A, b)
        return y0

    pass

class Trapezoidal(Integrator):
    def __call__(self, y0, iv):
        for ignored in range(iv.n):
            b = np.copy(y0)
            b *= self.Id + 0.5 * iv.dt * self.A
            y0 = spsolve(self.Id - 0.5 * iv.dt * self.A, b)
        return y0
    

class Problem:
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
        
def interpolate(x_fine, x_coarse, fx_coarse):
    # To interpolate, need to add boundary values manually
    xc_bc  = np.insert(x_coarse, 0, 0.0)
    xc_bc  = np.append(xc_bc, 1.0)
    fxc_bc = np.insert(fx_coarse, 0, 0.0)
    fxc_bc = np.append(fxc_bc, 0.0)
    f = interp1d(xc_bc, fxc_bc, kind='linear')
    return f(x_fine)

if __name__ == '__main__':
    # Daniel's test problem
    
    nu = 0.1 # diffusivity
    nsteps = 1001 # number of time steps
    Tend = 1.1 # final time
    dt = Tend/float(nsteps) # time step length

    fine_params = {'nu': nu, 'nx': 1001}
    initial = lambda x:  np.sin(np.pi*x)
    time = Interval(dt, 0, nsteps)
    
    pf = Problem(fine_params, initial, time, 'Euler')
    y_ie = pf.run()
    pf = Problem(fine_params, initial, time, 'trap')
    y_tp = pf.run()
    
    # exact solution is sin(pi*x)*exp(-nu*pi^2*t)
    yex  = pf.y0 * np.exp(-nu * np.pi**2 * Tend)

    print ("==== Fine level ====")
    print ("Implicit Euler error:   %5.3e" % np.linalg.norm(yex - y_ie))
    print ("Trapezoidal rule error: %5.3e" % np.linalg.norm(yex - y_tp))

    coarse_params = {'nu': nu, 'nx': 501}

    pc = Problem(coarse_params, initial, time, 'Euler')
    y_ie_coarse = pc.run()
    pc = Problem(coarse_params, initial, time, 'trap')
    y_tp_coarse = pc.run()

    # interpolate coarse level solutions back to fine level
    y_ie_interp = interpolate(pf.xs, pc.xs, y_ie_coarse)
    y_tp_interp = interpolate(pf.xs, pc.xs, y_tp_coarse)

    assert pf.xs[1::2].shape == pc.xs.shape
    assert np.all(pc.xs == pf.xs[1::2])
    
    print ("==== Coarse level ====")
    print ("Implicit Euler error:   %5.3e" % np.linalg.norm(yex - y_ie_interp))
    print ("Trapezoidal rule error: %5.3e" % np.linalg.norm(yex - y_tp_interp))
