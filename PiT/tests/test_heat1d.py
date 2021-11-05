# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import pytest
import io
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from parareal.interfaces.heat1d import Problem, ImplicitEuler
from parareal.interval import Interval


@pytest.fixture(scope='module')
def daniels_results():
    def getA(N, dx, nu):
        stencil = [1.0, -2.0, 1.0]
        A = sp.diags(stencil, [-1, 0, 1], shape=(N,N), format='csc')
        A *= nu/(dx**2.0)
        return A

    def restrict(y):
        return y[1::2] # use injection, that is remove every second data point

    def implicitEuler(nsteps, dt, A, y0):
        Id = sp.identity(np.shape(A)[0])
        for i in range(nsteps):
            b = np.copy(y0)
            y0 = spsolve(Id - dt*A, b)
        return y0

    def trapezoidal(nsteps, dt, A, y0):
        Id = sp.identity(np.shape(A)[0])
        for i in range(nsteps):
            b = np.copy(y0)
            b *= (Id + 0.5*dt*A)
            y0 = spsolve(Id - 0.5*dt*A, b)
        return y0
    # parameter
    nu     = 0.1 # diffusivity
    N      = 1001 # number of total finite difference nodes on fine mesh (boundary nodes not included).. needs to be odd to work properly!
    nsteps = N   # number of time steps
    Tend   = 1.1 # final time
    dt    = Tend/float(nsteps) # time step length

    # Solution at endpoints is assumed to be zero, so create a mesh without 0 and 1
    mesh = np.linspace(0.0, 1.0, N)
    mesh = mesh[1:N-1]

    # create coarse mesh by removing every second point
    mesh_coarse = mesh[1::2]
    N_coarse = np.size(mesh_coarse)

    # Initial value is sin(pi*x)
    y0   = np.sin(np.pi*mesh)

    # get finite difference matrix on fine and coarse mesh
    A        = getA(np.size(mesh), mesh[1]-mesh[0], nu)
    A_coarse = getA(np.size(mesh_coarse), mesh_coarse[1]-mesh_coarse[0], nu)

    # Run both integrators on fine level
    y_ie = implicitEuler(nsteps, dt, A, y0)
    y_tp = trapezoidal(nsteps, dt, A, y0)

    # ...and run both integrators on coarse level
    y_ie_coarse = implicitEuler(nsteps, dt, A_coarse, restrict(y0))
    y_tp_coarse = trapezoidal(nsteps, dt, A_coarse, restrict(y0))

    # exact solution is sin(pi*x)*exp(-nu*pi^2*t)
    yex  = y0*np.exp(-nu*np.pi**2.0*Tend)
    return {
        'exact': yex,
        'euler': {
            'fine': y_ie,
            'coarse': y_ie_coarse
            },
        'trap': {
            'fine': y_tp,
            'coarse': y_tp_coarse
            }
        }

params = {
    'fine': {'nu': 0.1, 'nx': 1001},
    'coarse': {'nu': 0.1, 'nx': 501}
    }

def initial_cond(x):
    return np.sin(np.pi*x)

def time_interval():
    nsteps = 1001 # number of time steps
    Tend = 1.1 # final time
    dt = Tend/float(nsteps) # time step length    
    return Interval(dt, 0, nsteps)
time_interval = time_interval()

@pytest.mark.parametrize('solver', ('euler', 'trap'))
@pytest.mark.parametrize('resolution', ('coarse', 'fine'))
def test_vs_daniel(daniels_results, solver, resolution):
    p = Problem(params[resolution], initial_cond, time_interval, solver)
    y = p.run()
    assert np.all(y == daniels_results[solver][resolution])
    
probspec = u'''
time:
  start: 0.0
  dt: 1.1 / n
  n: 1001
params:
  nu: 0.1
  nx: 501
initial:
  sin(pi*x)
solver:
  euler
'''
def test_from_file():
    p = Problem.from_stream(io.StringIO(probspec))
    assert p.nx == 501
    assert p.nu == 0.1
    expected = np.sin(np.pi*np.linspace(0.0, 1.0, p.nx)[1:-1])
    assert np.all(p.y0 == expected)
    assert isinstance(p.integrator, ImplicitEuler)
