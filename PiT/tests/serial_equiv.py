# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import numpy as np
import matplotlib.pyplot as plt
import py.path
import luigi

from test_hm import test_44_equiv_serial
import interval


def analytic(t):
    ''' x(t) = (A sin wt + B cos wt) exp(-bt/2)

    where (to satisfy x(0) = 1, v(0) = 0):
        A = b/(2w)
        B = 1
    
    '''
    beta = 0.1
    omega0 = 1.0
    
    omega = 0.5 * np.sqrt(4.0 * omega0**2 - beta**2)
    
    A = beta / (2.0*omega)
    sin_wt = np.sin(omega*t)
    cos_wt = np.cos(omega*t)
    e_bt2 = np.exp(-0.5 * beta * t)
    x = (A * sin_wt + cos_wt) * e_bt2
    v = -0.5 * beta * x + (A * omega * cos_wt - omega * sin_wt) * e_bt2
    return x, v

class Run:
    def __init__(self, data):
        self.data = data
        self.t, self.x, self.v = data.transpose()

if __name__ == '__main__':
    rundir = py.path.local('serial_equivalence')
    if rundir.exists():
        rundir.remove()
        pass
    rundir.mkdir()

    test_44_equiv_serial(rundir, luigi, plot_mode=True)

    with rundir.as_cwd():
        serial = Run(np.loadtxt('serial.txt'))
        parareal = Run(np.loadtxt('parareal.txt'))

    time = interval.Interval(2.0 * np.pi/1000,0.0, 1000)
    sub_times = time.subdivide(4)
    
    assert np.all(parareal.t == serial.t)

    t = parareal.t
    exact_x, exact_v = analytic(t)
    delta_x = parareal.x - serial.x
    delta_v = parareal.v - serial.v

    class Plotter(object):
        def __init__(self, x):
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111)
            self.x = x
            
            self.ylim = np.array([np.inf, -np.inf])
            
        def _update(self, y):
            self.ylim[0] = min(self.ylim[0], np.min(y))
            self.ylim[1] = max(self.ylim[1], np.max(y))
            
        def plot(self, y, **kwargs):
            self.ax.plot(self.x, y, **kwargs)
            self._update(y)
        def scatter(self, y, **kwargs):
            self.ax.scatter(self.x, y, **kwargs)
            self._update(y)
            
        def show_or_save(self, fname=None):
            pows = np.log10(np.abs(self.ylim))
            pow10 = np.floor(pows.max())
            base = 10**pow10
            ylim = self.ylim / base
            ymin = np.floor(ylim[0]) * base
            ymax = np.ceil(ylim[1]) * base

            self.ax.set_ylim(ymin, ymax)
            
            if fname is None:
                self.fig.show()
            else:
                self.fig.savefig(fname)
    
    def plt_x(fname=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t, exact_x)
        ax.scatter(t, serial.x)
        ax.scatter(t, parareal.x)
        show_or_save(fig, fname)

    def plt_v(fname=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t, exact_v)
        ax.scatter(t, serial.v)
        ax.scatter(t, parareal.v)
        show_or_save(fig, fname)

    def plt_x_err(fname=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(t, serial.x-exact_x)
        ax.scatter(t, parareal.x-exact_x)
        show_or_save(fig, fname)

    def plt_v_err(fname=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(t, serial.v-exact_v)
        ax.scatter(t, parareal.v-exact_v)
        show_or_save(fig, fname)

    def plot2(func, fname=None):
        p = Plotter(t)
        p.scatter(func(serial))
        p.scatter(func(parareal))
        p.show_or_save(fname)
        
    def plt_err(fname=None):
        def err(run):
            return 
        p = Plotter(t)
        
        p.scatter(err(serial))
        p.scatter(err(parareal))
        p.show_or_save(fname)
        
    plot2(lambda run: np.sqrt((run.x - exact_x)**2 + (run.v - exact_v)**2),
              'rms_err.pdf')
    
# #aplot = ax.plot(t, exact_x)
# splot = ax.scatter(t, serial.x-exact_x)
# pplot = ax.scatter(t, parareal.x-exact_x)
# #dplot = ax.scatter(t, delta_x)
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111)

# #ax.plot(t, exact_v)
# ax.scatter(t, serial.v-exact_v)
# ax.scatter(t, parareal.v-exact_v)
# #dplot = ax.scatter(t, delta_x)
# plt.show()
