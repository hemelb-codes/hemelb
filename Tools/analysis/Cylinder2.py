import numpy as np
import quantities as pq
import os.path
import Gnuplot
import warnings

from generation.PoiseuilleFlowCylinderGeneration import PoiseuilleFlowCylinderGeneration

from .memo import memo_property
from .cache import cache_property
from .HemeLbRunResults import HemeLbRunResults

import pdb
        
class CylinderResults(HemeLbRunResults, PoiseuilleFlowCylinderGeneration):
    RunTypeName = 'Cylinder'    
    ConvergenceTestInterval = 1
    
    def ShowPressures(self):
        dSets = []
        for t in self.Line.times:
            data = self.Line.GetByTimeStep(t)
            dSets.append(Gnuplot.Data(data.position[:,2], data.pressure))
        g = Gnuplot.Gnuplot()
        g.plot(*dSets)
        pdb.set_trace()

    def ShowVelocity(self):
        diameter = int(self.ScaleToLatticeUnits(self.Diameter))
        dSets = []
        for t in self.Plane.times:
            data = self.Plane.GetByTimeStep(t)
            cline = data[np.where(np.logical_and(data.grid[:, 1] == diameter/2, data.grid[:, 2] == diameter))]
            if len(dSets) == 0:
                x = cline.position[:, 0]
                y = cline.position[:, 1]
                pred = self.AnalyticSolutionFunction(x * pq.metre, y * pq.metre)
                pred_mps = self.ScaleForOutput(pred)
                inds = x.argsort()
                dSets.append(Gnuplot.Data(x[inds], pred_mps[inds], **{'with': 'lines lt -1'}))
                
            dSets.append(Gnuplot.Data(cline.position[:, 0], cline.velocity[:, 2]))
        g = Gnuplot.Gnuplot()
        g.plot(*dSets)
        pdb.set_trace()

    @property
    def AnalyticSolutionFunction(self):
        """Return a function that will calculate the predicted velocity
        field at a given point coordinate.
        """
        def f(xyz):
            return self.GenericSolution(xyz, self.Axis, self.Radius, self.ChosenPressureGradient, self.DynamicViscosity)
        return f

    @staticmethod
    def GenericSolution(xyz, axis, radius, gradP, eta):
        z = np.dot(xyz, axis)
        rSq = np.sum(xyz**2, axis=-1)
        axisDist2 = rSq - z**2
        return (radius**2 - axisDist2)[:,np.newaxis] * axis[np.newaxis, :] * gradP / (4. * eta)
    
    
    @cache_property
    def L2Error(self):
        data = self.ConvergedVelocityPlane
        
        xyz = data.position * pq.metre
        u = data.velocity * (pq.metre / pq.second)
        
        U = self.AnalyticSolutionFunction(xyz)

        errSq = np.mean(np.sum((u - U)**2,axis=-1))

        return np.sqrt(errSq) / self.MaximumVelocity

    @memo_property
    def PressureVsZ(self):
        lineData = self.ConvergedPressureLine
        z = np.dot(lineData.position, self.Axis)
        return z, lineData.pressure
    
    @memo_property
    def MeasuredPressureGradient(self):
        z, p = self.PressureVsZ
        nPts = len(z)
        pt1 = nPts / 4
        pt2 = 3 * pt1
        return -((pressure[pt2] - pressure[pt1]) / (z[pt2] - z[pt1])) * (self.mmHg / pq.metre)
    
    def SumSquareResidual(self, gradPressure, x, u):
        n = self.ScaleForOutput(self.Axis)
        r = self.ScaleForOutput(self.Radius)
        eta = self.ScaleForOutput(self.DynamicViscosity)
        U = self.GenericSolution(x, n, r, gradPressure, eta)
        err = U - u
        return np.sum(err**2)

    def ResidualDerivative(self, gradPressure, x, u):
        n = self.ScaleForOutput(self.Axis)
        r = self.ScaleForOutput(self.Radius)
        eta = self.ScaleForOutput(self.DynamicViscosity)
        U = self.GenericSolution(x, n, r, gradPressure, eta)
        return np.sum(2* U * (U - u) / gradPressure)

    @cache_property
    def MinErrPressureGradient(self):
        from scipy import optimize
        data = self.GetVelocityPlaneByTimeStep(self.Plane.times[-1])
        
        ans = optimize.fmin_cg(self.SumSquareResidual, self.ScaleForOutput(self.PressureGradient), args=(data.position, data.velocity), fprime=self.ResidualDerivative)
        return ans[0] * pq.pascal / pq.metre

    @property
    def ChosenPressureGradient(self):
        return getattr(self, self.Mode)
    Mode = 'PressureGradient'
    
    pass

if __name__ == "__main__":
    import argparse
    import glob
    
    p = argparse.ArgumentParser()
    # p.add_argument('-o', dest='outputfile', default=None,
    #                help='file to output the graph to, none implies plot to screen')
    p.add_argument('-j', dest='showjunk', action='store_true',
                   help='flag to include Junk-Yang results')
    # p.add_argument('-1', dest='first_const', default=0.4, type=float,
    #                help='coefficient for first order eye guide')
    # p.add_argument('-2', dest='second_const', default=0.7, type=float,
    #                help='coefficient for second order eye guide')
    p.add_argument('fab_results_dir', help='path to fabric results directory')
    
    # p.add_argument('results', nargs='+', help='results directories to process')
    args = p.parse_args()
    
    results = []
    for run in glob.glob(os.path.join(args.fab_results_dir, 'Cylinder_*')):
        try:
            res = CylinderResults.LoadFromSummary(run)
        except IOError:
            print "WARNING: '{}' Missing summary file".format(run)
            pass
        
        # Skip Junk
        if not args.showjunk and res.BoundaryCondition == 'jy':
            continue
            
        try:
            times = res.Plane.times
            if len(times) > 0:
                if len(times) < ((res.NumberOfTimeSteps + res.WarmUpSteps) / res.ExtractionInterval):
                    print "WARNING: '{}' Simulation incomplete - likely crashed or unfinished".format(run)
                    
                results.append(res)
            else:
                print "WARNING: '{}' No data in extracted property file - likely crashed or unfinished".format(run)
                pass
        except IOError:
            print "WARNING: '{}' No extracted property file - likely still queued".format(run)
            pass
        # print res.PressureGradient, " ", res.MeasuredPressureGradient, " ", res.MinErrPressureGradient

    resultsByRe = {}

    for res in results:
        try:
            resultsForRe = resultsByRe[int(res.ReynoldsNumber)]
        except KeyError:
            resultsForRe = resultsByRe[int(res.ReynoldsNumber)] = {}
            pass
        
        try:
            resForBC = resultsForRe[res.BoundaryCondition]
        except KeyError:
            resForBC = resultsForRe[res.BoundaryCondition] = {}
            pass
        
        try:
            resForVsCo = resForBC[res.VsCo]
        except KeyError:
            resForVsCo = resForBC[res.VsCo] = []
            pass
        
        resForVsCo.append(res)
        # Waste to sort multiple times, but who cares when it's only got 5 elements?
        resForVsCo.sort(lambda a, b: -cmp(a.VoxelSize, b.VoxelSize))        

    # Plot error vs resolution for all data
    def plotstyle(bc, vsCo):
        lc = {
            'sbb': 0,
            'bfl': 1,
            'gzs': 3,
            'jy': 5
            }[bc]
        pt = {
            'd3q15+lbgk': 1,
            'd3q15+mrt': 2,
            'd3q19+lbgk': 8,
            'd3q19+mrt': 10,
            'd3q27+lbgk': 4
            }[vsCo]
        # pt = {
        #     1.: 1,
        #     30.: 4,
        #     100.: 10,
        #     300.: 6
        #     }[re]
        return 'points pt {pt:d} lc {lc:d}'.format(**locals())

    def first_const(re):
        return {
            1.: 0.3,
            30.: 0.5,
            100.: 0.6,
            300: 0.7
            }[re]
    def second_const(re):
         return {
            1.: 0.7,
            30.: 10,
            100.: 5,
            300: 3.5
            }[re]

    g = Gnuplot.Gnuplot()
    base = 'convergence'
    g(r'''set logscale
set xrange [10:200]
set yrange [1e-5:2e-1]

set xlabel "$D/\\Delta x$"
set ylabel "$\\epsilon^2_u$"

set format y "$10^{{%L}}$"
set format x "$10^{{%L}}$"

set pointsize 2
# set key Left bmargin center reverse maxrows 5
set key off

set term epslatex standalone color size 17.8cm,7cm font 8
set output "{base}.tex"

set multiplot
'''.format(base=base))
    
    nPlots = 4
    rmargin = 0.8/28.7
    lmargin = 2.2/28.7
    plotwidth = (1. - rmargin - lmargin) / nPlots
    
    sizes = nPlots * [plotwidth]
    sizes[0] += lmargin
    sizes[-1]+= rmargin
    lefts = map(lambda i: sum(sizes[:i]), xrange(nPlots))
    rights = map(lambda i: sum(sizes[:i+1]), xrange(nPlots))
    
    reynolds = sorted(resultsByRe.keys())
    for i, re in enumerate(reynolds):
        resForRe = resultsByRe[re]
        
        plots = []
        for bc in ['sbb', 'bfl', 'gzs']:
            resForBC = resForRe[bc]
            keys = sorted(resForBC.keys())
            for vsCo in ['d3q15+lbgk', 'd3q19+lbgk', 'd3q27+lbgk']:
                resForVsCo = resForBC[vsCo]
                diameterVoxels = []
                error = []
                
                for res in resForVsCo:
                    diameterVoxels.append(res.ScaleToLatticeUnits(res.Diameter))
                    error.append(res.ScaleForOutput(res.L2Error))
                    continue
                dset = Gnuplot.Data(diameterVoxels, error, **{'with': plotstyle(bc, vsCo),
                                                              'title': res.VelocitySet.upper() + ', ' + bc.upper()})
                plots.append(dset)
        
        plots.append(Gnuplot.Func('%f*x**-1' % first_const(re),
                                  **{'title': None,
                                     'with': 'lines lc -1 lt 0'}))
        plots.append(Gnuplot.Func('%f*x**-2' % second_const(re),
                                  **{'title': None,
                                     'with': 'lines lc -1 lt 1'}))
        
        g(r'set title "(%s) $\\mathrm{Re}=%d$"' % (chr(ord('a') + i),re))
        g(r'set bmargin at screen 0.3')
        if i == 0:
            g('set lmargin')
            keywidth = 13.5/28.7
            g('set key reverse Left vertical maxrows 3 center top at screen 0.5,0.15')
            # g('set rmargin')
        else:
            g('set ylabel')
            g('set key off')
            g('set format y ""')
            g('set lmargin at screen %f' % lefts[i])
            
        if i == nPlots-1:
            g('set rmargin')
        else:
            g('set rmargin at screen %f' % rights[i])
        
        g.plot(*plots)
    
    g('unset multiplot')
    g('set output')
    os.system("epstopdf {base}-inc.eps && pdflatex {base} && open {base}.pdf".format(base=base))
        
        
    
