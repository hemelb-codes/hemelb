import numpy as np
import quantities as pq
import os.path
import Gnuplot
import warnings
from xml.etree import ElementTree

from hemeTools.parsers.extraction import ExtractedProperty

import _fixpath
from generation.WomersleyFlowCylinderGeneration import WomersleyFlowCylinderGeneration
from .memo import memo_property
from .cache import cache_property
from .HemeLbRunResults import HemeLbRunResults

import pdb
        
class WomersleyResults(HemeLbRunResults, WomersleyFlowCylinderGeneration):
    RunTypeName = 'Womersley'
    ConvergenceTestInterval = 4
    
    @property
    def PoiseuilleScaleFactor(self):
        from scipy.special import jv
        a = float(self.WomersleyNumber)
        
        ja = jv(0, 1j**1.5 * a)
        return 4 * abs(1 - 1./ja) / a**2

    @property
    def FittingFunction(self):
        def f(xyzt, wo, mag, phi):
            xyz = xyzt[:3,:].transpose()
            t = xyzt[3]
            n = self.Axis.view(np.ndarray)
            R = float(self.Radius)
            omega = 2 * np.pi / float(self.Period)
            
            return self.FullSolutionFunction(xyz, t, wo, mag, phi, n, R, omega)
        return f
    
    @staticmethod
    def FullSolutionFunction(xyz, t, wo, mag, phi, n, R, omega):
        """IMPORTANT: works only with non-quantities!
        """
        from scipy.special import jv
        
        arg = 1j**1.5 * wo
        bessel_denominator = jv(0, arg)

        arg *= WomersleyResults.GetRadius(n, xyz) / R
        bessel_numerator = jv(0, arg)

        posn_part = (1.0 - bessel_numerator / bessel_denominator)
        phase = np.exp(1j * (omega * t - phi))

        return -np.real(mag * posn_part * phase)
    
    @staticmethod
    def GetAxialComponent(n, x):
        return np.dot(x, n)
    
    @staticmethod
    def GetRadius(n, xyz):
        z = WomersleyResults.GetAxialComponent(n, xyz)
        return (np.sum(xyz**2, axis=-1) - z**2)**0.5
    
    @cache_property
    def FittedParameters(self):
        """IMPORTANT - returns non-quantities, in SI base units.
        """
        times = self.Plane.times
        iT1 = times.searchsorted(self.ConvergedTimeStep)
        for i, t in enumerate(times[iT1:iT1 + self.ConvergenceTestInterval]):
            tPhys = self.ScaleForOutput(t * self.TimeStep)
            data = self.GetVelocityPlaneByTimeStep(t)
            nPts = len(data)
            if i==0:
                x = np.zeros((4, self.ConvergenceTestInterval, nPts))
                y = np.zeros((self.ConvergenceTestInterval, nPts))
                
                r = self.GetRadius(self.Axis, data.position).view(np.ndarray)
                order = r.argsort()
                rSort = r[order]
                pass
            
            x[:3, i, :] = data.position.transpose()
            x[3, i, :] = tPhys
            y[i, :] = self.GetAxialComponent(self.Axis.view(np.ndarray), data.velocity)
            continue
            
        x.shape = (4, self.ConvergenceTestInterval * nPts)
        y.shape = (self.ConvergenceTestInterval * nPts)
        
        from scipy.optimize import curve_fit
        # Initial parameters
        p0 = (
            float(self.WomersleyNumber),
            float(self.PressureGradient / ( (2 * np.pi / self.Period) * self.Density)),
            0.
            )
        popt, pcov = curve_fit(self.FittingFunction, x, y, p0=p0)

        x.shape = (4, self.ConvergenceTestInterval, nPts)
        y.shape = (self.ConvergenceTestInterval, nPts)
        if self.ShouldDisplay:
            g = Gnuplot.Gnuplot()
            g.title(os.path.basename(self.RunDirectory))
            plots = []
            for i in xrange(self.ConvergenceTestInterval):
                r = self.GetRadius(self.Axis, x[:3, i, :].transpose())
                plots.append(Gnuplot.Data(r, y[i]))
                plots.append(Gnuplot.Data(r, self.FittingFunction(x[:, i,:], *popt)))
            g.plot(*plots)
            pdb.set_trace()
        return popt

    @property
    def FittedWomersleyNumber(self):
        return self.FittedParameters[0]
    
    @property
    def AnalyticSolutionFunction(self):
        # Uses quantities.
        def f(xyz, t):
            wo, mag, phi = self.FittedParameters
            n = self.Axis.view(np.ndarray)
            R = self.ScaleForOutput(self.Radius)
            omega = 2 * np.pi / self.ScaleForOutput(self.Period)
            oneD = self.FullSolutionFunction(self.ScaleForOutput(xyz), self.ScaleForOutput(t), wo, mag, phi, n, R, omega)
            return oneD[:, np.newaxis] * self.Axis[np.newaxis, :] * (pq.metre /  pq.second)
        
        return f
    
    ShouldDisplay = False
    @memo_property
    def L2Error(self):
        if self.ShouldDisplay:
            g = Gnuplot.Gnuplot()
            g.title(os.path.basename(self.RunDirectory))
            z = np.dot(self.Line.GetByIndex(0).position*pq.metre, self.Axis)
            order = z.argsort()
            z = z[order]

            times = self.Line.times[4:]
            phase = (times * self.TimeStep / self.Period) % (2. * np.pi)
            timeOrder = phase.argsort()
            times = times[timeOrder]

            plots = []
            for t in self.Line.times:
                line = self.Line.GetByTimeStep(t)
                p = line.pressure*self.mmHg
                phase = (t * self.TimeStep) / self.Period
                
                plots.append(Gnuplot.Data(self.ScaleToLatticeUnits(z), self.ScaleToLatticeUnits(p[order]), title='%d'%t))
            plots.append(
                Gnuplot.Func('{gradP} * x + {const}'.format(gradP=self.ScaleToLatticeUnits(self.PressureGradient),
                             const=self.ScaleToLatticeUnits(-0.5*self.PressureGradient*self.Length)),
                    **{'with': 'lines lc -1 lw 3'})
                )
            plots.append(
                Gnuplot.Func('{gradP} * x + {const}'.format(gradP=self.ScaleToLatticeUnits(-self.PressureGradient),
                                                            const=self.ScaleToLatticeUnits(0.5*self.PressureGradient*self.Length)),
            **{'with': 'lines lc -1 lw 3'})
                    )
            g.plot(*plots)
            pdb.set_trace()
                

        times = self.Plane.times
        iT1 = times.searchsorted(self.ConvergedTimeStep)
        errSq = 0 * (pq.metre / pq.second)**2

        if self.ShouldDisplay:
            plots = []
            
        for t in times[iT1:iT1 + self.ConvergenceTestInterval]:
            tPhys = t * self.TimeStep
            data = self.GetVelocityPlaneByTimeStep(t)
            xyz = data.position * pq.metre
            if self.ShouldDisplay:
                z = np.dot(xyz, self.Axis)
                r = (np.sum(xyz**2,axis=-1) - z**2)**0.5
                order = r.argsort()
                rSort = r[order]
            
            u = data.velocity * (pq.metre / pq.second)
            
            U = self.AnalyticSolutionFunction(xyz, tPhys)
            if self.ShouldDisplay:
                plots.append(Gnuplot.Data(rSort, self.GetAxialComponent(self.Axis, u)[order]))
                plots.append(Gnuplot.Data(rSort, self.GetAxialComponent(self.Axis, U)[order]))

            errSq += np.mean(np.sum((u - U)**2, axis=-1))

        if self.ShouldDisplay:
            g.plot(*plots)
            pdb.set_trace()
            
        errSq /= self.ConvergenceTestInterval
        return np.sqrt(errSq) / (self.PoiseuilleScaleFactor * self.MaximumVelocity)

    @memo_property
    def ParsedReport(self):
        return ElementTree.parse(os.path.join(self.RunDirectory, 'results', 'report.xml')).getroot()

    @property
    def ReportSites(self):
        return int(self.ParsedReport.find("geometry/sites").text)

    @property
    def ReportThreads(self):
        return int(self.ParsedReport.find("nodes/threads").text)

    @property
    def ReportSteps(self):
        return int(self.ParsedReport.find("results/steps/total").text)

    @property
    def ReportLbTime(self):
        for timer in self.ParsedReport.findall("timings/timer"):
            if timer.find("name").text == "LB calc only":
                return float(timer.find("mean").text)

        raise ValueError("Can't find the 'LB calc only' timer")

    @property
    def SUPS(self):
        return (self.ReportSites * self.ReportSteps / self.ReportLbTime)

    @property
    def SupsPerCore(self):
        # -1 for the rank 0 master core that does no LB
        try:
            return self.SUPS / (self.ReportThreads - 1)
        except IOError:
            print "Missing report " + self.RunDirectory
            self.MissingReports.add(self.RunDirectory)
            return -np.inf
        
    @property
    def SingleThreadSingleStepTime(self):
        try:
            return (self.ReportLbTime * (self.ReportThreads - 1)) / self.ReportSteps
        except IOError:
            print "Missing report " + self.RunDirectory
            self.MissingReports.add(self.RunDirectory)
            return np.inf
        
    MissingReports = set()
    pass

if __name__ == "__main__":
    import argparse
    import glob
    
    p = argparse.ArgumentParser()
    p.add_argument('-j', dest='showjunk', action='store_true',
                   help='flag to include Junk-Yang results')
    p.add_argument('--table', action='store_true')
    p.add_argument('--plot', action='store_true')
    p.add_argument('--perf', action='store_true')
    p.add_argument('--prof', action='store_true')

    p.add_argument('fab_results_dir', help='path to fabric results directory')
    
    args = p.parse_args()
    
    results = []
    for run in glob.glob(os.path.join(args.fab_results_dir, 'Womersley_Re*_D*')):
        try:
            res = WomersleyResults.LoadFromSummary(run)
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

    results = list(set(results))
    resultsByBc = {}

    for res in results:
        try:
            resForBc = resultsByBc[res.BoundaryCondition]
        except KeyError:
            resForBc = resultsByBc[res.BoundaryCondition] = {}
            pass
        
        try:
            resForVsCo = resForBc[res.VsCo]
        except KeyError:
            resForVsCo = resForBc[res.VsCo] = []
            pass
        
        resForVsCo.append(res)
        # Waste to sort multiple times, but who cares when it's only got 5 elements?
        resForVsCo.sort(lambda a, b: cmp(a.ReynoldsNumber, b.ReynoldsNumber))

    # Plot error vs resolution for all data
    def plotstyle(bc, vsCo):
        lc = {
            'sbb': 0,
            'bfl': 1,
            'gzs': 3,
            'jy': 5
            }[bc]
        pt = {
            'd3q15+lbgk': 4,
            'd3q15+mrt': 1,
            'd3q19+lbgk': 8,
            'd3q19+mrt': 2,
            'd3q27+lbgk': 12
            }[vsCo]
        return 'points pt {pt:d} lc {lc:d}'.format(**locals())
    
    if args.plot:
        plots = []
        for bc in ['sbb', 'bfl', 'gzs']:
            resForBC  = resultsByBc[bc]
            keys = sorted(resForBC.keys())
            for vsCo in keys:
                resForVsCo = resForBC[vsCo]
                woNo = []
                error = []
                for res in resForVsCo:
                    woNo.append(res.ScaleForOutput(res.WomersleyNumber))
                    error.append(res.ScaleForOutput(res.L2Error))
                    continue
                dset = Gnuplot.Data(woNo, error, **{'with': plotstyle(bc, vsCo),
                                                    'title': ' '.join((res.VelocitySet[-3:],
                                                                        res.CollisionOperator,
                                                        res.BoundaryCondition)).upper()})
                plots.append(dset)
                
        g = Gnuplot.Gnuplot()
        
        g(r'''
set xlabel "$\\alpha$"
set ylabel "$\\epsilon^2_u$"

set pointsize 2
set bmargin screen 0.399
set key Left reverse maxrows 5 at screen 0.46, 0.12 center width -2
''')
        
        base = "womersley"
        g('set term epslatex standalone color size 8.6cm,6cm font 8')
        g('set output "%s.tex"' % base)
        
        g.plot(*plots)
        g('set output')
        os.system("epstopdf {base}-inc.eps && pdflatex {base} && open {base}.pdf".format(base=base))
        print "Now insert the missing D3's into the .tex and adjust the legend text left and align the GZS  key entries."
        
    if args.table:
        lines = []
        for bc in ['sbb', 'bfl', 'gzs']:
            resForBC  = resultsByBc[bc]
            keys = sorted(resForBC.keys())
            for vsCo in keys:
                resForVsCo = resForBC[vsCo]
                woNo = []
                error = []
                
                line = [bc.upper(), vsCo[6:].upper(), vsCo[:5].upper()]
                
                for res in resForVsCo:
                    line.append('%.4f' % (res.ScaleForOutput(res.L2Error)))
                    continue
                lines.append(line)
        
        for line in lines:
            print ' & '.join(line) + r' \\'

    if args.perf:
        class Missing(object):
            SupsPerCore = -1.
            SingleThreadSingleStepTime = -1.
            pass
        nSites = 347401
        n_W, n_B = (42340, 305061)
        cols = {}
        for res in results:
            d = cols.get(res.VsCo, {})
            resList = d.get(res.BoundaryCondition, [])
            resList.append(res)
            d[res.BoundaryCondition] = resList
            cols[res.VsCo] = d
            
        cols['d3q15+mrt']['gzs'] = [Missing()]
        cols['d3q19+mrt']['gzs'] = [Missing()]
        
        header = ['']
        lines = {}
        for vsco in ('d3q15+lbgk', 'd3q19+lbgk', 'd3q27+lbgk', 'd3q15+mrt', 'd3q19+mrt'):
            header.append(vsco)
            
            sbbList = cols[vsco]['sbb']
            # T = time to do a single step on a single thread (for all sites)
            T_SBB = min(sbb.SingleThreadSingleStepTime for sbb in sbbList)
            # t = time to do a single site update
            # t_B = t_SBB (B = BULK)
            t_B = T_SBB / (n_W + n_B)
            
            for bc in ['sbb', 'bfl', 'gzs']:
                resList = cols[vsco][bc]
                line = lines.get(bc, [bc])
                
                maxSups = max(res.SupsPerCore for res in resList)
                line.append('%.1f' % (maxSups /1e6))
                
                T_BC = min(res.SingleThreadSingleStepTime for res in resList)
                t_BC = (T_BC - t_B * n_B) / n_W
                line.append('(%.1f)' % (t_BC/t_B))
                lines[bc] = line
                
        alllines = [header]
        for bc in ['sbb', 'bfl', 'gzs']:
            alllines.append(lines[bc])
            
        jyData = ((49.2, 84.2, 554),
                  (42.3, 45.3, 42),
                  (42.3, 44.0, 26),
                  (42.2, 43.7, 22))
        
        T_JY = min((x[1]-x[0]) / x[2] for x in jyData) * 63
        
        sups_JY = nSites / T_JY
        t_JY = (T_JY - t_B * n_B) / n_W
        alllines.append(['JY', '%.1f' % (sups_JY /1e6), '(%.1f)' % (t_JY/t_B)])
        
        print ' \\\\\n'.join(' & '.join(line) for line in alllines)

        print "The following reports are missing:"
        for missing in sorted(WomersleyResults.MissingReports):
            print missing

        
    if args.prof:
        plots = []
        for bc in ['bfl']:
            for vsCo in ['d3q15+lbgk']:
                resForLb = resultsByBc[bc][vsCo]
                
                res = resForLb[-1]
                
                wo, mag, phi = res.FittedParameters
                n = res.ScaleForOutput(res.Axis)
                R = res.ScaleForOutput(res.Radius)
                omega = 2 * np.pi / res.ScaleForOutput(res.Period)
                umax = res.ScaleForOutput(res.MaximumVelocity)
                
                times = res.Plane.times
                iT1 = times.searchsorted(res.ConvergedTimeStep)
                
                xyz = res.Plane.GetByIndex(0).position
                r = res.GetRadius(n, xyz) / R
                
                # nonzero = np.where(r != 0)
                # r = r[nonzero]
                ordering = r.argsort()
                r.sort()
                styles = [
                    (1, 0),
                    (4, 1,),
                    (6, 3),
                    (8, 4)
                    ]
                for t in times[iT1:iT1 + res.ConvergenceTestInterval]:
                    tSec = res.ScaleForOutput(t * res.TimeStep)
                    frac = res.ScaleForOutput(((t * res.TimeStep) % res.Period) / res.Period)
                    
                    uvw = res.Plane.GetByTimeStep(t).velocity
                    u = res.GetAxialComponent(n, uvw) / umax
                    
                    # u = u[nonzero]
                    u = u[ordering]
                    
                    plots.append(Gnuplot.Data(r, u, **{'title': '$t = %.2f T$'%frac,
                                                       'with': 'points pt %d lc %d' % styles.pop()}))
                    
                    U = res.FullSolutionFunction(xyz, tSec, wo, mag, phi, n, R, omega) / umax
                    # U = U[nonzero]
                    U = U[ordering]
                    
                    plots.append(Gnuplot.Data(r, U, **{'title': None, 'with': 'lines lt -1'}))
                    
        g = Gnuplot.Gnuplot()
if True:
        g(r'''
set xlabel "$r / R$"
set ylabel "$u / U_\\mathrm{max}$"

set bmargin screen 0.25
set xtics 0, 0.5
set ytics -0.05, 0.025
set yrange [-0.05:0.05]

set key reverse Right width -5 maxrows 1 at screen 0.95, 0.08
''')
        
        base = "womersley-profile"
        g('set term epslatex standalone color size 8.6cm,6cm font 8')
        g('set output "%s.tex"' % base)
        
        g.plot(*plots)
        g('set output')
        os.system("epstopdf {base}-inc.eps && pdflatex {base} && open {base}.pdf".format(base=base))

