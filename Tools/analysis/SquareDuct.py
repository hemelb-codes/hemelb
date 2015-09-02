import numpy as np
import quantities as pq
import os.path
from hemeTools.parsers.extraction import ExtractedProperty
import Gnuplot
from generation.SquareDuctGeneration import SquareDuctGeneration
from .memo import memo_property
from .cache import cache_property
from .HemeLbRunResults import HemeLbRunResults

class SquareDuctResults(SquareDuctGeneration, HemeLbRunResults):
    @classmethod
    def LoadFromSummary(cls, run):
        smy = os.path.join(run, 'config.smy')
        inst = super(SquareDuctResults, cls).LoadFromSummary(smy)
        inst.RunDirectory = run
        return inst
    
    @memo_property
    def Plane(self):
        planeFileName = os.path.join(self.RunDirectory, 'results', 'Extracted', 'plane.xtr')
        return ExtractedProperty(planeFileName)

        return ans

    @memo_property
    def Line(self):
        planeFileName = os.path.join(self.RunDirectory, 'results', 'Extracted', 'line.xtr')
        return ExtractedProperty(planeFileName)

    def ShowPressures(self):
        dSets = []
        for t in self.Line.times:
            data = self.Line.GetByTimeStep(t)
            dSets.append(Gnuplot.Data(data.position[:,2], data.pressure))
        g = Gnuplot.Gnuplot()
        g.plot(*dSets)
        pdb.set_trace()

    def ShowVelocity(self):
        pdb.set_trace()
        sideLength = int(self.ScaleToLatticeUnits(self.SideLength))
        dSets = []
        for t in self.Plane.times:
            data = self.Plane.GetByTimeStep(t)
            cline = data[np.where(np.logical_and(data.grid[:, 1] == sideLength/2, data.grid[:, 2] == sideLength))]
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
        field at a given (y,z) coordinate.
        """
        L = self.SideLength / 2.
        eta = self.DynamicViscosity
        gradPressure = self.ChosenPressureGradient
        
        Power = lambda x, n: x**n
        Pi = np.pi
        Cos = np.cos
        Cosh = np.cosh
        Sech = lambda x: 1. / np.cosh(x)

        return lambda y, z : (gradPressure*Power(L,2)*(1 - Power(z,2)/Power(L,2) + 
            4*((-8*Cos((Pi*z)/(2.*L))*Cosh((Pi*y)/(2.*L))*Sech(Pi/2.))/Power(Pi,3) + 
               (8*Cos((3*Pi*z)/(2.*L))*Cosh((3*Pi*y)/(2.*L))*Sech((3*Pi)/2.))/
                (27.*Power(Pi,3)) - (8*Cos((5*Pi*z)/(2.*L))*Cosh((5*Pi*y)/(2.*L))*
                  Sech((5*Pi)/2.))/(125.*Power(Pi,3)) + 
               (8*Cos((7*Pi*z)/(2.*L))*Cosh((7*Pi*y)/(2.*L))*Sech((7*Pi)/2.))/
                (343.*Power(Pi,3)) - (8*Cos((9*Pi*z)/(2.*L))*Cosh((9*Pi*y)/(2.*L))*
                  Sech((9*Pi)/2.))/(729.*Power(Pi,3)) + 
               (8*Cos((11*Pi*z)/(2.*L))*Cosh((11*Pi*y)/(2.*L))*Sech((11*Pi)/2.))/
                (1331.*Power(Pi,3)) - (8*Cos((13*Pi*z)/(2.*L))*Cosh((13*Pi*y)/(2.*L))*
                  Sech((13*Pi)/2.))/(2197.*Power(Pi,3)) + 
               (8*Cos((15*Pi*z)/(2.*L))*Cosh((15*Pi*y)/(2.*L))*Sech((15*Pi)/2.))/
                (3375.*Power(Pi,3)) - (8*Cos((17*Pi*z)/(2.*L))*Cosh((17*Pi*y)/(2.*L))*
                  Sech((17*Pi)/2.))/(4913.*Power(Pi,3)) + 
               (8*Cos((19*Pi*z)/(2.*L))*Cosh((19*Pi*y)/(2.*L))*Sech((19*Pi)/2.))/
                (6859.*Power(Pi,3)))))/(2.*eta)

    @memo_property
    def ConvergedVelocityPlane(self):
        t = self.Plane.times[-1]
        return self.GetVelocityPlaneByTimeStep(t)
    
    def GetVelocityPlaneByTimeStep(self, t):
        sideLength = int(self.ScaleToLatticeUnits(self.SideLength))
        data = self.Plane.GetByTimeStep(t)
        singleplane = data[np.where(data.grid[:, 2] == sideLength)]
        ordering = np.argsort(sideLength * singleplane.grid[:, 0] + singleplane.grid[:, 1])
        singleplane = singleplane[ordering]
        singleplane.shape = (sideLength, sideLength)
        return singleplane
    def GetSecondVelocityPlaneByTimeStep(self, t):
        sideLength = int(self.ScaleToLatticeUnits(self.SideLength))
        data = self.Plane.GetByTimeStep(t)
        singleplane = data[np.where(data.grid[:, 2] == sideLength+1)]
        ordering = np.argsort(sideLength * singleplane.grid[:, 0] + singleplane.grid[:, 1])
        singleplane = singleplane[ordering]
        singleplane.shape = (sideLength, sideLength)
        return singleplane
        
    @property
    def L2Error(self):
        data = self.ConvergedVelocityPlane
        x = data.position[..., 0] * pq.metre
        y = data.position[..., 1] * pq.metre
        
        u = data.velocity * (pq.metre / pq.second)
        
        U = self.AnalyticSolutionFunction(x, y)
        
        errSq = np.mean(u[..., 0]**2 + u[..., 1]**2 + (u[..., 2] - U)**2)

        return np.sqrt(errSq) / self.MaximumVelocity

    @cache_property
    def VelocityConvergenceError(self):
        ultimate = self.GetVelocityPlaneByTimeStep(self.Plane.times[-1])
        penultimate = self.GetVelocityPlaneByTimeStep(self.Plane.times[-2])
        vDiffNormSq = np.sum((ultimate.velocity - penultimate.velocity)**2, axis=-1)
        return np.sqrt(np.mean(vDiffNormSq)) * (pq.metre / pq.second) / self.MaximumVelocity

    @cache_property
    def VelocityAxisError(self):
        ultimate = self.GetVelocityPlaneByTimeStep(self.Plane.times[-1])
        penultimate = self.GetSecondVelocityPlaneByTimeStep(self.Plane.times[-1])
        vDiffNormSq = np.sum((ultimate.velocity - penultimate.velocity)**2, axis=-1)
        return np.sqrt(np.mean(vDiffNormSq)) * (pq.metre / pq.second) / self.MaximumVelocity
    
    @cache_property
    def MeasuredPressureGradient(self):
        lineData = self.ConvergedPressureLine
        nPts = len(lineData)
        pt1 = lineData[nPts / 4]
        pt2 = lineData[3 * (nPts / 4)]
        return -((pt2.pressure - pt1.pressure) / (pt2.position[2] - pt1.position[2])) * (self.mmHg / pq.metre)
    
    @memo_property
    def ConvergedPressureLine(self):
        t = self.Line.times[-1]
        lineData = self.Line.GetByTimeStep(t)
        sideLength = int(self.ScaleToLatticeUnits(self.SideLength))
        select = np.where(np.logical_and(lineData.grid[:, 0] == sideLength / 2,
                                         lineData.grid[:, 1] == sideLength / 2))
        zs = lineData.grid[:, 2][select]
        ordering = (select[0][zs.argsort()],)
        return lineData[ordering]

    def FitFunc(self, gradPressure, xdata, ydata):
        L = self.ScaleForOutput(self.SideLength / 2.)
        eta = self.ScaleForOutput(self.DynamicViscosity)
        
        Power = lambda x, n: x**n
        Pi = np.pi
        Cos = np.cos
        Cosh = np.cosh
        Sech = lambda x: 1. / np.cosh(x)

        y = xdata[0, :]
        z = xdata[1, :]
        return np.sum((ydata - (gradPressure*Power(L,2)*(1 - Power(z,2)/Power(L,2) + 
            4*((-8*Cos((Pi*z)/(2.*L))*Cosh((Pi*y)/(2.*L))*Sech(Pi/2.))/Power(Pi,3) + 
               (8*Cos((3*Pi*z)/(2.*L))*Cosh((3*Pi*y)/(2.*L))*Sech((3*Pi)/2.))/
                (27.*Power(Pi,3)) - (8*Cos((5*Pi*z)/(2.*L))*Cosh((5*Pi*y)/(2.*L))*
                  Sech((5*Pi)/2.))/(125.*Power(Pi,3)) + 
               (8*Cos((7*Pi*z)/(2.*L))*Cosh((7*Pi*y)/(2.*L))*Sech((7*Pi)/2.))/
                (343.*Power(Pi,3)) - (8*Cos((9*Pi*z)/(2.*L))*Cosh((9*Pi*y)/(2.*L))*
                  Sech((9*Pi)/2.))/(729.*Power(Pi,3)) + 
               (8*Cos((11*Pi*z)/(2.*L))*Cosh((11*Pi*y)/(2.*L))*Sech((11*Pi)/2.))/
                (1331.*Power(Pi,3)) - (8*Cos((13*Pi*z)/(2.*L))*Cosh((13*Pi*y)/(2.*L))*
                  Sech((13*Pi)/2.))/(2197.*Power(Pi,3)) + 
               (8*Cos((15*Pi*z)/(2.*L))*Cosh((15*Pi*y)/(2.*L))*Sech((15*Pi)/2.))/
                (3375.*Power(Pi,3)) - (8*Cos((17*Pi*z)/(2.*L))*Cosh((17*Pi*y)/(2.*L))*
                  Sech((17*Pi)/2.))/(4913.*Power(Pi,3)) + 
               (8*Cos((19*Pi*z)/(2.*L))*Cosh((19*Pi*y)/(2.*L))*Sech((19*Pi)/2.))/
                (6859.*Power(Pi,3)))))/(2.*eta))**2)
    
    def FitGradientFunc(self, gradPressure, xdata, ydata):
        L = self.ScaleForOutput(self.SideLength / 2.)
        eta = self.ScaleForOutput(self.DynamicViscosity)
        
        Power = lambda x, n: x**n
        Pi = np.pi
        Cos = np.cos
        Cosh = np.cosh
        Sech = lambda x: 1. / np.cosh(x)

        y = xdata[0, :]
        z = xdata[1, :]

        profile = (Power(L,2)*(1 - Power(z,2)/Power(L,2) + 
            4*((-8*Cos((Pi*z)/(2.*L))*Cosh((Pi*y)/(2.*L))*Sech(Pi/2.))/Power(Pi,3) + 
               (8*Cos((3*Pi*z)/(2.*L))*Cosh((3*Pi*y)/(2.*L))*Sech((3*Pi)/2.))/
                (27.*Power(Pi,3)) - (8*Cos((5*Pi*z)/(2.*L))*Cosh((5*Pi*y)/(2.*L))*
                  Sech((5*Pi)/2.))/(125.*Power(Pi,3)) + 
               (8*Cos((7*Pi*z)/(2.*L))*Cosh((7*Pi*y)/(2.*L))*Sech((7*Pi)/2.))/
                (343.*Power(Pi,3)) - (8*Cos((9*Pi*z)/(2.*L))*Cosh((9*Pi*y)/(2.*L))*
                  Sech((9*Pi)/2.))/(729.*Power(Pi,3)) + 
               (8*Cos((11*Pi*z)/(2.*L))*Cosh((11*Pi*y)/(2.*L))*Sech((11*Pi)/2.))/
                (1331.*Power(Pi,3)) - (8*Cos((13*Pi*z)/(2.*L))*Cosh((13*Pi*y)/(2.*L))*
                  Sech((13*Pi)/2.))/(2197.*Power(Pi,3)) + 
               (8*Cos((15*Pi*z)/(2.*L))*Cosh((15*Pi*y)/(2.*L))*Sech((15*Pi)/2.))/
                (3375.*Power(Pi,3)) - (8*Cos((17*Pi*z)/(2.*L))*Cosh((17*Pi*y)/(2.*L))*
                  Sech((17*Pi)/2.))/(4913.*Power(Pi,3)) + 
               (8*Cos((19*Pi*z)/(2.*L))*Cosh((19*Pi*y)/(2.*L))*Sech((19*Pi)/2.))/
                (6859.*Power(Pi,3)))))/(2.*eta)
        
        return -2.0 * np.sum((ydata - gradPressure * profile) * profile)
    
    @cache_property
    def MinErrPressureGradient(self):
        from scipy import optimize
        data = self.GetVelocityPlaneByTimeStep(self.Plane.times[-1])
        x = data.position[..., 0:2]
        x.shape = (x.shape[0] * x.shape[1], x.shape[2])
        x = x.transpose()
        
        y = data.velocity[..., 2].flatten()
        
        ans = optimize.fmin_cg(self.FitFunc, self.ScaleForOutput(self.PressureGradient), args=(x, y), fprime=self.FitGradientFunc)
        return ans[0] * pq.pascal / pq.metre

    @property
    def ChosenPressureGradient(self):
        return getattr(self, self.Mode)
    Mode = 'MinErrPressureGradient'
    
    pass

if __name__ == "__main__":
    import pdb
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument('-o', dest='outputfile', default=None,
                   help='file to output the graph to, none implies plot to screen')
    p.add_argument('results', nargs='+', help='results directories to process')
    args = p.parse_args()
    
    results = [SquareDuctResults.LoadFromSummary(run) for run in args.results]
    results.sort(lambda a, b: cmp(a._SideLengthVoxels, b._SideLengthVoxels))

    sideLengthVoxels = []
    error = []
    conv = []
    axis = []
    for res in results:
        sideLengthVoxels.append(res.ScaleForOutput(res.SideLength / res.VoxelSize))
        error.append(res.ScaleForOutput(res.L2Error))
        conv.append(res.ScaleForOutput(res.VelocityConvergenceError))
        axis.append(res.ScaleForOutput(res.VelocityAxisError))
        continue
    
    # Sort them by Reynolds no
    resultsByRe = {}
    for res in results:
        try:
            resultsByRe[res.ScaleForOutput(res.ReynoldsNumber)].append(res)
        except KeyError:
            resultsByRe[res.ScaleForOutput(res.ReynoldsNumber)] = [res]
            pass
    
    g = Gnuplot.Gnuplot()
    g('set logscale')
    
    dSets = []
    for re, reResults in resultsByRe.iteritems():
        sideLengthVoxels = []
        error = []
        for res in reResults:
            sideLengthVoxels.append(res.ScaleForOutput(res.SideLength / res.VoxelSize))
            error.append(res.ScaleForOutput(res.L2Error))
            continue
        dSets.append(Gnuplot.Data(sideLengthVoxels, error, title=str(re)))
    dSets.append(Gnuplot.Func('10 * x**-2', title=None))
    if args.outputfile is not None:
        g('set term post')
        g('set output "%s"' % args.outputfile)
        
    g.plot(*dSets)
    
        
        
