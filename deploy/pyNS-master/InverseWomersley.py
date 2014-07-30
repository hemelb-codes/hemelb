#!/usr/bin/env python

## Program:   PyNS
## Module:    InverseWomersley.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

import os, sys, inspect, subprocess
from math import pi, cos, sin
from numpy.core.numeric import arange, zeros
from numpy.core.fromnumeric import mean
from numpy.lib.scimath import sqrt
from numpy.ma.core import exp
from numpy.lib.function_base import linspace


def add_bessel(name):
    cmd_folder = os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0])
    cmd_folder = cmd_folder+'/includes/%s' % name
    if cmd_folder not in sys.path:
            sys.path.insert(0, cmd_folder)

def intepreterArchitecture():
    is64 = False
    if sys.maxsize > 2**32:
        is64 = True
    return is64

try:
    from scipy.special import jn
except:
    if sys.version_info[0:2] == (2, 6):    
        if sys.platform == 'darwin' and intepreterArchitecture() is False:
            add_bessel('py26_macOsX_i386')
        if sys.platform == 'darwin' and intepreterArchitecture() is True:  #not available at the moment
            add_bessel('py26_macOsX_x8664')
        if sys.platform == 'linux2' and intepreterArchitecture() is False:
            add_bessel('py26_linux_i386')
        if sys.platform == 'linux2' and intepreterArchitecture() is True: #not available at the moment
            add_bessel('py26_linux_x8664')
        if sys.platform == 'win32':
            sys.exit("Windows users have to install scipy package (http://www.scipy.org/Download)")
           
    if sys.version_info[0:2] == (2, 7):    
        if sys.platform == 'darwin' and intepreterArchitecture() is False:
            add_bessel('py27_macOsX_i386')
        if sys.platform == 'darwin' and intepreterArchitecture() is True:
            add_bessel('py27_macOsX_x8664')
        if sys.platform == 'linux2' and intepreterArchitecture() is False:
            add_bessel('py27_linux_i386')
        if sys.platform == 'linux2' and intepreterArchitecture() is True:
            add_bessel('py27_linux_x8664')
        if sys.platform == 'win32':
            sys.exit("Windows users have to install scipy package (http://www.scipy.org/Download)")
                 
    try:
        from Bessel import Bessel as jn
    except ImportError:
        sys.exit("Bessel function is not supported. Please install scipy (http://www.scipy.org/Download) or install cython and launch python setup.py build_ext --inplace. \n More info at http://docs.cython.org/")

class InverseWomersley(object):
    '''
    This class computes Womersley Wall Shear Stress from pressure signal,
    according to Inverse Womersley Method explained in Cezeaux et al. 1997.
    This class provides the following methods:
    SetSimulationContext : a method for setting simulation context.
    SetNetworkGraph: a method for setting NetworkGraph input.
    SetNetworkMesh : a method for setting NetworkMesh.
    SetNetworkSolution: a method for setting NetworkSolution.
    SetFlowSignal: a method for setting Flow Signal for specific mesh.
    GetVelFromQ : a method for computing velocity from flow volume using the Inverse Womersley Method.
    GetTaoFromQ: a method for computing wall shear stress from flow volume using the Inverse Womersley Method.
    GetWssPeaks: a method for computing wall shear stress peak values (on radius array) from flow volume using the Inverse Womersley Method.
    SaveVelocityProfile: a method for showing velocity profile animation over the fractional radius using WX library (you have to switch the import statement in line 19)
    ShowVelocityProfile: a method for saving velocity profile movie over the fractional radius using MenCoder library.
    PlotWss: a method for plotting Wss signal.
    PlotFlow: a method for plotting Flow signal.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.fourierModes = []
        self.signalF = []
        self.signal = []
        self.t = 0.0
        self.tPeriod = None
        self.nSteps = None
        self.dt = None
        self.viscosity = None
        self.density = 1.05e3
        self.samples = 500.0
        self.dtPlot = None
        self.nHarmonics = 10
        self.radius = None
        self.q = 0.0
        self.tau = 0.0
        self.yTao = 0.0
        self.y = 0.0
    
    def SetSimulationContext(self, context):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = context
        try:
            self.viscosity = self.SimulationContext.Context['dynamic_viscosity']
        except KeyError:
            print "Error, Please set Dynamic Viscosity[Pa*s] in Boundary Conditions XML File"
            raise
        try:
            self.density = self.SimulationContext.Context['blood_density']
        except KeyError:
            print "Error, Please set Blood Density[kg*m^3] in Boundary Conditions XML File" 
            raise  
        try:
            self.tPeriod = self.SimulationContext.Context['period']
        except KeyError:
            print "Error, Please set period in Boundary Conditions XML File"
            raise            
        
    def SetNetworkMesh(self, networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetNetworkGraph(self, networkGraph):
        '''
        Setting NetworkMesh
        '''
        self.NetworkGraph = networkGraph
    
    def SetNetworkSolutions(self, networkSolutions):
        '''
        Setting Solutions
        '''
        self.NetworkSolutions = networkSolutions
        
    
    def SetFlowSignal(self, el, flowsig):
        '''
        Setting Flow Signal for specific mesh
        '''
        for sig in flowsig:
            self.signal.append(float(sig)) #flow in m3/s
        self.nSteps = arange(0,len(self.signal),1)
        self.dt = self.tPeriod/(len(self.nSteps)-1)
        self.dtPlot = self.tPeriod/self.samples
    
    def GetVelFromQ(self,el):
        '''
        Computing velocity profile in terms of the flow rate,
        using inverse womersley method of Cezeaux et al.1997
        '''
        self.radius = mean(el.Radius)        
        self.Res = el.R
        self.length = el.Length
        self.Name = el.Name
        Flow = mean(self.signal)
        
        #WOMERSLEY NUMBER
        self.alpha = self.radius * sqrt((2.0 *pi*self.density)/(self.tPeriod*self.viscosity))
        self.Wom = self.alpha
        self.Re = (2.0*Flow*self.SimulationContext.Context['blood_density'])/(pi*self.radius*self.SimulationContext.Context['dynamic_viscosity'])
        
        #FOURIER SIGNAL
        k = len(self.signal)
        n = 0
        while n < (self.nHarmonics):
            An = 0
            Bn = 0
            for i in arange(k):
                An += self.signal[i] * cos(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
                Bn += self.signal[i] * sin(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
            An = An * (2.0/k)
            Bn = Bn * (2.0/k)
            self.fourierModes.append(complex(An, Bn))
            n+=1
            
        self.fourierModes[0] *= 0.5   #mean Flow, as expected. It's defined into xml input file.
        
        self.Steps = linspace(0,self.tPeriod,self.samples)
        
        self.VelRadius = {}
        self.VelRadiusSteps = {}
        self.VelocityPlot = {}
        for step in self.Steps:
            self.Velocity = {}
            y = -1 # raggio da -1 a 1, 200 punti.
            while y <=1.:
                self.VelRadius[y] = 2*(1.0**2 - y**2)*self.fourierModes[0]
                y+=0.01
                
            k=1
            while k < self.nHarmonics: 
                cI = complex(0.,1.)
                cA = (self.alpha * pow((1.0*k),0.5)) * pow(cI,1.5)      
                c1 = 2.0 * jn(1, cA)
                c0 = cA * jn(0, cA)
                cT = complex(0, -2.0*pi*k*self.t/self.tPeriod)  
                y=-1 #da -1 a 1  #y=0 #centerline
                while y<=1.0:
                    '''vel computation'''
                    c0_y = cA * jn(0, (cA*y))
                    vNum = c0-c0_y
                    vDen = c0-c1
                    vFract = vNum/vDen
                    cV = self.fourierModes[k] * exp(cT) * vFract
                    self.VelRadius[y] += cV.real #valore di velocity riferito al raggio adimensionalizzato
                    self.Velocity[y] = self.VelRadius[y].real
                    y+=0.01
                k+=1
            
            unsortedRadii = []
            for rad, vel in self.Velocity.iteritems():
                unsortedRadii.append(rad)
            radii = sorted(unsortedRadii)
            
            self.VelPlot = []
            for x in radii:
                for rad, vel in self.Velocity.iteritems():
                    if x == rad:
                        self.VelPlot.append(vel*(100.0/(self.radius**2*pi)))                   
            self.VelocityPlot[step] = self.VelPlot
            self.t += self.dtPlot
            
        
    def GetTaoFromQ(self,el):
        '''
        Computing wall shear stress in terms of the flow rate,
        using inverse womersley method of Cezeaux et al.1997
        '''
        self.radius = mean(el.Radius)        
        self.Res = el.R
        self.length = el.Length
        self.Name = el.Name
        
        #WOMERSLEY NUMBER
        self.alpha = self.radius * sqrt((2.0 *pi*self.density)/(self.tPeriod*self.viscosity))
        
        #FOURIER SIGNAL
        k = len(self.signal)
        n = 0
        while n < (self.nHarmonics):
            An = 0
            Bn = 0
            for i in arange(k):
                An += self.signal[i] * cos(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
                Bn += self.signal[i] * sin(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
            An = An * (2.0/k)
            Bn = Bn * (2.0/k)
            self.fourierModes.append(complex(An, Bn))
            n+=1
        
        self.Steps = linspace(0,self.tPeriod,self.samples)
        self.WssSignal = []  
        self.Tauplot = []
       
        for step in self.Steps:
            self.tao = -self.fourierModes[0].real * 2.0 
            
            k=1
            while k < self.nHarmonics:  
                cI = complex(0.,1.)
                cA = (self.alpha * pow((1.0*k),0.5)) * pow(cI,1.5)  
                c1 = 2.0 * jn(1, cA)
                c0 = cA * jn(0, cA)
                cT = complex(0, -2.0*pi*k*self.t/self.tPeriod)  
                '''tao computation'''
                taoNum = self.alpha**2*cI**3*jn(1,cA)
                taoDen = c0-c1
                taoFract = taoNum/taoDen
                cTao = self.fourierModes[k] * exp(cT) * taoFract
                self.tao += cTao.real
                k+=1

            self.tao *= -(self.viscosity/(self.radius**3*pi))
            self.Tauplot.append(self.tao*10) #dynes/cm2
            self.WssSignal.append(self.tao)
            self.t += self.dtPlot
            
        return self.WssSignal #Pascal
     
    def GetWssPeaks(self,el, flowsig):
        '''
        This method returns Wss peak along the element.
        Wss in s=0 and s=1 is computed.
        '''
        r0 = el.Radius[0]
        r1 = el.Radius[len(el.Radius)-1]
        r01Signal = []
        
        for sig in flowsig:
            r01Signal.append(sig)
        
        self.nSteps = arange(0,len(r01Signal),1)
        self.dt = self.tPeriod/(len(self.nSteps)-1)
        self.dtPlot = self.tPeriod/self.samples
        fourierModes = []
        
        #Computing for s=0
        r0WssSignal = []
        #WOMERSLEY NUMBER
        r0Alpha = r0 * sqrt((2.0 *pi*self.density)/(self.tPeriod*self.viscosity))
        
        #Computing for s=1
        r1WssSignal = []
        #WOMERSLEY NUMBER
        r1Alpha = r1 * sqrt((2.0 *pi*self.density)/(self.tPeriod*self.viscosity))
        
        k01 = len(r01Signal)
        n = 0
        while n < (self.nHarmonics):
            An = 0
            Bn = 0
            for i in arange(k01):
                An += r01Signal[i] * cos(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
                Bn += r01Signal[i] * sin(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
            An = An * (2.0/k01)
            Bn = Bn * (2.0/k01)
            fourierModes.append(complex(An, Bn))
            n+=1
            
        self.Steps = linspace(0,self.tPeriod,self.samples)
        for step in self.Steps:
            tao0 = -fourierModes[0].real * 2.0
            tao1 = -fourierModes[0].real * 2.0
            k=1
            while k < self.nHarmonics: 
                cI = complex(0.,1.)
                cA_0 = (r0Alpha * pow((1.0*k),0.5)) * pow(cI,1.5)
                c1_0 = 2.0 * jn(1, cA_0)  
                c0_0 = cA_0 * jn(0, cA_0)
                cA_1 = (r1Alpha * pow((1.0*k),0.5)) * pow(cI,1.5)
                c1_1 = 2.0 * jn(1, cA_1)  
                c0_1 = cA_1 * jn(0, cA_1) 
                cT = complex(0, -2.0*pi*k*self.t/self.tPeriod)
                '''R0: Wall shear stress computation'''
                taoNum_0 = r0Alpha**2*cI**3*jn(1,cA_0)
                taoDen_0 = c0_0-c1_0
                taoFract_0 = taoNum_0/taoDen_0
                cTao_0 = fourierModes[k] * exp(cT) * taoFract_0
                tao0 += cTao_0.real
                '''R1: Wall shear stress computation'''
                taoNum_1 = r1Alpha**2*cI**3*jn(1,cA_1)
                taoDen_1 = c0_1-c1_1
                taoFract_1 = taoNum_1/taoDen_1
                cTao_1 = fourierModes[k] * exp(cT) * taoFract_1
                tao1 += cTao_1.real
                   
                k+=1
                
            tao0 *= -(self.viscosity/(r0**3*pi))    
            r0WssSignal.append(tao0)
            tao1 *= -(self.viscosity/(r1**3*pi))    
            r1WssSignal.append(tao1)
            self.t += self.dtPlot
            
        r0Peak = max(r0WssSignal)
        r1Peak = max(r1WssSignal)
        return r0Peak,r1Peak
    
    def SaveVelocityProfile(self, meshid, daystr):
        '''
        This method plots velocity profile into png files and makes 
        an avi file from png set. Mencoder is required.
        '''
        try:
            import matplotlib
            matplotlib.use('Agg') 
            from matplotlib.pyplot import plot, xlabel, ylabel, title, savefig, ylim, axis, clf
        except:
            sys.exit("VelocityProfile methods require matplotlib package (http://matplotlib.sourceforge.net.\n")
        
        #Create temporary image and videos directories'''
        if not os.path.exists ('tmp/'):
            os.mkdir('tmp/')
        if not os.path.exists ('videos/'):
            os.mkdir('videos/')
        if not os.path.exists ('videos/%s' % daystr):
            os.mkdir('videos/%s' % daystr)
        
        not_found_msg = """
        The mencoder command was not found;
        mencoder is used by this script to make an avi file from a set of pngs.
        It is typically not installed by default on linux distros because of
        legal restrictions, but it is widely available.
        """
        try:
            subprocess.check_call(['mencoder'])
        except subprocess.CalledProcessError:
            print "mencoder command was found"
            pass # mencoder is found, but returns non-zero exit as expected
            # This is a quick and dirty check; it leaves some spurious output
            # for the user to puzzle over.
        except OSError:
            print not_found_msg
            sys.exit("quitting\n")
        
        self.count = 0
        orderingStep = []
        for step, vel in self.VelocityPlot.iteritems():
            orderingStep.append(step)
            lenStep = len(self.VelocityPlot)
            lenVel = len(vel)
        orderedStep = sorted(orderingStep)
        orderedVel = []
        for timeStep in orderedStep:
            for step, vel in self.VelocityPlot.iteritems():
                if timeStep == step:
                    orderedVel.append(vel)
        maxY = 0
        minY = 0
        for vY in orderedVel:
            if max(vY) > maxY:
                maxY = max(vY)
            if min(vY) < minY:
                minY = min(vY)
        x = linspace(-1.0,1.0,lenVel)   # Values to be plotted on the x-axis.
        ylim(ymax=maxY)
        ylim(ymin=minY)
        i = 0
        print 'Computing Images for velocity profile...'
        while i<lenStep:
            plot(x,orderedVel[i],'r-',linewidth = 3)
            axis((x[1],x[-1],minY,maxY))
            xlabel('Fractional radius')
            ylabel('Velocity profile ($cm/s$)')
            title (str(self.Name)+' radius($mm$)= '+str(round(self.radius*1e3,1))+' Reynolds N.= '+str(round(self.Re,0))+' Womersley N.= '+str(round(self.Wom,2)))
            filename = str('%04d' % i) + '.png'
            savefig('tmp/'+filename, dpi=100)
            clf()
            i+=1
        print 'Making movie from images..'    
        command = ('mencoder',
           'mf://tmp/*.png',
           '-mf',
           'type=png:w=800:h=600:fps=25',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'videos/%s%s.avi' % (daystr,self.Name))
        print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
        subprocess.check_call(command)
        print "\n\n The movie was written"

    def ShowVelocityProfile(self, meshid):
        '''
        This method plots an animated representation of the velocity profile
        evolving in time using wx python library.
        '''
        try:
            import matplotlib
            matplotlib.use('WXAgg') 
            from matplotlib.pyplot import xlabel, ylabel, title, close, figure, ylim, show
        except:
            sys.exit("VelocityProfile methods require matplotlib package (http://matplotlib.sourceforge.net.\n")
        try:
            from wx import GetApp,EVT_CLOSE, EVT_IDLE
        except:
            sys.exit("ShowVelocityProfile method requires wxpython package (http://www.wxpython.org).\n")
        self.count = 0
        orderingStep = []
        for step, vel in self.VelocityPlot.iteritems():
            orderingStep.append(step)
            lenVel = len(vel)
            if step == 0.0:
                firstVel = vel
        orderedStep = sorted(orderingStep)
        orderedVel = []
        for timeStep in orderedStep:
            for step, vel in self.VelocityPlot.iteritems():
                if timeStep == step:
                    orderedVel.append(vel)
        orderedVel.reverse()
        maxY = 0
        minY = 0
        for vY in orderedVel:
            if max(vY) > maxY:
                maxY = max(vY)
            if min(vY) < minY:
                minY = min(vY)
        fig = figure()
        ax = fig.add_subplot(111)
        t = linspace(-1.0,1.0,lenVel)
        line, = ax.plot(t, firstVel,'r-',linewidth = 3)
        ylim(ymax=maxY)
        ylim(ymin=minY)
        xlabel('Fractional radius')
        ylabel('Velocity profile ($cm/s$)')
        title ('Mean radius($mm$)= '+str(round(self.radius*1e3,0))+' Reynolds N.= '+str(round(self.Re,0))+' Womersley N.= '+str(round(self.Wom,2))) 
        
        '''WX ANIMATION'''
        def update_line(idleevent):
            if orderedVel == []:
                for timeStep in orderedStep:
                    for step, vel in self.VelocityPlot.iteritems():
                        if timeStep == step:
                            orderedVel.append(vel)
                orderedVel.reverse()
                self.count +=1
            if self.count >=2:
                self.count = 0
                EVT_CLOSE(GetApp(), update_line)
                close()
            line.set_ydata(orderedVel.pop())
            fig.canvas.draw_idle() 
        EVT_IDLE(GetApp(), update_line)
        show()
        
    def PlotWss(self, meshid, imagpath):
        '''
        This method plots Wss signal and returns peak wss.
        '''
        try:
            import matplotlib
            matplotlib.use('Agg') #switch to matplotlib.use('WXAgg') if you want to show and not save velocity profile.
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close, ylim
        except:
            sys.exit("PlotWss method requires matplotlib package (http://matplotlib.sourceforge.net.\n")
        
        tplot = linspace(0, self.tPeriod, len(self.Tauplot))
        plot(tplot, self.Tauplot,'g-',linewidth = 3, label = 'WSS')
        minY = 0
        for w in self.Tauplot:
            if w < minY:
                minY = w
        
        if minY != 0:
            plot(tplot, zeros(len(self.Tauplot)),':',linewidth = 1)
            
        ylim(ymin=minY)
        
        xlabel('Time ($s$)')
        ylabel('Wall shear stress ($dyne/cm^2$)')
        title ('Wss'+' peak:'+str(round(max(self.Tauplot),1))+' mean:'+str(round(mean(self.Tauplot),1))+' min:'+str(round(min(self.Tauplot),1)))    
        legend()
        savefig(imagpath+str(meshid)+'_'+str(self.Name)+'_wss.png')
        print "Wss, MeshId", meshid, self.Name, "=", str(round(max(self.Tauplot),1)), "$dyne/cm^2$"
        close()
        return (round(max(self.Tauplot),1))