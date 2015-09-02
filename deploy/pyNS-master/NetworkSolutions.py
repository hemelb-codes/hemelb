#!/usr/bin/env python

## Program:   PyNS
## Module:    NetworkSolutions.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

import csv
from DofMap import DofMap
from InverseWomersley import InverseWomersley
from numpy.core.fromnumeric import mean
from numpy.core.numeric import array, zeros
from math import pi
from numpy.lib.function_base import linspace
from numpy.core.numeric import arange
from numpy.lib.scimath import sqrt
from numpy.ma.core import ceil
from json import dump 
from xml.etree import ElementTree as etree
import sys

class NetworkSolutions(object):
    '''
    NetworkSolutions elaborates results for post processing.
    This class provides he following methods:
    SetNetworkGraph: a method for setting NetworkGraph input
    SetNetworkMesh: a method for setting NetworkMesh input
    SetSolutions: a method for setting Solutions input
    SetSimulationContext: a method for setting SimulationContext input.
    #Json  Methods:
    WriteJsonInfo: a method for creating the json which includes simulation information (elements and adaptation steps if adaptation is active).
    WriteJson: a method for creating a json for a specified mesh which includes all information needed for postprocessing.
    WriteJsonAdapt: a method for creating a json for a specified mesh which includes information about vascular adaptation.
    #General Methods:
    PlotBrachial: a method for plotting Brachial flows, pressure and wss (for each segment) in the same figure. (Requires matplotlib package)
    PlotRadial: a method for plotting Radial flows, pressure and wss (for each segment) in the same figure. (Requires matplotlib package)
    PlotCephalic: a method for plotting Cephalic flows, pressure and wss (for each segment) in the same figure. (Requires matplotlib package)
    WriteToXML: a method for writing solutions in XML Solutions File.
    GetSolution: a method for plotting flow, pressure and WSS for specific entity.
    #Flow Methods:
    GetInflow: a method for plotting input flow function.
    PlotFlow: a method for plotting mean flow for a single mesh. (Requires matplotlib package)
    PlotFlowComparative: a method for plotting brachial, radial and ulnar mean flow. (Requires matplotlib package)
    GetFlowSignal: a method for returning flow signal for specific mesh.
    WriteFlowOutput: a method for writing flow output values for a specific mesh in a .txt file.
    PlotReynolds: a method for plotting reynolds number for a single mesh. (Requires matplotlib package)
    WriteReynoldsOutput: a method for writing reynolds number output values for a specific mesh in a .txt file.
    #Pressure Methods:
    PlotPressure: a method for plotting mean pressure for a single mesh. (Requires matplotlib package)
    PlotPressureTwo: a method for plotting pressures for a couple of meshes. (Requires matplotlib package)
    PlotPressureComparative: a method for plotting brachial, radial and ulnar mean pressure. (Requires matplotlib package)
    GetPressureSignal: a method for returning pressure signal for specific mesh.
    PlotPressureDrop : a method for plotting pressure drop for specific mesh. (Requires matplotlib package)
    WritePressureInput: a method for writing pressure input values for a specific mesh in a .txt file.
    WritePressureOutput: a method for writing pressure output values for a specific mesh in a .txt file.
    #Wall Shear Stress Methods:
    PlotPWSS: a method for plotting mean WSS(Poiseuille) for a single mesh. (Requires matplotlib package)
    PlotPWSSComparative: a method for plotting brachial, radial and ulnar mean WSS (Poiseuille). (Requires matplotlib package)
    GetPWSSSignal: a method for returning WSS signal(Poiseuille) for specific mesh.
    PlotWSS: a method for plotting mean WSS for a single mesh. (Requires matplotlib package)
    GetWSSSignal: a method for returning WSS signal for specific mesh.
    WriteWSSOutput: a method for writing WSS output values for a specific mesh in a .txt file.
    #Other Methods:
    PulseWaveVelocity: a method for computing Pulse Wave Velocity(m/s).
    If SuperEdge2 is specified, PWV is computed between first and second superedge,
    otherwise PWV is computed over a single superedge.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkGraph = None
        self.NetworkMesh = None
        self.SimulationContext = None
        self.DofMap = None
        self.Solutions = None
        self.t = None
        self.TimeStep = None
        self.Period = None
        self.CardiacFreq = None
        self.Cycles = None
        self.images = None
        self.dayFlow = {} #{element.Id:Q}
        self.dayWssP = {} #{element.Id:taoPeak}
        self.dayPressure = {} #{element.Id:pressure}
        self.dayDiameter = {} #{element.Id:diameter}
        self.flowDirection = {'(+)':"from node1 to node2", '(-)':"from node2 to node1"} #Flow direction is assumed positive if flow goes from node1 to node2.
        
    def SetNetworkGraph(self,networkGraph):
        '''
        Setting NetworkMesh
        '''
        self.NetworkGraph = networkGraph
        
    def SetNetworkMesh(self,networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        self.DofMap = DofMap()
        self.DofMap.SetNetworkMesh(networkMesh)
        self.DofMap.Build()
    
    def SetSimulationContext(self, simulationContext):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = simulationContext
        try:
            self.TimeStep = simulationContext.Context['timestep']
        except KeyError:
            print "Error, Please set timestep in Simulation Context XML File"
            raise
        try:
            self.Period = self.SimulationContext.Context['period']
            self.CardiacFreq = int(self.Period/self.TimeStep)
        except KeyError:
            print "Error, Please set period in Simulation Context XML File"
            raise
        try:
            self.Cycles = self.SimulationContext.Context['cycles']
        except KeyError:
            print "Error, Please set cycles number in Simulation Context XML File"
            raise 
        self.t = linspace(self.TimeStep,self.Period,self.CardiacFreq)
    
    def SetSolutions(self, solutions):
        '''
        Setting Solutions
        '''
        self.Solutions = solutions
        
    def SetImagesPath(self, imagDict):
        '''
        Setting images directory
        '''
        for name, path in imagDict.iteritems():
            if name == 'im':
                self.images = path
                self.f_images = path
                self.p_images = path
                self.w_images = path
        for name, path in imagDict.iteritems():
            if name == 'f':
                self.f_images = path
            if name == 'p':
                self.p_images = path
            if name == 'w':
                self.w_images = path
    
    #INPUT GNUID
    
    def GetGnuidInformation(self, PatientId, mesh_radius):
        '''This method provides information to be used as Gnuid input.'''
        ro = self.SimulationContext.Context['blood_density']/1000
        mu = self.SimulationContext.Context['dynamic_viscosity']*10
        period_pyns = self.SimulationContext.Context['period']
        frequency_pyns = 1./period_pyns
        
        #find mesh
        
        for element in self.NetworkMesh.Elements:
            if element.Type == 'Anastomosis':
                anastomosis = element
                prox_artery = anastomosis.Proximal
                
                dofs = prox_artery.GetPoiseuilleDofs()
                prox_artery_p1 = self.Solutions[(self.DofMap.DofMap[prox_artery.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] #pressure [Pa]
                prox_artery_p2 = self.Solutions[(self.DofMap.DofMap[prox_artery.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):] #pressure [Pa]
                prox_artery_q = mean(((prox_artery_p1-prox_artery_p2)/prox_artery.R))*6e7  #flow [mL/min] 
                prox_artery_radius = (max(prox_artery.Radius))*100  # radius [cm]
                      
                dist_artery = anastomosis.Distal
                dofs = dist_artery.GetPoiseuilleDofs()
                dist_artery_p1 = self.Solutions[(self.DofMap.DofMap[dist_artery.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] #pressure [Pa]
                dist_artery_p2 = self.Solutions[(self.DofMap.DofMap[dist_artery.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):] #pressure [Pa]
                dist_artery_q = mean(((dist_artery_p1-dist_artery_p2)/dist_artery.R))*6e7  #flow [mL/min] 
                
                
        prox_artery_mass_flow = (prox_artery_q * ro)/60  #mass flow [g/s]
        dist_artery_mass_flow = (dist_artery_q * ro)/60  #mass flow [g/s]
        womersley = prox_artery_radius*sqrt((2*pi*frequency_pyns*ro)/mu) #womersley number
        gnuid_frequency = ((womersley**2)*mu)/((float(mesh_radius)**2)*2*pi*ro) #frequency for gnuid
        gnuid_period = 1./gnuid_frequency  #period for gnuid
        gnuid_scaling_factor = dist_artery_mass_flow/prox_artery_mass_flow #gnuid scaling factor
        
        txtpath = 'Results/' + PatientId + '/gnuid_parameters.txt'
        text_file = open(txtpath, "w")
        text_file.write("GNUID Parameters\n")
        text_file.write("Proximal artery mass flow " + str(prox_artery_mass_flow) + "\n")
        text_file.write("Distal artery mass flow " +  str(dist_artery_mass_flow) + "\n")
        text_file.write("Womersley number " +  str(womersley) + "\n")
        text_file.write("Blood Density " +  str(ro) + "\n")
        text_file.write("Blood Viscosity " +  str(mu) + "\n")
        text_file.write("Gnuid frequency " +  str(gnuid_frequency) + "\n")
        text_file.write("Gnuid period " +  str(gnuid_period) + "\n")
        text_file.write("Gnuid scaling factor " + str(gnuid_scaling_factor) + "\n")
        text_file.close()      
                   
    # JSON METHODS
    
    def WriteJsonInfo(self, days, elements, PatientId):
        '''
        ''' 
        info = {}
        info['time'] = []
        info['elements'] = []
        
        if days > 0:
            info['adaptation']=True
        else:
            info['adaptation']=False
        for d in range(-1,days+1):
            if d >0:
                info['time'].append(d*10)
            else:
                info['time'].append(d)
        info['time'].sort()
        
        for el in elements:
            if el.Type == 'WavePropagation' or el.Type == 'Resistance':
                info['elements'].append(el.Name)
                info['elements'].sort()
        path = 'Results/' + PatientId + '/json/info.json'
        f = open(path,'w')
        dump(info, f)
        f.close()
    
    def WriteJson(self, meshid, time, excludeWss, PatientId):
        '''
        This method writes a json file for each mesh.
        '''
        meshInfo = {}
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()     
                p1 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                Flow = (p1-p2)/element.R
                Wss = []
                Reynolds = []
                elName = element.Name
                meshInfo['meshId']=str(meshid)
                meshInfo['name']=str(elName)
                meshInfo['mean_pressure']=str(round(mean(p1/133.32),2))+' mmHg'
                meshInfo['min_pressure'] = str(0)
                min_p = round(min(p1/133.32),2)
                if min_p < 0:
                    meshInfo['min_pressure'] = str(min_p) 
                meshInfo['mean_flow']=str(round(mean(Flow*6e7),2))+' mL/min'
                meshInfo['min_flow'] = str(0)
                min_q = round(min(Flow*6e7),2)
                if min_q < 0:
                    meshInfo['min_flow'] = str(min_q)
                    
                meshInfo['items'] = []
                timeValues = {}
                timeValues['flow'] = []
                timeValues['pressure'] = []
                
                self.dayFlow[element.Name] = (round(mean(Flow*6e7),1))
                self.dayPressure[element.Name] = (round(mean(p1/133.32),1))
                
                if element.Type != "Resistance":
                    Radius = element.Radius
                    Reynolds = (2.0*Flow*self.SimulationContext.Context['blood_density'])/(pi*max(Radius)*self.SimulationContext.Context['dynamic_viscosity'])
                    if excludeWss is False:
                        Wss = self.GetWSSSignal(element)
                        tWss = linspace(0, self.Period, len(Wss))
                    meshInfo['length']=str(round(element.Length*1e2,2))+' cm'
                    meshInfo['diameter_min']=str(round((min(element.Radius)*2e3),2))+' mm'
                    meshInfo['diameter_max']=str(round((max(element.Radius)*2e3),2))+' mm'
                    if excludeWss is False:    
                        meshInfo['mean_wss']=str(round(mean(Wss*10),2))+' dynes/cm<sup>2</sup>'
                        meshInfo['min_wss'] = str(0)
                        min_wss = round(min(Wss*10),2)
                        if min_wss < 0:
                            meshInfo['min_wss'] = str(min_wss)
                    meshInfo['mean_re']=str(round(mean(Reynolds),1))
                    meshInfo['min_re'] = str(0)
                    min_re = round(min(Reynolds),1)
                    if min_re < 0:
                        meshInfo['min_re'] = str(min_re)
                    timeValues['wss'] = []
                    timeValues['re'] = []
                    if excludeWss is False:
                        self.dayWssP[element.Name] = (round(max(Wss*10),2))
                    try:
                        self.dayDiameter[element.Name] = (round(element.dayRadius[time][0]*2e3,2))
                    except:
                        self.dayDiameter[element.Name] = (round(element.Radius[0]*2e3,2))
                
                
                             
        i=0
        for q in Flow:
            timeValues['flow'].append([self.t[i],(round(q*6e7,2))])
            i+=1
        i=0
        for p in p1:
            timeValues['pressure'].append([self.t[i],(round(p/133.32,2))])
            i+=1
        if excludeWss is False:
            i=0
            for w in Wss:
                timeValues['wss'].append([tWss[i],(round(w*10,2))])
                i+=1
        i=0
        for re in Reynolds:
            timeValues['re'].append([self.t[i],(round(re,2))])
            i+=1
        
        meshInfo['items'].append(timeValues)
        if time > 0: 
            path = 'Results/' + PatientId + '/json/'+str(time*10)+'_'+str(elName)+'.json'
        else:
            path = 'Results/' + PatientId + '/json/'+str(time)+'_'+str(elName)+'.json'
        f = open(path,'w')
        dump(meshInfo, f)
        f.close()
        
        
    def WriteJsonAdapt(self, adaptation, PatientId):
        '''
        If adaptation algorithm is active, this method 
        will write a json file for each mesh including
        information about vessel adaptation.
        '''
            
        meshInfo = {}
        
        for element in self.NetworkMesh.Elements:
            if element.Type == 'WavePropagation':
                elName = element.Name
                meshInfo['meshId']=str(element.Id)
                meshInfo['name']=str(elName)
                meshInfo['items'] = []
                timeValues = {}
                timeValues['flow'] = []
                timeValues['pressure'] = []
                timeValues['wssP'] = []
                timeValues['diameter'] = []
                
                for day,sol in adaptation.solutions.iteritems():
                    if day == -1:
                        try:
                            timeValues['flow'].append([day,sol.dayFlow[element.Name]])
                            timeValues['pressure'].append([day,sol.dayPressure[element.Name]])
                            timeValues['wssP'].append([day,sol.dayWssP[element.Name]])
                            timeValues['diameter'].append([day,sol.dayDiameter[element.Name]])
                        except KeyError:
                            pass
                    if day != -1:
                        timeValues['flow'].append([day*10,sol.dayFlow[element.Name]])
                        timeValues['pressure'].append([day*10,sol.dayPressure[element.Name]])
                        timeValues['wssP'].append([day*10,sol.dayWssP[element.Name]])
                        timeValues['diameter'].append([day*10,sol.dayDiameter[element.Name]])
                   
                timeValues['flow'].sort()
                min_q = 0
                for q in timeValues['flow']:
                    if q[1] < min_q:
                        min_q = q[1]
                meshInfo['min_q'] = str(min_q)
                timeValues['pressure'].sort()
                timeValues['wssP'].sort()
                timeValues['diameter'].sort()
                meshInfo['items'].append(timeValues)
                path = 'Results/' + PatientId + '/json/adapt_'+str(elName)+'.json'
                f = open(path,'w')
                dump(meshInfo, f)
                f.close()
        
    # GENERAL METHODS
    
    def PlotBrachial(self):
        '''
        This method plots Brachial flows, pressure
        and wss (for each segment) in the same figure.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find('rachial') != -1:
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    Flow = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])   
                    WSSPoiseuille = ((4.0*element.eta)/pi) * (Flow/mean(element.Radius)**3)
                    print "Name: ", element.Name, " PressureIN(mmHg): ", mean(PressureIN)/133.3223684211, " MeanFlow(mL/min): ", mean(Flow)*6.0e7,  " MeanWSS(Pa): ", mean(WSSPoiseuille)
                    plot(self.t, Flow*6e7, colourvector[indexcolour],linewidth = 3, label = 'Flow '+element.Name)
                    xlabel('Time ($s$)')
                    ylabel('Flow ($mL/min$)')
                    title ('Brachial Flow')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'brachial_flow.png')
                close()
                indexcolour = 0 
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211  
                    plot(self.t, PressureIN, colourvector[indexcolour],linewidth = 3, label = 'Pressure '+element.Name)
                    xlabel('Time ($s$)')
                    ylabel('Pressure ($mmHG$)')
                    title ('Brachial Pressure')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'brachial_pressure.png')
                close()
                indexcolour = 0                        
                
    def PlotRadial(self):
        '''
        This method plots Radial flows, pressure
        and wss (for each segment) in the same figure.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
            
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find('radial') != -1:
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    Flow = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])   
                    WSSPoiseuille = ((4.0*element.eta)/pi) * (Flow/mean(element.Radius)**3)
                    print "Name: ", element.Name, " PressureIN(mmHg): ", mean(PressureIN)/133.3223684211, " MeanFlow(mL/min): ", mean(Flow)*6.0e7,  " MeanWSS(Pa): ", mean(WSSPoiseuille)
                    plot(self.t, Flow*6e7, colourvector[indexcolour],linewidth = 3, label = 'Flow '+element.Name)
                    xlabel('Time ($s$)')
                    ylabel('Flow ($mL/min$)')
                    title ('Radial Flow')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'radial_flow.png')
                close()
                indexcolour = 0 
                
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211  
                    plot(self.t, PressureIN, colourvector[indexcolour],linewidth = 3, label = 'Pressure '+element.Name)
                    xlabel('Time ($s$)')
                    ylabel('Pressure ($mmHG$)')
                    title ('Radial Pressure')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'radial_pressure.png')
                close()
                indexcolour = 0   
        
    def PlotCephalic(self):
        '''
        This method plots Cephalic flows, pressure
        and wss (for each segment) in the same figure.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find('cephalic') != -1:
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    Flow = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])   
                    WSSPoiseuille = ((4.0*element.eta)/pi) * (Flow/mean(element.Radius)**3)
                    print "Name: ", element.Name, " PressureIN(mmHg): ", mean(PressureIN)/133.3223684211, " MeanFlow(mL/min): ", mean(Flow)*6.0e7,  " MeanWSS(Pa): ", mean(WSSPoiseuille)
                    plot(self.t, Flow*6e7, colourvector[indexcolour],linewidth = 3, label = element.Name)
                    xlabel('Time ($s$)')
                    ylabel('Flow ($mL/min$)')
                    title ('Cephalic Flow')    
                    legend(loc=0)
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'cephalic_flow.png')
                close()
                indexcolour = 0 
                
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211  
                    plot(self.t, PressureIN, colourvector[indexcolour],linewidth = 3, label = 'Pressure '+element.Name)
                    xlabel('Time ($s$)')
                    ylabel('Pressure ($mmHG$)')
                    title ('Cephalic Pressure')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'cephalic_pressure.png')
                close()
                indexcolour = 0   
    
    # FLOW METHODS
    
    def GetInflow(self, flow):
        '''
        This method plots inflow function
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        t = linspace(0.0,self.Period+self.TimeStep,self.CardiacFreq).reshape((1,ceil(self.Period/self.TimeStep+1.0)))        
        plot(t[0,:], flow[0,:]*6e7, 'r-',linewidth = 3, label = 'Inlet Flow')   #red line
        xlabel('Time ($s$)')
        ylabel('Flow ($mL/min$)')
        title ('InFlow: '+str(mean(flow)*6.0e7)+' $mL/min$')    
        legend()
        savefig(self.images+'inflow.png')
        close()
    
    def PlotDaysBrachial(self, wdir, adaptation):
        '''
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        Flow = {-1:self.SimulationContext.Context['brachial_flow']}
        Diameter = {}
        t = [-1]
        x = 0
        yFlows=[]
        yDiams=[]
        
        
        for el in self.NetworkMesh.Elements:
            if el.Name ==  'brachial_prox_1':
                for day,sol in adaptation.solutions.iteritems():
                    if day != -1:
                        Flow.update({day:sol.dayFlow[el.Id]})
                        t.append(x)
                        x+=1
                        try:
                            Diameter.update({day:el.dayRadius[day][0]*2e3})
                        except:
                            Diameter.update({day:el.Radius[0]*2e3})
                            
        Diameter[-1]=Diameter[0]           
        tPlot = []
        for day in sorted(t):
            yFlows.append(Flow[day])
            yDiams.append(Diameter[day])
            tPlot.append(day*10)
            
       
        plot(tPlot,yDiams,marker='o', linestyle='-',linewidth = 2)
        minY = 0
        ylim(ymin=minY)
        xlabel('Time ($days$)')
        ylabel('Diameter ($mm$)')
        brachDiam = wdir+'/Brachial_days_diam.png'
        savefig(wdir+'/Brachial_days_diam.png')
        close()
        
        plot(tPlot,yFlows,marker='o', linestyle='-',linewidth = 2)
        ylim(ymin=minY)
        xlabel('Time ($days$)')
        ylabel('Flow ($mL/min$)')
        brachFlow = wdir+'/Brachial_days_flow.png'
        savefig(wdir+'/Brachial_days_flow.png')
        close()
        meanBrachFlow = Flow[4]
        
        return brachDiam, brachFlow, meanBrachFlow
    
    
    def GetDistalPressure(self, adaptation):
        '''
        This method returns distal pressure in the index finger (mmHg)
        '''
        for element in self.NetworkMesh.Elements:
            if element.Name == 'index_finger':
                for day,sol in adaptation.solutions.iteritems():
                    if day == 4:
                        meshid = element.Id
                        dofs = element.GetPoiseuilleDofs()
                        distalPressure = mean(self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])
               
        return round(distalPressure/133.32,2)
        
    def PlotDistalPressure(self, wdir, distalPressures):
        '''
        '''
        try:
            from matplotlib.pyplot import plot, savefig, close, ylim, figure
        except ImportError:
            MatPlotLibError()
        pressure = []
        ind = arange(4)
        pos = arange(4)
        width = 0.5
        for va, p in distalPressures.iteritems():
            for x in ind:
                if va == x:
                    pressure.append(p)
                    
        fig = figure()
        ax = fig.add_subplot(111)
        rects = ax.bar(0.2+ind, pressure, width, color='c')
        ax.set_title('Distal pressure prediction in different AVFs')
        ax.set_ylabel('Index finger pressure at 40 days ($mmHg$)')
        ax.set_xticks(0.45+pos)
        ax.set_xticklabels( ('RCEE', 'RCES', 'BCES', 'BBES') )
        plot([0,4], [60,60],'--', color='k', linewidth=2 )
        maxY = max(pressure)
        maxY+=0.02*maxY
        ylim(ymax=maxY)

        pressureImage = wdir+'/distalPressure.png'
        savefig(wdir+'/distalPressure.png')
        close()
        
        return pressureImage
    
    
    def Compare(self, wdir, compare):
        '''
        '''
        try:
            from matplotlib.pyplot import plot, savefig, close, ylim, figure
        except ImportError:
            MatPlotLibError()
        flow = []
        std = []
        ind = arange(4)
        pos = arange(4)
        width = 0.5
        for va, q in compare.iteritems():
            for x in ind:
                if va == x:
                    flow.append(q)
                    std.append(q*0.2)
        
        fig = figure()
        ax = fig.add_subplot(111)
        rects = ax.bar(0.2+ind, flow, width, color='c', yerr = std, ecolor='r')
        ax.set_title('Comparing flow volume in different AVFs')
        ax.set_ylabel('Brachial artery flow volume at 40 days ($mL/min$)')
        ax.set_xticks(0.45+pos)
        ax.set_xticklabels( ('RCEE', 'RCES', 'BCES', 'BBES') )
        plot([0,2], [400,400],'--', color='k', linewidth=2 )
        maxY = max(flow)+max(std)
        maxY+=0.02*maxY
        ylim(ymax=maxY)

        compareImage = wdir+'/compare.png'
        savefig(wdir+'/compare.png')
        close()
        
        return compareImage
        
    
    def Plots(self, wdir, day, name):
        '''
        '''
        try: 
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        timeImages = []     
        for element in self.NetworkMesh.Elements:
            if element.Name == name:
                meshid = element.Id
                el = element
                dofs = element.GetPoiseuilleDofs()
                        
                p1 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                
                if mean(p1)>=mean(p2):
                    flowDirection = self.flowDirection['(+)']
                    Flow = (p1-p2)/element.R
                else:
                    flowDirection = self.flowDirection['(-)']
                    Flow = (p2-p1)/element.R
               
        self.dayFlow[meshid] = (round(mean(Flow*6e7),1))
        
        plot(self.t, Flow*6e7, 'r-',linewidth = 3, label = 'Volumetric flow rate')   #red line
        minY = 0
        for q in Flow*6e7:
            if q < minY:
                minY = q       
        if minY != 0:
            plot(self.t, zeros(len(Flow)),':',linewidth = 1)
        ylim(ymin=minY)
        xlabel('Time ($s$)')
        ylabel('Volumetric flow rate ($mL/min$)')
        if flowDirection == self.flowDirection['(+)']:
            title ('Flow (+)'+' peak:'+str(round(max(Flow*6e7),1))+' mean:'+str(round(mean(Flow*6e7),1))+' min:'+str(round(min(Flow*6e7),1)))    
        if flowDirection == self.flowDirection['(-)']:
            title ('Flow (-)'+' peak:'+str(round(max(Flow*6e7),1))+' mean:'+str(round(mean(Flow*6e7),1))+' min:'+str(round(min(Flow*6e7),1)))    
        legend()
            
        if name.__contains__('brachial'):
            timeImages.append(wdir+'/%s_brachflow.png' % day)
            savefig(wdir+'/%s_brachflow.png' % day)
            close()
            
        if name.__contains__('radial'):
            timeImages.append(wdir+'/%s_radflow.png' % day)
            savefig(wdir+'/%s_radflow.png' % day)
            close()
            
        
        plot(self.t, p1/133.32, 'b-', linewidth = 3, label = 'Pressure Signal')   #blue line
        minY = 0
        for p in p1/133.32:
            if p < minY:
                minY = p
        if minY != 0:
            plot(self.t, zeros(len(p1)),':',linewidth = 1)
        ylim(ymin=minY)
        xlabel('Time ($s$)')
        ylabel('Pressure ($mmHg$)')
        title ('Pressure'+' peak:'+str(round(max(p1/133.32),1))+' mean:'+str(round(mean(p1/133.32),1))+' min:'+str(round(min(p1/133.32),1)))    
        legend()
        
        if name.__contains__('brachial'):
            timeImages.append(wdir+'/%s_brachpressure.png' % day)
            savefig(wdir+'/%s_brachpressure.png' % day)
            close()
            
        if name.__contains__('radial'):
            timeImages.append(wdir+'/%s_radpressure.png' % day)
            savefig(wdir+'/%s_radpressure.png' % day)
            close()
            
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        inverseWomersley.SetFlowSignal(self.GetFlowSignal(el))
        inverseWomersley.GetTaoFromQ(el)
        tplot = linspace(0, inverseWomersley.tPeriod, len(inverseWomersley.Tauplot))
        plot(tplot, inverseWomersley.Tauplot,'g-',linewidth = 3, label = 'WSS')
        minY = 0
        for w in inverseWomersley.Tauplot:
            if w < minY:
                minY = w
        if minY != 0:
            plot(tplot, zeros(len(inverseWomersley.Tauplot)),':',linewidth = 1)     
        ylim(ymin=minY)
        xlabel('Time ($s$)')
        ylabel('Wall shear stress ($dyne/cm^2$)')
        title ('Wss'+' peak:'+str(round(max(inverseWomersley.Tauplot),1))+' mean:'+str(round(mean(inverseWomersley.Tauplot),1))+' min:'+str(round(min(inverseWomersley.Tauplot),1)))    
        legend()
        
        if name.__contains__('brachial'):
            timeImages.append(wdir+'/%s_brachwss.png' % day)
            savefig(wdir+'/%s_brachwss.png' % day)
            close()
            
        if name.__contains__('radial'):
            timeImages.append(wdir+'/%s_radwss.png' % day)
            savefig(wdir+'/%s_radwss.png' % day)
            close()
                    
        return timeImages
    
    def PlotDaysRadial(self, wdir, adaptation):
        '''
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        Flow = {-1:self.SimulationContext.Context['radial_flow']}
        Diameter = {}
        t = [-1]
        x = 0
        yFlows=[]
        yDiams=[]
        
        for el in self.NetworkMesh.Elements:
            if el.Name ==  'radial_prox_3':
                for day,sol in adaptation.solutions.iteritems():
                    if day != -1:
                        Flow.update({day:sol.dayFlow[el.Id]})
                        t.append(x)
                        x+=1
                        try:
                            Diameter.update({day:el.dayRadius[day][0]*2e3})
                        except:
                            Diameter.update({day:el.Radius[0]*2e3})
                            
        Diameter[-1]=Diameter[0]       
        tPlot = []
        for day in sorted(t):
            yFlows.append(Flow[day])
            yDiams.append(Diameter[day])
            tPlot.append(day*10)
            
        minY = 0
        plot(tPlot,yDiams,marker='o', linestyle='-',linewidth = 2)
        ylim(ymin=minY)
        xlabel('Time ($days$)')
        ylabel('Diameter ($mm$)')
        radDiam = wdir+'/Radial_days_diam.png'
        savefig(wdir+'/Radial_days_diam.png')
        close()
        
        plot(tPlot,yFlows,marker='o', linestyle='-',linewidth = 2)
        ylim(ymin=minY)
        xlabel('Time ($days$)')
        ylabel('Flow ($mL/min$)')
        radFlow = wdir+'/Radial_days_flow.png'
        savefig(wdir+'/Radial_days_flow.png')
        close()
        
        return radDiam, radFlow
        
    def PlotFlow(self, meshid):
        '''
        This method plots mean flow for a single mesh.
        Flow volume is considered as absolute value.
        Direction of the flow is assumed + if goes from node1 to node2. Instead
        flowDirection is assumed to be -.
        '''
        try:  
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                name = self.NetworkGraph.Edges[self.NetworkMesh.MeshToGraph[meshid]].Name
                dofs = element.GetPoiseuilleDofs()
                        
                p1 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                
                if mean(p1)>=mean(p2):
                    flowDirection = self.flowDirection['(+)']
                    Flow = (p1-p2)/element.R
                else:
                    flowDirection = self.flowDirection['(-)']
                    Flow = (p2-p1)/element.R
               
                print "Flow, MeshId ", element.Id, ' ', element.Name, " = " , mean(Flow)*6e7, "mL/min ", flowDirection
            
        self.dayFlow[meshid] = (round(mean(Flow*6e7),1))
        
        plot(self.t, Flow*6e7, 'r-',linewidth = 3, label = 'Flow Output')   #red line
        minY = 0
        for q in Flow*6e7:
            if q < minY:
                minY = q
                   
        if minY != 0:
            plot(self.t, zeros(len(Flow)),':',linewidth = 1)
        ylim(ymin=minY)
        xlabel('Time ($s$)')
        ylabel('Volumetric flow rate ($mL/min$)')
        if flowDirection == self.flowDirection['(+)']:
            title ('Flow (+)'+' peak:'+str(round(max(Flow*6e7),1))+' mean:'+str(round(mean(Flow*6e7),1))+' min:'+str(round(min(Flow*6e7),1)))    
        if flowDirection == self.flowDirection['(-)']:
            title ('Flow (-)'+' peak:'+str(round(max(Flow*6e7),1))+' mean:'+str(round(mean(Flow*6e7),1))+' min:'+str(round(min(Flow*6e7),1)))    
        legend()
        savefig(self.f_images + meshid + '_' + name +'_flow.png')
        close()
    
    def PlotVelocity(self, meshid):
        '''
        This method plots mean flow for a single mesh
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                name = self.NetworkGraph.Edges[self.NetworkMesh.MeshToGraph[meshid]].Name
                dofs = element.GetPoiseuilleDofs()
                
                p1 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                
                if mean(p1)>=mean(p2):
                    flowDirection = self.flowDirection['(+)']
                    Flow = (p1-p2)/element.R
                else:
                    flowDirection = self.flowDirection['(-)']
                    Flow = (p2-p1)/element.R
                
                Radius = mean(element.Radius)
                print element.Name, Radius
                Velocity = Flow/(pi*Radius**2)
                
        if flowDirection == self.flowDirection['(+)']:
            plot(self.t, Velocity, 'r-',linewidth = 3, label = 'Velocity (+)')   #red line
        if flowDirection == self.flowDirection['(-)']:
            plot(self.t, Velocity, 'r-',linewidth = 3, label = 'Velocity (-)')   #red line
        xlabel('Time ($s$)')
        ylabel('Velocity ($m^3/s$)')
        title ('Velocity')    
        legend()
        savefig(self.images + meshid + '_' + name +'_vel.png')
        close()
       
    def PlotFlowComparative(self):
        '''
        This method plots brachial, radial and ulnar mean flow.
        If cycle is not specified, default cycle is the last one
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        FlowBrachial = 0
        FlowRadial = 0
        FlowUlnar = 0
        i = 0
        j = 0
        k = 0
                    
        for element in self.NetworkMesh.Elements:
            if element.Type is not 'Windkessel' and element.Name.find('brachial') != -1:
                dofs = element.GetPoiseuilleDofs()
                p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                if mean(p1)>=mean(p2):
                    FlowBrachial += (p1-p2)/element.R
                else:
                    FlowBrachial += (p2-p1)/element.R
                i += 1               
        
            if element.Type is not 'Windkessel' and element.Name.find('radial') != -1:   
                dofs = element.GetPoiseuilleDofs()
                p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                if mean(p1)>=mean(p2):
                    FlowRadial += (p1-p2)/element.R
                else:
                    FlowRadial += (p2-p1)/element.R
                j += 1
                
            if element.Type is not 'Windkessel' and element.Name.find('ulnar') != -1: 
                dofs = element.GetPoiseuilleDofs()
                p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                if mean(p1)>=mean(p2):
                    FlowUlnar += (p1-p2)/element.R
                else:
                    FlowUlnar += (p2-p1)/element.R
                k += 1

        if i != 0:
            FlowBrachial = (FlowBrachial*6e7)/i
            print "Brachial Flow = ", mean(FlowBrachial), "mL/min"
            plot(self.t, FlowBrachial, 'r-', linewidth = 3, label = 'brachial')   #red line
        if i == 0:
            FlowBrachial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if j != 0:
            FlowRadial = (FlowRadial*6e7)/j 
            print "Radial Flow = ", mean(FlowRadial), "mL/min"
            plot(self.t, FlowRadial, 'b-', linewidth = 3, label = 'radial')   #blue line
        if j == 0:
            FlowRadial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if k != 0:
            FlowUlnar = (FlowUlnar*6e7)/k
            print "Ulnar Flow = " , mean(FlowUlnar), "mL/min"
            plot(self.t, FlowUlnar, 'g-', linewidth = 3, label = 'ulnar')   #green line
        if k == 0:
            FlowUlnar = 0
            print "There isn't any element named brachial, check your network xml file"       

        xlabel('Time ($s$)')
        ylabel('Flow ($mL/min$)')
        title ('Brachial, Radial and Ulnar Flow Output')    
        legend()
        savefig(self.images+'brach_rad_uln_flow.png')
        close()
    
    def PlotFlows(self, meshIds):
        '''
        This method plots different flow volume signals in the same figure.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
            
        for meshid in meshIds:
            meshid = str(meshid)
            
            for element in self.NetworkMesh.Elements:
                if element.Id == meshid:
                    
                    dofs = element.GetPoiseuilleDofs()
                    p1 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                    p2 = self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                    
                    if mean(p1)>=mean(p2):
                        flowDirection = self.flowDirection['(+)']                
                        Flow = (p1-p2)/element.R
                    else:
                        flowDirection = self.flowDirection['(-)']
                        Flow = (p2-p1)/element.R
                    
                    if flowDirection == self.flowDirection['(+)']:                
                        plot(self.t, Flow*6e7, colourvector[indexcolour],linewidth = 3, label = element.Name)   
                    if flowDirection == self.flowDirection['(-)']:
                        plot(self.t, Flow*6e7, '--', colourvector[indexcolour],linewidth = 3, label = element.Name)   
                    
                    minY = 0
                    for q in Flow*6e7:
                        if q < minY:
                            minY = q
                    if minY != 0:
                        plot(self.t, zeros(len(Flow)),':',linewidth = 1)
                    ylim(ymin=minY)
                    
                    xlabel('Time ($s$)')
                    ylabel('Volumetric flow rate ($mL/min$)')
                    title ('Flows')   
                    legend(loc=0)
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
  
        savefig(self.f_images + ' Flows.png')
        close()
        indexcolour = 0
    
    def GetMeanFlow(self, el):
        '''
        This method returns flow signal for specific mesh
        '''
        dofs = el.GetPoiseuilleDofs()
        p1 = self.Solutions[(self.DofMap.DofMap[el.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
        p2 = self.Solutions[(self.DofMap.DofMap[el.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
        
        if mean(p1)>=mean(p2):                
            Flow = (p1-p2)/el.R
        else:
            Flow = (p2-p1)/el.R
        
        return mean(Flow) #m3/s
    
    def GetFlowSignal(self, el):
        '''
        This method returns flow signal for specific mesh
        '''
        dofs = el.GetPoiseuilleDofs()
        p1 = self.Solutions[(self.DofMap.DofMap[el.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
        p2 = self.Solutions[(self.DofMap.DofMap[el.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
        
        if mean(p1)>=mean(p2):               
            Flow = (p1-p2)/el.R
        else:
            Flow = (p2-p1)/el.R
        
        return Flow
    
    def WriteFlowTot(self, txtpath):
        '''
        This method writes in a txt file total flow over the network
        '''
        i = 0
        Flow = 0
        for element in self.NetworkMesh.Elements:
            dofs = element.GetPoiseuilleDofs()
            p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
            p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
            
            if mean(p1)>=mean(p2):                        
                Flow += (p1-p2)/element.R
            else:
                Flow += (p2-p1)/element.R            
            i+=1
            
        Flow = Flow*6e7 
        text_file = open(txtpath, "w")
        text_file.write("Flow Output (mL/min)\n")
        for word in Flow:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close()          
         
    def WriteFlowOutput(self, meshid, txtpath):
        '''
        This method writes flow output values in a .txt file.
        '''
        
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                if mean(p1)>=mean(p2):                        
                    Flow = (p1-p2)/element.R
                else:
                    Flow = (p2-p1)/element.R
                Flow = Flow*6e7     
                minrad = str(round((min(element.Radius))*1000,10))+' mm'
                maxrad = str(round((max(element.Radius))*1000,10))+' mm'
        text_file = open(txtpath, "w")
        text_file.write("Flow Output (mL/min)\n" +  "Min Radius: " +minrad +"\n" + "Max Radius: " + maxrad+ "\n")
        for word in Flow:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close()         
    
    def WriteReynolds(self, meshid, txtpath):
        '''
        This method writes reynolds number output values in a .txt file.
        '''
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                if mean(p1)>=mean(p2):                        
                    Flow = (p1-p2)/element.R
                else:
                    Flow = (p2-p1)/element.R
                Radius = element.Radius
                Reynolds = (2.0*Flow*self.SimulationContext.Context['blood_density'])/(pi*max(Radius)*self.SimulationContext.Context['dynamic_viscosity'])
                    
        text_file = open(txtpath, "w")
        text_file.write("Reynolds Number\n")
        for word in Reynolds:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close() 
    
    def PlotReynolds(self, meshid):
        '''
        This method plots reynolds number for a single mesh
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                name = self.NetworkGraph.Edges[self.NetworkMesh.MeshToGraph[meshid]].Name
                dofs = element.GetPoiseuilleDofs()
                p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                if mean(p1)>=mean(p2):                        
                    Flow = (p1-p2)/element.R
                else:
                    Flow = (p2-p1)/element.R

                Radius = element.Radius
                Reynolds = (2.0*Flow*self.SimulationContext.Context['blood_density'])/(pi*mean(Radius)*self.SimulationContext.Context['dynamic_viscosity'])
      
        plot(self.t, Reynolds, 'r-',linewidth = 3, label = 'Reynolds Number')   #red line
        xlabel('Time ($s$)')
        ylabel('Reynolds Number')
        title ('Reynolds N.'+' peak:'+str(round(max(Reynolds),1))+' mean:'+str(round(mean(Reynolds),1))+' min:'+str(round(min(Reynolds),1)))    
        legend()
        savefig(self.images + meshid + '_' + name +'_reynoldsN.png')
        close()
        
    # PRESSURE METHODS
    
    def PlotPressure(self, meshid):
        '''
        This method plots pressures for a single mesh
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                name = self.NetworkGraph.Edges[self.NetworkMesh.MeshToGraph[meshid]].Name
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 0]),self.CardiacFreq*(self.Cycles-1):])
                Pressure2 = (self.Solutions[(self.DofMap.DofMap[meshid, 2]),self.CardiacFreq*(self.Cycles-1):])
                
                
        print name, round(mean(Pressure/133.32),1),"--",round(mean(Pressure2/133.32),1)        
        plot(self.t, Pressure/133.32, 'b-', linewidth = 3, label = 'Pressure Signal')   #blue line
        
        minY = 0
        for p in Pressure/133.32:
            if p < minY:
                minY = p
        
        if minY != 0:
            plot(self.t, zeros(len(Pressure)),':',linewidth = 1)
        
        ylim(ymin=minY)
        xlabel('Time ($s$)')
        ylabel('Pressure ($mmHg$)')
        title ('Pressure'+' peak:'+str(round(max(Pressure/133.32),1))+' mean:'+str(round(mean(Pressure/133.32),1))+' min:'+str(round(min(Pressure/133.32),1)))    
        legend()
        savefig(self.p_images + meshid + '_' + name +'_pressure.png')
        close()
    
    def PlotPressures(self, meshIds):
        '''
        This method plots different pressure signals in the same figure.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
            
        for meshid in meshIds:
            meshid = str(meshid)
            for element in self.NetworkMesh.Elements:
                if element.Id == meshid:
                    Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 0]),self.CardiacFreq*(self.Cycles-1):])
                    plot(self.t, Pressure/133.32, colourvector[indexcolour],linewidth = 3, label = element.Name)   
                    
                    minY = 0
                    for p in Pressure/133.32:
                        if p < minY:
                            minY = p
                    
                    if minY != 0:
                        plot(self.t, zeros(len(Pressure)),':',linewidth = 1)
                    
                    ylim(ymin=minY)
                    
                    xlabel('Time ($s$)')
                    ylabel('Pressure ($mmHg$)')
                    title ('Pressures')   
                    legend(loc=0)
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
  
        savefig(self.p_images + ' Pressures.png')
        close()
        indexcolour = 0
    
    
    def PlotPressureTwo(self, meshid, meshid2):
        '''
        This method plots pressures for a couple of meshes
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211
                print element.Name
        
        meshid2 = str(meshid2)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid2:
                Pressure2 = (self.Solutions[(self.DofMap.DofMap[meshid2, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211
                print element.Name
        
        plot(self.t, Pressure, 'b-', linewidth = 3, label = 'Pressure Signal')   #blue line
        plot(self.t, Pressure2, 'r-', linewidth = 3, label = 'Pressure Signal')   #red line
        xlabel('Time ($s$)')
        ylabel('Pressure ($mmHg$)') 
        title ('Pressure')    
        legend()
        savefig(self.images + meshid + meshid2 +'_pressure.png')
        close()
        
    def PlotPressureComparative(self):
        '''
        This method plots brachial, radial and ulnar pressure signal.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        PressureBrachial = 0
        PressureRadial = 0
        PressureUlnar = 0
        i = 0
        j = 0
        k = 0
        for element in self.NetworkMesh.Elements:
            if element.Type is not 'Windkessel' and element.Name.find('brachial') != -1:
                PressureBrachial += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211
                i += 1               
        
            if element.Type is not 'Windkessel' and element.Name.find('radial') != -1:
                PressureRadial += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211
                j += 1
                
            if element.Type is not 'Windkessel' and element.Name.find('ulnar') != -1:
                PressureUlnar += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211
                k += 1
        
        if i != 0:
            PressureBrachial = PressureBrachial/i
            plot(self.t, PressureBrachial, 'r-', linewidth = 3, label = 'brachial')   #red line
        if i == 0:
            PressureBrachial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if j != 0:
            PressureRadial = PressureRadial/j
            plot(self.t, PressureRadial, 'b-', linewidth = 3, label = 'radial')   #blue line
        if j == 0:
            PressureRadial = 0
            print "There isn't any element named radial, check your network xml file"
        if k != 0:
            PressureUlnar = PressureUlnar/k
            plot(self.t, PressureUlnar, 'g-', linewidth = 3, label = 'ulnar')   #green line
        if k == 0:
            PressureUlnar = 0
            print "There isn't any element named ulnar, check your network xml file"
        
        xlabel('Time ($s$)')
        ylabel('Pressure ($mmHg$)')
        title ('Brachial, Radial and Ulnar Pressure Signal')    
        legend()
        savefig(self.images+'brach_rad_uln_pressure.png')
        close()
        
    def GetPressureSignal(self, meshid):
        '''
        This method returns pressures signal for specific mesh
        '''
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dof = element.GetNodeLocalDofs()
                PressureIN = (self.Solutions[(self.DofMap.DofMap[meshid, dof[0]]),self.CardiacFreq*(self.Cycles-1):])
                PressureOUT = (self.Solutions[(self.DofMap.DofMap[meshid, dof[1]]),self.CardiacFreq*(self.Cycles-1):])

        return PressureIN,PressureOUT
    
    def PlotPressureDrop(self, meshid):
        '''
        This method returns pressure drop for specific mesh
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                length = element.Length
                dof = element.GetNodeLocalDofs()
                PressureIN = (self.Solutions[(self.DofMap.DofMap[meshid, dof[0]]),self.CardiacFreq*(self.Cycles-1):])
                PressureOUT = (self.Solutions[(self.DofMap.DofMap[meshid, dof[1]]),self.CardiacFreq*(self.Cycles-1):])
        
        plot(self.t, (PressureIN-PressureOUT), 'b-', linewidth = 3, label = 'p1-p2')   #blue line
        plot(self.t, (PressureIN-PressureOUT)/length, 'r-', linewidth = 3, label = '(p1-p2)/L')   #blue line
        xlabel('Time ($s$)')
        ylabel('P ($Pa$)')
        title ('Pressure')    
        legend()
        savefig(self.images + meshid +'_pressureDrop.png')
        close()
        
        return (PressureIN-PressureOUT)/length
    
    def WritePressureDrop(self, meshid, txtpath):
        '''
        This method writes pressure drop values in a .txt file.
        '''
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                length = element.Length
                dof = element.GetNodeLocalDofs()
                PressureDrop = ((self.Solutions[(self.DofMap.DofMap[meshid, dof[0]]),self.CardiacFreq*(self.Cycles-1):])-(self.Solutions[(self.DofMap.DofMap[meshid, dof[1]]),self.CardiacFreq*(self.Cycles-1):]))/length
                meanP = mean(PressureDrop)
        text_file = open(txtpath, "w")
        text_file.write("Pressure Drop over length (Pa)\n")
        for word in PressureDrop:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.write("MEAN\n")
        text_file.write(str(meanP))
        text_file.close()
        
    def WritePressureInput(self, meshid, txtpath):
        '''
        This method writes pressure output values in a .txt file.
        '''
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 0]),self.CardiacFreq*(self.Cycles-1):])
                meanP = mean(Pressure)

        text_file = open(txtpath, "w")
        text_file.write("Pressure Input (Pa)\n")
        for word in Pressure:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.write("MEAN\n")
        text_file.write(str(meanP))
        text_file.close()
    
    def WritePressureOutput(self, meshid,  txtpath):
        '''
        This method writes pressure output values in a .txt file.
        '''
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 2]),self.CardiacFreq*(self.Cycles-1):])
                meanP = mean(Pressure)
                
        text_file = open(txtpath, "w")
        text_file.write("Pressure Output (Pa)\n")
        for word in Pressure:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.write("MEAN\n")
        text_file.write(str(meanP))
        text_file.close()
        
    # Wall Shear Stress methods    
    def PlotPWSS(self, meshid):
        '''
        This method plots mean WSS (POISEUILLE) for a single mesh 
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                if mean(p1)>=mean(p2):                        
                    Flow = (p1-p2)/element.R
                else:
                    Flow = (p2-p1)/element.R
                
                Wss = ((4.0*element.eta)/pi) * (Flow/(mean(element.Radius))**3)
                print element.Name, "Wss(mean) = ", mean(Wss), "Pa\n", "Wss(peak) = ", max(Wss), "Pa"
                
        plot(self.t, Wss, 'r-',linewidth = 3, label = 'wss')   #red line
        xlabel('Time ($s$)')
        ylabel('Wss ($Pa$)')
        title ('Wall Shear Stress')    
        legend()
        savefig(self.w_images + meshid +'_Pwss.png')
        close()
    
    def PlotPWSSs(self, meshIds):
        '''
        This method plots different poiseuille's wss signals in the same figure.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
            
        for meshid in meshIds:
            meshid = str(meshid)
            for element in self.NetworkMesh.Elements:
                if element.Id == meshid:
                    dofs = element.GetPoiseuilleDofs()
                    p1 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):]
                    p2 = self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):]
                    if mean(p1)>=mean(p2):                        
                        Flow = (p1-p2)/element.R
                    else:
                        Flow = (p2-p1)/element.R

                    Wss = ((4.0*element.eta)/pi) * (Flow/(mean(element.Radius))**3)
                    plot(self.t, Wss, colourvector[indexcolour],linewidth = 3, label = element.Name)   
                    xlabel('Time ($s$)')
                    ylabel('Wss ($Pa$)')
                    title ('Wall shear stress')   
                    legend(loc=0)
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
  
        savefig(self.w_images + ' Wss.png')
        close()
        indexcolour = 0
    
        
    def PlotPWSSComparative(self):
        '''
        This method plots brachial, radial and ulnar WSS (POISEUILLE) signal.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        FlowBrachial = 0
        FlowRadial = 0
        FlowUlnar = 0
        WSSBrachial = 0
        WSSRadial = 0
        WSSUlnar = 0
        i = 0
        j = 0
        k = 0
        for element in self.NetworkMesh.Elements:
            if element.Type is not 'Windkessel' and element.Name.find('rachial') != -1:
                dofs = element.GetPoiseuilleDofs()
                FlowBrachial = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                WSSBrachial += ((4.0*element.eta)/pi) * (FlowBrachial/(mean(element.Radius))**3)
                i += 1               
        
            if element.Type is not 'Windkessel' and element.Name.find('adial') != -1:
                dofs = element.GetPoiseuilleDofs()
                FlowRadial = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                WSSRadial += ((4.0*element.eta)/pi) * (FlowRadial/(mean(element.Radius))**3)
                j += 1
                
            if element.Type is not 'Windkessel' and element.Name.find('lnar') != -1:
                dofs = element.GetPoiseuilleDofs()
                FlowUlnar = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                WSSUlnar += ((4.0*element.eta)/pi) * (FlowUlnar/(mean(element.Radius))**3)
                k += 1
        
        if i != 0:
            WSSBrachial = WSSBrachial/i
            print "Brachial"
            print "Wss(mean) = ", mean(WSSBrachial), " Pa", " Wss(peak) = ", max(WSSBrachial), " Pa"
            plot(self.t, WSSBrachial, 'r-', linewidth = 3, label = 'brachial')   #red line
        if i == 0:
            WSSBrachial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if j != 0:
            WSSRadial = WSSRadial/j
            print "Radial"
            print "Wss(mean) = ", mean(WSSRadial), " Pa", " Wss(peak) = ", max(WSSRadial), " Pa"
            plot(self.t, WSSRadial, 'b-', linewidth = 3, label = 'radial')   #blue line
        if j == 0:
            WSSRadial = 0
            print "There isn't any element named radial, check your network xml file"
        if k != 0:
            WSSUlnar = WSSUlnar/k
            print "Ulnar"
            print "Wss(mean) = ", mean(WSSUlnar), " Pa", " Wss(peak) = ", max(WSSUlnar), " Pa"
            plot(self.t, WSSUlnar, 'g-', linewidth = 3, label = 'ulnar')   #green line
        if k == 0:
            WSSUlnar = 0
            print "There isn't any element named ulnar, check your network xml file"
        
        xlabel('Time ($s$)')
        ylabel('Wall Shear Stress ($Pa$)')
        title ('Brachial, Radial and Ulnar Wall Shear Stress')    
        legend()
        savefig(self.images+'brach_rad_uln_wss.png')    
        close()
        
    def GetPWSSSignal(self, meshid):
        '''
        This method returns Wall Shear Stress signal (POISEUILLE) for specific mesh
        '''
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                WSS = ((4.0*element.eta)/pi) * (Flow/(mean(element.Radius))**3)
                
        return WSS
        
    
    def WritePWSSOutput(self, meshid, txtpath):
        '''
        This method writes Wall Shear Stress output values (POISEUILLE) in a .txt file
        '''
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                WSS = ((4.0*element.eta)/pi) * (Flow/(mean(element.Radius))**3)
        
        text_file = open(txtpath, "w")
        text_file.write("Wall Shear Stress Output (Pa)\n")
        for word in WSS:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close()
        
    def WriteWSSOutput(self, el, txtpath):
        '''
        This method writes Wall Shear Stress output values (POISEUILLE) in a .txt file
        If cycle is not specified, default cycle is the last one.
        '''     
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        
        inverseWomersley.SetFlowSignal(el, self.GetFlowSignal(el))
        inverseWomersley.GetTaoFromQ(el)
        
        text_file = open(txtpath, "w")
        text_file.write("Wall Shear Stress Output (dyne/cm2)\n")
        for word in inverseWomersley.Tauplot:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close() 
        
    
    def GetWssPeak(self, mesh):
        '''
        This method returns Wall Shear Stress peak value for specific mesh
        '''
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        rPeaks = inverseWomersley.GetWssPeaks(mesh, self.GetFlowSignal(mesh))
        return rPeaks
       
    def ShowVelocityProfile(self, el):
        '''
        This method show velocity profile signal (over fractional radius) using WX python library.
        '''
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        inverseWomersley.SetFlowSignal(el, self.GetFlowSignal(el))
        inverseWomersley.GetVelFromQ(el)
        inverseWomersley.ShowVelocityProfile(el.Id)
        
    def SaveVelocityProfile(self, el, daystr):
        '''
        This method save velocity profile signal (over fractional radius) in a .avi file using menCoder library.
        '''
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        inverseWomersley.SetFlowSignal(el, self.GetFlowSignal(el))
        inverseWomersley.GetVelFromQ(el)
        inverseWomersley.SaveVelocityProfile(el.Id, daystr)
        
    def PlotWSS(self, el):
        '''
        This method plots mean WSS for a single mesh 
        If cycle is not specified, default cycle is the last one 
        ''' 
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        
        inverseWomersley.SetFlowSignal(self.GetFlowSignal(el))
        inverseWomersley.GetTaoFromQ(el)
        peak = inverseWomersley.PlotWss(el.Id, self.w_images)
        
        self.dayWssP[el.Id] = peak
        
    def PlotManyWSS(self, els):
        '''
        This method plots different womersley wss signals ino the same figure.
        ''' 
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close, ylim
        except ImportError:
            MatPlotLibError()
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
        
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        
        for el in els:
        
            inverseWomersley.SetFlowSignal(el, self.GetFlowSignal(el))
            inverseWomersley.GetTaoFromQ(el)
            
            tplot = arange(0,inverseWomersley.tPeriod,inverseWomersley.dtPlot)
            plot(tplot, inverseWomersley.Tauplot, colourvector[indexcolour] ,linewidth = 3, label = el.Name)
            minY = 0
            for w in inverseWomersley.Tauplot:
                if w < minY:
                    minY = w
        
            if minY != 0:
                plot(tplot, zeros(len(inverseWomersley.Tauplot)),':',linewidth = 1)
                
            ylim(ymin=minY)
            
            xlabel('Time ($s$)')
            ylabel('Wall shear stress ($dyne/cm^2$)')
            title ('Wss')    
            legend(loc=0)
            if indexcolour == 6:
                indexcolour = 0
            else:
                indexcolour+=1
            
            
        savefig(self.w_images + ' Wss.png')   
        close()
        indexcolour = 0
        
        
    def GetWSSSignal(self, el):
        '''
        This method returns Wall Shear Stress signal for specific mesh
        '''
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        
        inverseWomersley.SetFlowSignal(el, self.GetFlowSignal(el))
        WssSignal = array(inverseWomersley.GetTaoFromQ(el))
    
        return WssSignal
    
    #OTHER METHODS
    
    def PulseWaveVelocity(self, superEdge, superEdge2 = None):
        '''
        This method computes Pulse Wave Velocity.(m/s)
        If SuperEdge2 is specified, PWV is computed between first and second superedge,
        otherwise PWV is computed over a single superedge.
        '''
        distance = 0 
        startSE = str(superEdge)
        if superEdge2 == None:
            endSE = str(superEdge)
        else:
            endSE = str(superEdge2)
        if startSE == endSE:
            for se in self.NetworkGraph.SuperEdges.itervalues():
                if se.Name == startSE:
                    if len(se.Edges) == 1:
                        startE = se.Edges.values()[0]
                        endE = startE
                        nodoStart = startE.NodeIds[0]
                        nodoEnd = startE.NodeIds[1]
                        distance = startE.Length['value']
                    else:
                        for ed in se.Edges.itervalues():
                            node0 = ed.NodeIds[0]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[1] == node0:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoStart = node0
                        for ed in se.Edges.itervalues():
                            distance+=ed.Length['value']
                            if nodoStart == ed.NodeIds[0]:
                                startE = ed                
                        for ed in se.Edges.itervalues():
                            node1 = ed.NodeIds[1]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[0] == node1:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoEnd = node1
                        for ed in se.Edges.itervalues():
                            if nodoEnd == ed.NodeIds[1]:
                                endE = ed    
        else: 
            for se in self.NetworkGraph.SuperEdges.itervalues():
                if se.Name == startSE:
                    if len(se.Edges) == 1:
                        startE = se.Edges.values()[0]
                        nodoStart = startE.NodeIds[0]
                        distance+=startE.Length['value']
                    else:
                        for ed in se.Edges.itervalues():
                            node0 = ed.NodeIds[0]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[1] == node0:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoStart = node0
                        for ed in se.Edges.itervalues():
                            distance+=ed.Length['value']
                            if nodoStart == ed.NodeIds[0]:
                                startE = ed                             
                if se.Name == endSE:
                    if len(se.Edges) == 1:
                        endE = se.Edges.values()[0]
                        nodoEnd = endE.NodeIds[1]
                        distance+=endE.Length['value']
                    else:
                        for ed in se.Edges.itervalues():
                            node1 = ed.NodeIds[1]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[0] == node1:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoEnd = node1
                        for ed in se.Edges.itervalues():
                            distance+=ed.Length['value']
                            if nodoEnd == ed.NodeIds[1]:
                                endE = ed        
        #mesh
        meshNodeStart = self.NetworkMesh.meshToEdges[nodoStart]
        meshNodeEnd = self.NetworkMesh.meshToEdges[nodoEnd]
        for el in self.NetworkMesh.Elements:
            if el.NodeIds[0] == meshNodeStart:
                meshStart = el
            if el.NodeIds[1] == meshNodeEnd:
                meshEnd = el
        #finding peak pressure and corresponding time.
        PeakStartAll = (self.Solutions[(self.DofMap.DofMap[meshStart.Id, 0]),self.CardiacFreq*(self.Cycles-1):])
        PeakStart = max(PeakStartAll)   
        step1 = 0
        for p1 in PeakStartAll:
            if p1 == PeakStart:
                t1 = self.t[step1]
                break
            step1+=1        
        PeakEndAll = (self.Solutions[(self.DofMap.DofMap[meshEnd.Id, 1]),self.CardiacFreq*(self.Cycles-1):])
        PeakEnd = max(PeakEndAll)      
        step2 = 0
        for p2 in PeakEndAll:
            if p2 == PeakEnd:
                t2 = self.t[step2]
                break
            step2+=1        
        deltat = abs(t1-t2)        
        PWV = float(distance)/deltat       
        print PWV, 'm/s'
        return PWV
     
    # ENTITY METHOD
    
    def GetSolution(self, entity):
        '''
        This method gets and plots flow, pressure and WSS for specific entity.
        '''
        try:
            from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
        except ImportError:
            MatPlotLibError()
        Flow = 0
        WSS = 0
        WSSW = 0
        PressureIN = 0
        PressureOUT = 0
        num_el = 0
        
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find(entity) != -1:
                for element in el:
                    num_el+=1
                    dofs = element.GetPoiseuilleDofs()
                    Flow += (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/element.R
                    WSS += ((4.0*element.eta)/pi) * (Flow/pow(mean(element.Radius), 3))
                    WSSW += self.GetWSSSignal(element)
                    PressureIN += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),self.CardiacFreq*(self.Cycles-1):])   
                    PressureOUT += (self.Solutions[(self.DofMap.DofMap[element.Id, 2]),self.CardiacFreq*(self.Cycles-1):])
                
        PressureIN = ((PressureIN)/num_el)/133.3223684211
        PressureOUT = ((PressureOUT)/num_el)/133.3223684211           
        print entity, " PressureIN = ", mean(PressureIN), "mmHg"
        print entity, " PressureOUT = ", mean(PressureOUT), "mmHg"          
        plot(self.t, PressureIN, 'b-', linewidth = 3, label = 'Pressure IN')   #blue line            
        plot(self.t, PressureOUT, 'r-', linewidth = 3, label = 'Pressure OUT')   #red line   
        xlabel('Time ($s$)')
        ylabel('Pressure ($mmHg$)')
        title (entity + ' Pressure Signal')    
        legend()
        savefig(self.images+entity+'_pressure.png')
        close() 
                   
        Flow = (Flow*6e7)/num_el
        print entity, " Flow = ", mean(Flow), "mL/min"
        plot(self.t, Flow, 'r-', linewidth = 3, label = entity)   #red line
        xlabel('Time ($s$)')
        ylabel('Flow ($mL/min$)')
        title (entity + ' Flow Output')    
        legend()
        savefig(self.images+entity+'_flow.png')
        close()
        
        WSS = WSS/num_el
        WSSW = WSSW/num_el
        print entity, "Wss(mean) = ", mean(WSS), " Pa", " Wss(peak) = ", max(WSSW)
        plot(self.t, WSS, 'g-', linewidth = 3, label = entity)   #green line
        xlabel('Time ($s$)')
        ylabel('Wall Shear Stress ($Pa$)')
        title (entity + ' WSS Output')    
        legend()
        savefig(self.images+entity+'_wss.png')
        close()               
    
    # SOLUTIONS IN XML NETWORK GRAPH FORMAT
    
    def WriteToXML(self, xmlsolutionspath):
        '''
        This method writes solutions in XML MeshSolutions File
        '''
        print "Writing xml Solution file..."
        root = etree.Element("Solutions", id=self.NetworkGraph.Id, version="2.0")
        xmlsolutions = etree.ElementTree(root)
        
        #CASE
        case = etree.SubElement(root, "case")
        patId = etree.SubElement(case, "patient_id")
        patId.text = self.NetworkGraph.PatientId
        visit = etree.SubElement(case, "visit")
        visit.text = self.NetworkGraph.Visit
        
        #NODES
        nodes_list = []
        nodes = etree.SubElement(root, "nodes")
        for node in self.NetworkGraph.Nodes.itervalues():
            nodes_list.append(int(node.Id))
        nodes_list.sort()
        
        for id in nodes_list:
            for nodeG in self.NetworkGraph.Nodes.itervalues():
                if int(nodeG.Id) == id:
                    if nodeG.Name:
                        n_name = nodeG.Name
                        n_type = nodeG.Type
                        node = etree.SubElement(nodes, "node", id = str(id), type = str(n_type), name = str(n_name))
                        if nodeG.Name == 'anastomosis':
                            prop = etree.SubElement(node, "properties")
                            conn = etree.SubElement(prop, "connections")
                            etree.SubElement(conn, "proximal_artery", edge_id=str(nodeG.Properties['proximal']))
                            try:
                                etree.SubElement(conn, "distal_artery", edge_id=str(nodeG.Properties['distal']))
                            except:
                                pass
                            etree.SubElement(conn, "proximal_vein", edge_id=str(nodeG.Properties['vein']))
                    else:
                        etree.SubElement(nodes, "node", id = str(id))
            
        #SUPEREDGES
        superedges_list = []
        superedges = etree.SubElement(root, "superedges")
        for sedges in self.NetworkGraph.SuperEdges.iterkeys():
            superedges_list.append(int(sedges))
        superedges_list.sort()
        
        for sedg in superedges_list:
            for s in self.NetworkGraph.SuperEdges.itervalues():
                if s.Id == str(sedg):
                    if s.SuperEdges != {}:
                        superedge = etree.SubElement(superedges, "superedge", id = str(s.Id), name = str(s.Name))
                        superedges2 = etree.SubElement(superedge, "superedges")
                    if s.SuperEdges == {}:
                        try:
                            superedge2 = etree.SubElement(superedges2,"superedge", id = str(s.Id), name = str(s.Name))
                        except:
                            superedge2 = etree.SubElement(superedges,"superedge", id = str(s.Id), name = str(s.Name))
                        edgeIdsel = etree.SubElement(superedge2, "edgesIds")
                        for edgeIds in s.Edges.iterkeys():
                            etree.SubElement(edgeIdsel, "edgeIds", edge_id = str(edgeIds))
                            
        #EDGES
        edges_list = []
        edges = etree.SubElement(root, "edges")
        for edge in self.NetworkGraph.Edges.iterkeys():
            edges_list.append(int(edge))
        edges_list.sort()
        
        for edg in edges_list:
            for e in self.NetworkGraph.Edges.itervalues():
                if e.Id == str(edg):
                    edge = etree.SubElement(edges, "edge", id = str(e.Id), name = str(e.Name), side = str(e.Side), node1_id = str(e.NodeIds[0]), node2_id = str(e.NodeIds[1]))
                    
                    Flow = 0
                    i = 0                 
                    for el in self.NetworkMesh.Elements:
                        if el.NodeIds[0] == self.NetworkMesh.s_mesh[(0.0,e)]:
                            P1 = (self.Solutions[(self.DofMap.DofMap[el.Id, 0]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211
                            P1Mean = mean(P1)
                            P1Max = max(P1)
                            P1Min = min(P1)         
                            break
                    for el in self.NetworkMesh.Elements:
                        if el.NodeIds[1] == self.NetworkMesh.s_mesh[(1.0,e)]:
                            P2 = (self.Solutions[(self.DofMap.DofMap[el.Id, 2]),self.CardiacFreq*(self.Cycles-1):])/133.3223684211
                            P2Mean = mean(P2)
                            P2Max = max(P2)     
                            P2Min = min(P2)           
                            break                
                    for meshId, edgeId in self.NetworkMesh.MeshToGraph.iteritems():
                        if edgeId == e.Id:
                            for el in self.NetworkMesh.Elements:
                                if str(el.Id) == str(meshId):
                                    dofs = el.GetPoiseuilleDofs()
                                    Flow += (self.Solutions[(self.DofMap.DofMap[el.Id, dofs[0]]),self.CardiacFreq*(self.Cycles-1):] - self.Solutions[(self.DofMap.DofMap[el.Id, dofs[1]]),self.CardiacFreq*(self.Cycles-1):])/el.R
                                    Wss = ((4.0*el.eta)/pi) * (Flow/mean(el.Radius)**3)
                                    WssMax = max(self.GetWssPeak(el))
                                    i += 1
                                         
                    Flow = (Flow*6e7) / i
                    Wss = (Wss) / i
                    FlowMean = mean(Flow)
                    FlowMax = max(Flow)
                    FlowMin = min(Flow)
                    WssMean = mean(Wss)*10
                    WssMax = WssMax*10
                    WssMin = min(Wss)*10
                    
                    solution = etree.SubElement(edge, "solution")
                    
                    solPmean = etree.SubElement(solution, "pressure_timemean", unit = "mmHg")
                    solPmean_s1 = etree.SubElement(solPmean, "value", s="0.0")
                    solPmean_s1_v = etree.SubElement(solPmean_s1, "scalar")
                    solPmean_s1_v.text = str(P1Mean)
                    solPmean_s2 = etree.SubElement(solPmean, "value", s="1.0")
                    solPmean_s2_v = etree.SubElement(solPmean_s2, "scalar")
                    solPmean_s2_v.text = str(P2Mean)
                    
                    solPmax = etree.SubElement(solution, "pressure_timemax", unit = "mmHg")
                    solPmax_s1 = etree.SubElement(solPmax, "value", s="0.0")
                    solPmax_s1_v = etree.SubElement(solPmax_s1, "scalar")
                    solPmax_s1_v.text = str(P1Max)
                    solPmax_s2 = etree.SubElement(solPmax, "value", s="1.0")
                    solPmax_s2_v = etree.SubElement(solPmax_s2, "scalar")
                    solPmax_s2_v.text = str(P2Max)
                    
                    solPmin = etree.SubElement(solution, "pressure_timemin", unit = "mmHg")
                    solPmin_s1 = etree.SubElement(solPmin, "value", s="0.0")
                    solPmin_s1_v = etree.SubElement(solPmin_s1, "scalar")
                    solPmin_s1_v.text = str(P1Min)
                    solPmin_s2 = etree.SubElement(solPmin, "value", s="1.0")
                    solPmin_s2_v = etree.SubElement(solPmin_s2, "scalar")
                    solPmin_s2_v.text = str(P2Min)
                    
                    solQmean = etree.SubElement(solution, "flow_mean", unit = "mL/min")
                    solQmean_value = etree.SubElement(solQmean, "scalar")
                    solQmean_value.text = str(FlowMean)
                    
                    solQmax = etree.SubElement(solution, "flow_max", unit = "mL/min")
                    solQmax_value = etree.SubElement(solQmax, "scalar")
                    solQmax_value.text = str(FlowMax)
                    
                    solQmin = etree.SubElement(solution, "flow_min", unit = "mL/min")
                    solQmin_value = etree.SubElement(solQmin, "scalar")
                    solQmin_value.text = str(FlowMin)
                         
                    solWssmean = etree.SubElement(solution, "wss_mean", unit = "dyne/cm2")
                    solWssmean_value = etree.SubElement(solWssmean, "scalar")
                    solWssmean_value.text = str(WssMean)
                    
                    solWssmax = etree.SubElement(solution, "wss_max", unit = "dyne/cm2")
                    solWssmax_value = etree.SubElement(solWssmax, "scalar")
                    solWssmax_value.text = str(WssMax)
                    
                    solWssmin = etree.SubElement(solution, "wss_min", unit = "dyne/cm2")
                    solWssmin_value = etree.SubElement(solWssmin, "scalar")
                    solWssmin_value.text = str(WssMin)
        
        indent(root)            
        xmlsolutions.write (xmlsolutionspath,encoding='iso-8859-1')
        
    def WriteToCsv(self, adaptation, solutionType, cycle = None):
        '''
        This method writes in .csv files mean values for diameters, flows, pressures and wss.
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles
        
        path = solutionType+'.csv' 
        ofile  = open(path, "wb")
        writer = csv.writer(ofile, delimiter=",", quoting=csv.QUOTE_ALL)
        header_list = []
        header_list.append("Id")
        header_list.append("Name")
        header_list.append("IN/OUT")
        
        for d in adaptation.solutions.keys():
            if d != -1:
                header_list.append(d)
        header = [[],header_list]
        
        writer.writerows(header)
            
        if solutionType == 'Pressure':
            print "Writing Pressures Csv file..."
            for el in self.NetworkMesh.Elements:
                if el.Type ==  'WavePropagation':
                    el_row_list_in = [el.Id,el.Name,"IN"]
                    el_row_list_out = [el.Id,el.Name,"OUT"]
                    for day,sol in adaptation.solutions.iteritems():
                        if day != -1:
                            P1 = (sol.Solutions[(self.DofMap.DofMap[el.Id, 0]),:])/133.3223684211
                            P1Mean = mean(P1[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))])
                            P2 = (sol.Solutions[(self.DofMap.DofMap[el.Id, 2]),:])/133.3223684211
                            P2Mean = mean(P2[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))])
                            el_row_list_in.append(str(P1Mean))
                            el_row_list_out.append(str(P2Mean))
                    el_row = [el_row_list_in,el_row_list_out]
                    writer.writerows(el_row)
        
        if solutionType == 'Diameter':
            print "Writing Diameters Csv file..."
            for el in self.NetworkMesh.Elements:
                if el.Type ==  'WavePropagation':
                    el_row_list_in = [el.Id,el.Name,"IN"]
                    el_row_list_out = [el.Id,el.Name,"OUT"]
                    for day,sol in adaptation.solutions.iteritems():
                        if day != -1:
                            try:
                                d1 = el.dayRadius[day][0]*2e3
                                d2 = el.dayRadius[day][1]*2e3
                            except:
                                d1 = el.Radius[0]*2e3
                                d2 = el.Radius[len(el.Radius)-1]*2e3
                            el_row_list_in.append(str(d1))
                            el_row_list_out.append(str(d2))
                    el_row = [el_row_list_in,el_row_list_out]
                    writer.writerows(el_row)
                
        if solutionType == 'Flow':
            print "Writing Flows Csv file..."
            for el in self.NetworkMesh.Elements:
                if el.Type ==  'WavePropagation':
                    el_row_list = [el.Id,el.Name,""]
                    for day,sol in adaptation.solutions.iteritems():
                        if day != -1:
                            Flow = sol.dayFlow[el.Id]
                            el_row_list.append(str(Flow))
                    el_row = [el_row_list]
                    writer.writerows(el_row)
        
        if solutionType == 'Wss':
            print "Writing WssPeaks Csv file..."
            for el in self.NetworkMesh.Elements:
                if el.Type ==  'WavePropagation':
                    el_row_list = [el.Id,el.Name,""]
                    for day,sol in adaptation.solutions.iteritems():
                        if day != -1:
                            tao = sol.dayWssP[el.Id]
                            el_row_list.append(str(tao))
                    el_row = [el_row_list]
                    writer.writerows(el_row)

class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass

class MatPlotLibError(Error):
    '''
    Exception raised if matplotlib package is not installed.
    '''
    def __init__(self):
        sys.exit('Matplotlib package is required for plotting solutions in .png files.\nPlease download matplotlib from matplotlib.sourceforge.net.')        
              
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
