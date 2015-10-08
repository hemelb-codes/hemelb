#!/usr/bin/env python

## Program:   PyNS
## Module:    Elements.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from numpy.core.numeric import arange, array, zeros, dot
from numpy.numarray.numerictypes import Int32
from math import pi
from numpy.core.fromnumeric import sum, mean
from numpy.lib.scimath import sqrt

class Element(object):
    '''
    A general Element, each specific element must be referred to it.
    This class provides the following methods:
    Initialize: a method for setting parameters for simulation from SimulationContext.
    IsNonLinear: a method which returns true if the current element is a non linear one, otherwise it returns false.
    GetZeroOrderMatrix, GetFirstOrderMatrix, GetSecondOrderMatrix -> this methods calculate local zero, first or second order matrix.
    GetNumberOfNodes: a method for calculating and returning element's number of nodes.
    GetNumberOfDofs: a method for calculating and returning element's number of dofs (degrees of freedom).
    SetParameterInHistory: a methodP for setting parameter History dictionary. Each value (for each parameters name) is associated to different timesteps
    RegisterElementType: a method for registering element's class according to its type.
    NewElement: a method for creating a new element according to its type.
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.simulationContext = None
        self.nonLinear = False
        self.ParameterInfo = {}
        self.LastIncrementNumber = 0
        self.Type = None
               
    def Initialize(self, simulationContext):
        '''
        Setting SimulationContext.
        '''
        self.simulationContext = simulationContext
                 
    def InputParameters(self):
        '''
        General method for parameters input.
        '''
        pass

    def IsNonLinear(self):
        '''
        This method returns True if the current element
        is a non linear one, otherwise it returns false.
        '''
        return self.nonLinear

    def GetZeroOrderMatrix(self):
        '''
        This method returns element's local zero order matrix
        '''
        circuitMatrix = self.GetCircuitMatrix()      
        numberOfEdges = circuitMatrix.shape[0]
        elementNodes = array(circuitMatrix[:,:2].astype(Int32))
        elementData = circuitMatrix[:,2:]              
        zeroOrderMatrix = zeros((numberOfEdges+1, numberOfEdges+1))
        reverse = array([[1,-1],[-1,1]])
        for i in arange(0,numberOfEdges):
            inductance = elementData[i,2]
            if inductance > 0.0:
                inductance = 1.0/inductance
            intermediateMatrix = dot(inductance,reverse)
            ind1 = array([elementNodes[i,:],elementNodes[i,:]])
            ind2 = ind1.transpose()
            zeroOrderMatrix[ind2, ind1] = zeroOrderMatrix[ind2, ind1] + intermediateMatrix
        self.zeroOrderMatrix = zeroOrderMatrix
        return zeroOrderMatrix
        
    def GetFirstOrderMatrix (self):
        '''
        This method returns element's local first order matrix
        '''
        circuitMatrix = self.GetCircuitMatrix()       
        numberOfEdges = circuitMatrix.shape[0]
        elementNodes = array(circuitMatrix[:,:2].astype(Int32))
        elementData = circuitMatrix[:,2:]              
        firstOrderMatrix = zeros((numberOfEdges+1, numberOfEdges+1))
        reverse = array([[1,-1],[-1,1]])
        for i in arange(0,numberOfEdges):
            resistance = elementData[i,1]
            if resistance > 0.0:
                resistance = 1.0/resistance
            intermediateMatrix = dot(resistance,reverse)  
            ind1 = array([elementNodes[i,:],elementNodes[i,:]])
            ind2 = ind1.transpose()
            firstOrderMatrix[ind2, ind1] = firstOrderMatrix[ind2, ind1] + intermediateMatrix
        self.firstOrderMatrix = firstOrderMatrix
        return firstOrderMatrix
            
    def GetSecondOrderMatrix (self):
        '''
        This method returns element's local second order matrix
        '''
        circuitMatrix = self.GetCircuitMatrix()   
        numberOfEdges = circuitMatrix.shape[0]
        elementNodes = array(circuitMatrix[:,:2].astype(Int32))
        elementData = circuitMatrix[:,2:]              
        secondOrderMatrix = zeros((numberOfEdges+1, numberOfEdges+1))
        reverse = array([[1,-1],[-1,1]])
        for i in arange(0,numberOfEdges):
            compliance = elementData[i,0]
            intermediateMatrix = dot(compliance,reverse)  
            ind1 = array([elementNodes[i,:],elementNodes[i,:]])
            ind2 = ind1.transpose()
            secondOrderMatrix[ind2, ind1] = secondOrderMatrix[ind2, ind1] + intermediateMatrix
        self.secondOrderMatrix = secondOrderMatrix
        return secondOrderMatrix
        
    def GetNumberOfNodes(self):
        '''
        This method returns element's number of nodes
        '''
        numberOfNodes = len(self.NodeIds)
        return numberOfNodes
        
    def GetNumberOfDofs(self):
        '''
        This method returns element's number of dof (degrees of freedom)
        '''
        numberOfDofs = len(self.dof)
        return numberOfDofs

    def SetParameterInHistory(self,parameterName,parameterValue,currentIncrementNumber):
        '''
        This method sets Parameter History dictionary. Each value (for each parameters name) is associated to:
        0->Current Timestep
        1->Previous Timestep
        2->2 steps back
        '''
        self.MaxHistorySize = 3
        if parameterName not in self.ParameterInfo:
            self.ParameterInfo[parameterName] = []
        if self.LastIncrementNumber < currentIncrementNumber:
            self.ParameterInfo[parameterName].insert(0,parameterValue)
        if len(self.ParameterInfo[parameterName]) > self.MaxHistorySize:
            self.ParameterInfo[parameterName] = self.ParameterInfo[parameterName][0:self.MaxHistorySize]
                  
elementFactory = dict()

def RegisterElementType(elementType,elementClass):
    '''
    Registering the element's class according to its type
    '''
    elementFactory[elementType] = elementClass

def NewElement(elementType, elId, nodeIds, elementParameters, side=None, name=None):
    '''
    This method creates a new element according to its mesh type.
    '''
    return elementFactory[elementType](elId,nodeIds,elementParameters,side,name)
    
class WavePropagationElement(Element):
    '''
    Each Element is marked by n nodes and unique Id.
    Side: Arterial or Venous side.
    Name: Vessel Name.
    Each element (wave propagation element) is modeled with RCL electrical circuit.
    (1/2*C+Rleakage/2) + L + R + (1/2*C+Rleakage/2).
    C =  Capacitor, used to model the storage capacity of the vessel.
    L =  Inductance, used to represent the inertia dominated impedance in the central core.
    R =  Resistance, used to represent the resistance to flow through the vessel.
    LeakageR = Linear resistance, used to represent  the resistance to flow through small side-branches.
    This class provide the following methods:
    SetWallThickness: a method for setting wall thickness from radius fixed ratio.
    SetRadius: a method for setting radius.
    SetResistance: a method for setting non linear resistance.
    SetCompliance: a method for setting non linear compliance.
    SetLeakage: a method for setting leakage elements number.
    SetLinearValues: a method for setting linear values computed after a pre-run steady non linear simulation.
    Womersley: a method for building functions from Womersley Model used for calculating R1 and L.
    InputParameters: a method for calculating C, R and L from input parameters.
    GetCircuitMatrix: a method for building local circuit matrix.
    GetExternalPressureLocalDofs: a method for setting Transmural pressure in the correct local dofs.
    GetVenousPressureLocalDofs:
    GetPoiseuilleDofs: a method for getting Poiseuille's resistance local dofs.
    GetFlow: a method for calculating volumetric flow rate on the poiseuille resistance.
    GetWss: a method for calculating poiseuille wall shear stress on the poiseuille resistance.
    GetPressure: a method for calculating pressure over the element.
    GetArea: a method for calculation vessel's cross sectional area.
    GetLength: a method for returning element's length.
    GetRadius: a method for returning element's radius.
    GetRadiusAtRest: a method for returning element's radius at rest.
    GetRadius_a: a method for returning element's x axis radius.
    GetRadius_b: a method for returning element's y axis radius.
    GetWallThickness: a method for returning element's wall thickness.
    GetYoungModulus: a method for returning element's young's modulus.
    GetLocalDof and GetNodeLocalDofs: two methods for mapping element's NodeIds in local dof (if possible).
    GetDofNodes: a method for mapping local dof numbers in his NodeIds.
    '''

    def __init__(self, elId, nodeIds, elementParameters, side=None, name=None):
        '''
        Constructor
        '''
        Element.__init__(self)
        
        self.Type = "WavePropagation"
        self.Side = side
        self.Id = elId
        self.Name = name
        self.NodeIds = []
        self.NodeIds[:] = nodeIds[:]  
        self.nonLinearParameter = {}
        try:
            self.Resistance = elementParameters["resistance"]
        except KeyError:
            self.Resistance = None
        try:
            self.Compliance = elementParameters["compliance"]
        except KeyError:
            self.Compliance = None
        try:
            self.QLeakage = elementParameters["leakage"]
        except KeyError:
            self.QLeakage = None   
        
        self.s1 = elementParameters["s1"]
        self.s2 = elementParameters["s2"]
        self.Length = elementParameters["length"]
        self.Radius = elementParameters["radius"]
        self.dayRadius = {} #day:radius
        try:
            self.RadiusAtRest = elementParameters["radiusAtRest"]
            self.RadiusExp = self.Radius
            self.nonLinearParameter['radiusExp'] = True      
            self.Radius = self.RadiusAtRest         
        except KeyError:
            self.RadiusAtRest = None
            self.RadiusExp = None
        self.xRadius = elementParameters["xradius"]
        self.yRadius = elementParameters["yradius"]
        self.WallThickness = elementParameters["wall_thickness"]
        self.YoungModulus = elementParameters["young_modulus"]
        self.dz = self.Length/1.0e5
        for name in elementParameters:
            if name != 'wall_thickness' and name != 'leakage':
                if type(elementParameters[name]) is str:
                    self.nonLinearParameter[name] = True
                else:
                    self.nonLinearParameter[name] = False
        for val in self.nonLinearParameter.itervalues():
            if val == True:
                self.nonLinear = True
                break
        self.C = 0.0
        self.R = 0.0
        self.L = 0.0
        self.Leakages = None
        self.LeakageR = 0.0
        self.dof = [0,1,2,3,4]
        self.Flow = None
        self.Pressure = None
        self.Wss = None
        self.Initialized = False

    def SetWallThickness(self, wallthickness, info):
        '''
        This method sets WallThickness.
        '''
        self.WallThickness = wallthickness
        if info['history'] != []:
            self.LastIncrementNumber = info['history'][0]
            self.SetParameterInHistory('WallThickness', self.WallThickness, info['incrementNumber'])
        
    def SetRadius(self, radius, info=None, timeIndex=0):
        '''
        This method sets Radius.
        '''
        self.Radius = radius
        if info['history'] != []:
            self.LastIncrementNumber = info['history'][0]
            self.SetParameterInHistory('Radius', self.Radius, info['incrementNumber'])
    
    def SetResistance(self, resistance, info=None, timeIndex=0):
        '''
        This method sets non linear resistance.
        '''
        self.R = resistance
        if info['history'] != []:
            self.LastIncrementNumber = info['history'][0]
            self.SetParameterInHistory('Resistance', self.R, info['incrementNumber'])   
    
    def SetCompliance(self, compliance, info=None, timeIndex=0):
        '''
        This method sets non linear compliance.
        '''
        self.C2 = compliance
        self.C = compliance
        if info['history'] != []:
            self.LastIncrementNumber = info['history'][0]
            self.SetParameterInHistory('Compliance', self.C, info['incrementNumber'])
        
    def SetQLeakage(self, qleakage, info):
        '''
        This method sets Leakage Resistance.
        '''
        self.LeakageR = qleakage*self.Leakages
        
    def SetLinearValues(self,parameters):
        '''
        This method sets computed linear values after a pre-run non linear steady simulation.
        '''
        for el in parameters:
            if el == 'Compliance':
                self.Compliance = self.C/self.Length
                self.C = 0.0
                self.nonLinearParameter['compliance'] = False
            if el == 'Radius':
                self.R = 0.0
                self.L = 0.0
                self.nonLinearParameter['radiusExp'] = False
        self.Initialized = False
        self.Flow = None
        self.Pressure = None
        self.Wss = None
        self.nonLinear = False       
        for val in self.nonLinearParameter.itervalues():
            if val == True:
                self.nonLinear = True
                break
    
    def Womersley (self, r):
        '''
        Calculating falfa and galfa. Functions used for calculating RAlpha and LAlpha,
        parameters from the momentum equation. Cp and Cq are evaluated for the characteristic
        frequencyl.Alpha is the Womersley number corresponding to the characteristic frequency.
        '''
        self.nu=self.eta/self.rho               #cinematic viscosity    
        r = mean(r)
        self.alpha=r*sqrt(self.omega/self.nu)   #Womersley number    
        if (self.alpha <= sqrt(2)):         
            Cp = (3.0/2.0)
            Cq = (1.0/2.0)  
        else:    
            Cp = 1.0+(sqrt(2.0)/self.alpha)*(1.0-(sqrt(2.0)/(2.0*self.alpha)))
            Cq = (self.alpha/(4.0*sqrt(2.0)))*pow(1.0-(sqrt(2.0)/(2.0*self.alpha)),-1)     
        falpha = Cq/(2.0-Cp)
        galpha = 1.0/(2.0-Cp)
        return falpha, galpha
    
    def InputParameters(self, evaluator=None):
        '''
        This method calculates C, R and L from element's parameters:
        R and L are computed from the momentum equation (boundary layer theory)
        Veins can have elliptical cross-sectional area and therefore different expressions
        for L and R are needed for the venous segments.(Haslam et al. 1998)
        We assume a linear relation between pressure and area.
        Compliance (C) is computed from the vessel radius, the vessel wall thickness
        and the Young's modulus (Bessems et al. 2007) assuming that the artery is a 
        thick-walled linear elastic tube.
        Veins can have an elliptical cross-sectional area and therefore, instead of radius
        a weighted average radius ( sqrt(ao*bo) ) is used. Because veins are much thinner than
        arteries.
        If resistance and/or compliance are expressed with non-linear equations,
        non linear values will overwrite linear ones.
        All parameters are calculated integrating C, R and L over the segment length.
        '''
        Element.InputParameters(self)    
        try:
            self.eta = self.simulationContext.Context['dynamic_viscosity']
        except KeyError:
            print "Error, Please set Dynamic Viscosity[Pa*s] in Boundary Conditions XML File"
            raise
        try:
            self.rho = self.simulationContext.Context['blood_density']
        except KeyError:
            print "Error, Please set Blood Density[kg*m^3] in Boundary Conditions XML File" 
            raise  
        try:
            self.mu = self.simulationContext.Context['poisson_ratio']
        except KeyError:
            print "Error, Please set Poisson Ratio in Boundary Conditions XML File"
            raise
        try:
            self.freq = 1.0/(self.simulationContext.Context['period'])
            self.omega = 2*pi*self.freq    #phase
        except KeyError:
            print "Error, Please set Frequency[Hz] in Boundary Conditions XML File"
            raise
       
        #Computing dz steps and finding curvilinear abscissas
        z = arange(0.0,self.Length,self.dz)
        s1 = self.s1
        s2 = self.s2
        
        #Element is not initialized. Computing linear and non-linear parameters.
        if self.Initialized == False:
            #Radius
            if self.Radius is not None:
                if type(self.Radius) is dict:    
                    r1 = ((self.Radius[s2] - self.Radius[s1])/self.Length)
                    r2 = self.Radius[s1]
                    r_z = r2+(r1*z)
                    self.Radius = r_z
                    self.RadiusAtRest = self.Radius
               
            #xRadius and yRadius, for elliptical geometry
            if self.xRadius is not None:
                xr1 = ((self.xRadius[s2] - self.xRadius[s1])/self.Length)
                xr2 = self.xRadius[s1]
                xr_z = xr2+(xr1*z)
                self.xRadius = xr_z
            if self.yRadius is not None:
                yr1 = ((self.yRadius[s2] - self.yRadius[s1])/self.Length) 
                yr2 = self.yRadius[s1]
                yr_z = yr2+(yr1*z)
                self.yRadius = yr_z
                self.Radius = (self.xRadius*self.yRadius)**0.5
                
            #Wall-Thickness
            if type(self.WallThickness) is not str:
                h1 = (self.WallThickness[s2] - self.WallThickness[s1])/self.Length
                h2 = self.WallThickness[s1]
                h_z = h2+(h1*z)
                self.WallThickness = h_z
            else:
                evaluator.SetAbscissa(self.s1+((self.s2-self.s1)/2))
                evaluator.Evaluate(self.WallThickness)
                
            #Young's Modulus
            E1 = (self.YoungModulus[s2] - self.YoungModulus[s1])/self.Length
            E2 = self.YoungModulus[s1]
            E_z = E2+(E1*z)
            self.YoungModulus = E_z 
            
            #Resistance
            if self.Resistance == None:
                #Resistance without elliptical geometry
                if self.xRadius == None and self.yRadius == None:
                    R = (8.0*self.eta*self.dz)/(pi*self.Radius**4)
                    R = float(sum(R))
                else:
                    #Resistance with elliptical geometry
                    R = (8.0*self.eta*self.dz*((self.xRadius*self.xRadius)+(self.yRadius*self.yRadius)))/(2.0*pi*self.xRadius**3*self.yRadius**3)
                    R = float(sum(R))
            else:
                #Specific Linear Resistance
                if type(self.Resistance) is not str:
                    self.R = self.Resistance
                else:
                    #Specific Non-Linear Resistance
                    evaluator.Evaluate(self.Resistance)
                    R = self.R
                    R = float(mean(R))
                    
            Ralpha = float(R * self.Womersley(self.Radius)[0])
            
            #Inductance
            #Inductance without elliptical geometry
            if self.xRadius == None and self.yRadius == None:
                L = (self.rho*self.dz)/(pi*self.Radius**2)
                L = float(sum(L))
            else:
                #Inductance with elliptical geometry
                L = (self.rho*self.dz)/(pi*self.xRadius*self.yRadius)
                L = float(sum(L))
            
            Lalpha = float(L * self.Womersley(self.Radius)[1])
            
            #Compliance
            if self.Compliance == None:
                #Arterial Compliance computed from Young's Modulus and Wall-thickness 
                if self.Side == "arterial":
                    self.C = ((2.0*pi*self.Radius**2)*(((2.0*self.Radius**2*(1.0-self.mu**2))/(self.WallThickness**2))+((1.0+self.mu)*(((2.0*self.Radius)/self.WallThickness)+1.0)))*self.dz)/(self.YoungModulus*(((2.0*self.Radius)/self.WallThickness)+1.0))         
                    self.C = float(sum(self.C))
                if self.Side == "venous":
                    #Venous Compliance computed from Young's Modulus and Wall-thickness without elliptical geometry
                    if self.xRadius == None and self.yRadius == None:
                        self.C = ((2.0*pi*(sqrt(self.Radius*self.Radius))**3)*(1.0-self.mu**2)*self.dz)/(self.YoungModulus*self.WallThickness)
                    #Venous Compliance computed from Young's Modulus and Wall-thickness with elliptical geometry
                    else:
                        self.C = ((2.0*pi*(sqrt(self.xRadius*self.yRadius))**3)*(1.0-self.mu**2)*self.dz)/(self.YoungModulus*self.WallThickness)
                    self.C = float(sum(self.C))
            else:
                #Specific Linear Compliance
                if type(self.Compliance) is not str:
                    self.C = self.Compliance*self.Length
                else:
                    #Specific Non-Linear Compliance
                    evaluator.SetAbscissa(self.s1+((self.s2-self.s1)/2))
                    evaluator.Evaluate(self.Compliance)
                    self.C = self.C*self.dz
                    self.C = float(sum(self.C))
                    
            #Leakage Resistance
            if self.QLeakage is None:
                self.LeakageR = 1.0e25
            else:
                evaluator.Evaluate(self.QLeakage)
            
            #Element Initialized
            self.C = float(self.C)
            self.R = float(Ralpha)    
            self.L = float(Lalpha)
            self.Initialized = True   
            
        if self.Initialized == True:
            
            #Radius
            for name, value in self.nonLinearParameter.iteritems():
                if value == True and name == 'radiusExp':
                    evaluator.SetAbscissa(self.s1+((self.s2-self.s1)/2))
                    evaluator.Evaluate(self.RadiusExp[s1])
                    #Resistance and Inductance are radius dependents
                    R = (8.0*self.eta*self.dz)/(pi*self.Radius**4)
                    R = float(sum(R))
                    L = (self.rho*self.dz)/(pi*self.Radius**2)
                    L = float(sum(L))
                    Ralpha = float(R * self.Womersley(self.Radius)[0])
                    Lalpha = float(L * self.Womersley(self.Radius)[1])
                    self.L = Lalpha
                    self.R = Ralpha 
                    
            #xRadius and yRadius for elliptical geometry
            for name, value in self.nonLinearParameter.iteritems():
                if value == True and name == 'xradius':
                    pass
            for name, value in self.nonLinearParameter.iteritems():
                if value == True and name == 'yradius':
                    pass
            
            #Wall-Thickness
            for name, value in self.nonLinearParameter.iteritems():
                if value == True and name == 'wall_thickness':
                    evaluator.SetAbscissa(self.s1+((self.s2-self.s1)/2))
                    evaluator.Evaluate(self.WallThickness)              
            
            #Resistance        
            for name, value in self.nonLinearParameter.iteritems():              
                if value == True and name == 'resistance':         
                    evaluator.Evaluate(self.Resistance)
                    R = self.R
                    R = float(mean(R))
                    Ralpha = float(R * self.Womersley(self.Radius)[0])
                    self.R = float(Ralpha)
                    
            #Compliance
            for name, value in self.nonLinearParameter.iteritems():
                if value == True and name == 'compliance':
                    evaluator.SetAbscissa(self.s1+((self.s2-self.s1)/2))
                    evaluator.Evaluate(self.Compliance)
                    self.C = self.C *self.dz
                    self.C = float(sum(self.C))

        return self.C, self.R, Ralpha, Lalpha, self.LeakageR
        
    def GetCircuitMatrix(self):
        '''
        This method builds element's circuit matrix
        Each Row is an edge, Node1 - Node2 - C - R - L
        '''
        
        CircuitMatrix = array ([[self.dof[0], self.dof[1], 0, 0, self.L],                        # Inductance
                               [self.dof[1], self.dof[2], 0, self.R, 0],                         # Resistance R
                               [self.dof[0], self.dof[3], 0.5*self.C, self.LeakageR, 0],         # Left C/2 // Leakage
                               [self.dof[2], self.dof[4], 0.5*self.C, self.LeakageR, 0]])        # Right C/2 // Leakage
        
        return CircuitMatrix
                
    def GetExternalPressureLocalDofs(self):
        '''
        Setting Transmural pressure in the correct local dofs.
        '''
        return [self.dof[3], self.dof[4]]
    
    def GetVenousPressureLocalDofs(self):
        '''
        Setting Nodal Output Prescribed Pressure
        '''
        return self.dof[2]
        
    def GetPoiseuilleDofs(self):
        '''
        This method return Poiseuille's resistance local dofs
        '''
        return [self.dof[1], self.dof[2]]
    
    def GetFlow(self, info, timeIndex=0):
        '''
        This method returns volumetric flow rate calculated on the poiseuille resistance.(mL/min)
        If cycle is not specified, default cycle is the last one.
        '''
        # t=0, no flow.
        if info['solution'] is None:
            self.Flow = 1.0e-25
            return self.Flow
        try:
            self.Period = self.simulationContext.Context['period']
        except KeyError:
            print "Error, Please set period in Boundary Conditions XML File"
            raise
        try:
            self.Cycles = self.simulationContext.Context['cycles']
        except KeyError:
            print "Error, Please set cycles number in Boundary Conditions XML File"
            raise
        try:
            self.TimeStep = self.simulationContext.Context['timestep']
        except KeyError:
            print "Error, Please set timestep in Boundary Conditions XML File"
            raise
        try:
            solution = info['solution'][timeIndex]
        except KeyError:
            print "Error, Please provide Solution"
            raise
        try:
            dofmap = info['dofmap']
        except KeyError:
            print "Error, Please provide Dofmap"
            raise
        try:
            Cycle = info['cycle']
        except KeyError:
            Cycle = self.Cycles
        
        dofs = self.GetPoiseuilleDofs()
        self.Flow = (solution[(dofmap.DofMap[self.Id, dofs[0]]),:] - solution[(dofmap.DofMap[self.Id, dofs[1]]),:])/self.R
        if len(self.Flow) != 1:
            self.Flow = mean(self.Flow[(int(self.Period/self.TimeStep)*(Cycle-1)):(int(self.Period/self.TimeStep)*(Cycle))])*6.0e7
        else:
            self.Flow = self.Flow[0]*6.0e7 
        return self.Flow
    
    def GetWss(self, info, timeIndex=0):
        '''
        This method returns Wall Shear Stress on the specified element.(Pa)
        Wall Shear Stress is computed on the Poiseuille Resistance.
        If element is tapered, Radius is considered as mean value over segment length.
        '''
        self.Wss = ((4.0*self.eta)/6.0e7*pi) * (self.GetFlow(info,timeIndex)/(mean(self.Radius)**3))
        return self.Wss
    
    def GetPressure (self, info, timeIndex=0):
        '''
        This method returns pressure on the specified element's first node.
        If cycle is not specified, default cycle is the last one.
        '''
        # t=0, no pressure.
        if info['solution'] is None:
            self.Pressure = 1e-12
            return self.Pressure 
        try:
            self.Period = self.simulationContext.Context['period']
        except KeyError:
            print "Error, Please set period in Boundary Conditions XML File"
            raise
        try:
            self.Cycles = self.simulationContext.Context['cycles']
        except KeyError:
            print "Error, Please set cycles number in Boundary Conditions XML File"
            raise
        try:
            self.TimeStep = self.simulationContext.Context['timestep']
        except KeyError:
            print "Error, Please set timestep in Boundary Conditions XML File"
            raise  
        try:
            solution = info['solution'][timeIndex]
        except KeyError:
            print "Error, Please provide Solution"
            raise
        try:
            dofmap = info['dofmap']
        except KeyError:
            print "Error, Please provide Dofmap"
            raise
        try:
            Cycle = info['cycle']
        except KeyError:
            Cycle = self.Cycles
        
        dofs = self.GetPoiseuilleDofs()
        self.Pressure = (solution[(dofmap.DofMap[self.Id, dofs[0]]),:]) 
        if len(self.Pressure) != 1:
            self.Pressure = mean(self.Pressure[(int(self.Period/self.TimeStep)*(Cycle-1)):(int(self.Period/self.TimeStep)*(Cycle))])
        else:
            self.Pressure = self.Pressure[0] 
        if self.Pressure <= 0.0:
            self.Pressure = 1e-12
        return self.Pressure

    def GetArea(self,info, timeIndex=0):
        '''
        This method returns vessel's Cross-Sectional Area
        '''   
        if self.Side == 'arterial':
            if type(self.Radius) == dict:
                Radius1 = self.Radius[self.s1]           
                Radius2 = self.Radius[self.s2]
                Areas1 = pi*(Radius1**2)
                Areas2 = pi*(Radius2**2)       
            else:
                Radius = self.Radius[len(self.Radius)-1]
                Radius1 = Radius2 = Radius
                Areas1 = pi*(Radius**2)
                Areas2 = pi*(Radius**2)
        if self.Side == 'venous':
            if type(self.Radius) == dict:
                Radius1 = self.Radius[self.s1]
                Radius2 = self.Radius[self.s2]
                Areas1 = pi*(Radius1**2)
                Areas2 = pi*(Radius2**2)
            else:
                Radius = self.Radius[0]
                Areas1 = pi*(Radius**2)
                Areas2 = pi*(Radius**2)          
        return Areas1, Areas2
    
    def GetLength(self,info):
        '''
        This method returns Length
        '''
        Length = self.Length
        return Length
    
    def GetRadius(self, info, timeIndex=0):
        '''
        This method returns Radius
        '''
        try:
            Radius = self.ParameterInfo['Radius'][timeIndex] 
        except:
            Radius = self.Radius
        return Radius
    
    def GetRadiusAtRest(self, info, timeIndex=0):
        '''
        This method returns Radius
        '''
        RadiusAtRest = self.RadiusAtRest
        return RadiusAtRest
    
    def GetRadius_a(self, info, timeIndex=0):
        '''
        This method returns Radius_a axis
        '''
        Radius_a = self.xRadius
        return Radius_a
    
    def GetRadius_b(self, info, timeIndex=0):
        '''
        This method returns Radius_b axis
        '''
        Radius_b = self.yRadius
        return Radius_b
    
    def GetWallThickness(self,info):
        '''
        This method returns WallThickness
        '''
        WallThickness = self.WallThickness
        return WallThickness
    
    def GetYoungModulus(self,info):
        '''
        This method returns YoungModulus
        '''
        YoungModulus = self.YoungModulus
        return YoungModulus
    
    def GetLocalDof (self, NodeId):
        '''
        This method returns Local dof number corresponding to specific NodeId 
        '''
        if NodeId == self.NodeIds[0]:
            LocalDof = 0
        if NodeId == self.NodeIds[1]:
            LocalDof = 2
        return LocalDof
    
    def GetNodeLocalDofs(self):
        '''
        This method returns local dof number corresponding to its NodeId (if exist)
        '''
        for dofs in self.NodeIds:
            if dofs == self.NodeIds[0]:
                NodeDof1 = 0
            if dofs == self.NodeIds[1]:
                NodeDof2 = 2 
        NodeDofs = [NodeDof1, NodeDof2]
        return NodeDofs
        
    def GetDofNodes(self):
        '''
        This method returns NodeId corresponding to local dof number (if exist)
        '''
        for dofs in self.dof:
            if dofs == 0:
                DofNodeId1 = self.NodeIds[0]
            if dofs == 2:
                DofNodeId2 = self.NodeIds[1]
        DofNodes = [DofNodeId1,DofNodeId2]
        return DofNodes

RegisterElementType(None,WavePropagationElement)
RegisterElementType("WavePropagation",WavePropagationElement)

class WindkesselElement(Element):
    '''
    Windkessel Element is an element used for downstream network.
    Windkessel Element is marked by 2 nodes and unique Id.
    Windkessel Element has 3 local dofs.
    Side: Arterial or Venous side (optional)
    Name: End Vessel Name.
    Each Windkessel Element is modeled like a 0D Lumped parameter Windkessel Element consisting of
    a resistance R1 in series with a parallel resistance R2 and a Compliance C.
    R1: Wave impedance. The value is chosen such that wave reflections at the interface between the Windkessel element
    and the wave propagation element are minimal.
    R2 : Peripheral Resistance.
    C : Compliance.
    R2 and C are determined by fitting flow as computed locally in the model with that measured experimentally.
    This class provides the following methods:
    SetLastElement: a method for setting the last element connected to the current end segment.
    SetWindkesselRel: a method for setting the peripheral resistance of the current end segment.
    SetLinearValues: a method for setting linear values computed after a pre-run steady non linear simulation.
    InputParameters: This methods calculates windkessel parameters for specific patient assuming general flow distributions.   
    GetCircuitMatrix: a method for building local circuit matrix.
    GetExternalPressureLocalDofs: a method for setting Transmural pressure in the correct local dofs.
    GetDofNodes: a method for mapping local dof numbers in his NodeIds.
    GetLocalDof: a method for mapping element's NodeIds in local dof (if possible).
    '''
    def __init__(self, elId, NodeIds, name, side=None):
        '''
        Constructor
        '''
        Element.__init__(self)
        self.Type = "Windkessel"
        self.Side = side
        self.Id = "E%s" % (elId)
        self.Name = name
        self.NodeIds = []
        self.NodeIds[:] = NodeIds[:]
        self.LastElement = None
        self.Leakage = False
        self.LeakageR = 0.0
        self.Rel = None
        self.RelExpression = None
        self.R1 = 0.0
        self.C = 0.0
        self.R2 = 0.0
        self.dof = [0, 1, 2, 3]
        self.Initialized = False

    def SetLastElement (self, element):
        '''
        This method sets the element to which is connected.
        '''
        self.LastElement = element
    
    def SetWindkesselRel(self, rel, info):
        '''
        This method sets windkessel peripheral resistance which
        is used for generating windkessel element.
        '''
        self.Rel = rel
        
    def SetLinearValues(self,parameters):
        '''
        This method is used for setting linear values computed after a pre-run steady non linear simulation.
        '''
        self.Initialized = False
        self.R1 = 0.0
        self.C = 0.0
        self.R2 = 0.0
        
    def InputParameters(self, evaluator=None):
        '''
        This methods calculates windkessel parameters for specific patient
        assuming general flow distributions
        ''' 
        Element.InputParameters(self)
        evaluator.Evaluate(self.RelExpression)
        try:
            self.R1 = sqrt(self.LastElement.L/self.LastElement.C)
        except:
            self.R1 = 1.0
        self.R2 = self.Rel-self.R1
        self.C = 1.1 / self.R2
        self.Initialized = True
        
    def GetCircuitMatrix(self):
        '''
        This method builds element's circuit matrix
        Each Row is an edge, Node1 - Node2 - C - R - L
        '''
        circuitMatrix = array ([[self.dof[0], self.dof[1], 0, self.R1, 0],
                                [self.dof[1], self.dof[2], 0, self.R2, 0],
                                [self.dof[1], self.dof[3], self.C, 0 , 0]])
        return circuitMatrix
      
    def GetLocalDof (self, NodeId):
        '''
        This method returns Local dof number corresponding to specific NodeId 
        '''
        if NodeId == self.NodeIds[0]:
            LocalDof = 0
        if NodeId == self.NodeIds[1]:
            LocalDof = 2    
        return LocalDof
    
    def GetExternalPressureLocalDofs(self):
        '''
        Setting Transmural pressure in the correct local dofs.
        '''
        return [self.dof[2], self.dof[3]]
        
    def GetDofNodes(self):
        '''
        This method returns NodeId corresponding to local dof number (if exists)
        '''
        for dofs in self.dof:
            if dofs == 0:
                dofNodeId1 = self.NodeIds[0]
            if dofs == 2:
                dofNodeId2 = self.NodeIds[1]
            
        dofNodes = [dofNodeId1,dofNodeId2]
        return dofNodes

RegisterElementType("Windkessel",WindkesselElement)


class Anastomosis(Element):
    '''
    Anastomosis Element is a 3-nodes Element marked by its NodeIds and its unique Id.
    This Element is composed by 2 resistances, one located between the first and the second Node
    and the second one located between the first and the third node.
    Side: Arterial or Venous side (optional)
    Name: Anastomosis Name (optional)
    Each Resistance is modeled like a resistance (linear or not) between two nodes.(R_0_1 and R_0_2)
    This class provides the following methods:
    SetResistance_0_1: a method for setting the first resistance between dof 0 and dof 1.
    SetResistance_0_2: a method for setting the second resistance between dof 0 and dof 2.
    SetProximal: a method for setting connections, proximal artery.
    SetDistal: a method for setting connections, distal artery.
    SetVein: a method for setting connections, vein.
    SetLinearValues: a method for setting linear values computed after a pre-run steady non linear simulation.
    InputParameters: a method for calculating R_0_1 and R_0_2 from element's parameters.
    GetFlowProximal: a method for computing volumetric flow rate calculated on the poiseuille resistance(mL/min) of the proximal artery.
    GetFlowVein: a method for computing volumetric flow rate calculated on the poiseuille resistance(mL/min) of the vein.
    GetFlowDistal: a method for computing volumetric flow rate calculated on the poiseuille resistance(mL/min) of the distal artery.
    GetRadiusProximal: a method for returning proximal artery radius.
    GetRadiusDistal: a method for returning distal artery radius.
    GetRadiusVein: a method for returning vein radius.
    GetAreaProximal: a method for computing proximal artery cross sectional area.
    GetAreaDistal: a method for computing distal artery cross sectional area.
    GetAreaVein: a method for computing vein cross sectional area.
    GetFlowRatio: a method for computing the ratio between vein and proximal artery volumetric flow rates.
    GetRadiusRatio: a method for computing the ratio between vein and proximal artery radii.
    GetAreaRatio: a method for computing the ratio between vein and proximal artery cross sectional areas.
    GetCircuitMatrix: a method for building local circuit matrix.
    GetExternalPressureLocalDofs: a method for setting Transmural pressure in the correct local dofs.
    GetDofNodes: a method for mapping local dof numbers in his NodeIds.
    GetLocalDof: a method for mapping element's NodeIds in local dof (if possible).
    '''
    def __init__(self, elId, nodeIds, elementParameters, side=None, name=None):
        '''
        Constructor
        '''
        Element.__init__(self)
        self.Type = "Anastomosis"
        self.Side = side
        self.Id = elId
        self.Name = name
        self.NodeIds = []
        self.NodeIds[:] = nodeIds[:]
        self.Proximal = None
        self.Distal = None
        self.Vein = None

        try:
            self.Resistance_0_1 = elementParameters["resistance_0_1"]
        except KeyError:
            self.Resistance_0_1 = None      
        try:
            self.Resistance_0_2 = elementParameters["resistance_0_2"]
        except KeyError:
            self.Resistance_0_2 = None
        
        for name in elementParameters:
            if type(elementParameters[name]) is str:
                self.nonLinear = True
        self.R_0_1 = None
        self.R_0_2 = None
        self.dof = [0,1,2]
        self.Flow = None    
        self.Initialized = False
        
    def SetResistance_0_1(self, resistance, info):
        '''
        This method sets Resistance_0_1.
        '''
        self.R_0_1 = resistance
        if info['history'] != []:
            self.LastIncrementNumber = info['history'][0]
            self.SetParameterInHistory('Resistance_0_1', self.R_0_1, info['incrementNumber'])
        
    def SetResistance_0_2(self, resistance, info):
        '''
        This method sets Resistance_0_2.
        '''
        self.R_0_2 = resistance
        if info['history'] != []:
            self.LastIncrementNumber = info['history'][0]
            self.SetParameterInHistory('Resistance_0_2', self.R_0_2, info['incrementNumber'])
            
    def SetProximal(self,proximal):
        '''
        Setting connections, proximal artery.
        '''
        self.Proximal = proximal
    
    def SetDistal(self,distal):
        '''
        Setting connections, distal artery.
        '''
        self.Distal = distal
    
    def SetVein(self,vein):
        '''
        Setting connections, proximal vein.
        '''
        self.Vein = vein
        
    def SetLinearValues(self,parameters):
        '''
        This method sets linear values computed after a pre-run steady non linear simulation.
        '''
        self.R_0_1 = None
        self.R_0_2 = None
        self.Flow = None    
        self.Initialized = False
          
    def InputParameters(self, evaluator=None):
        '''
        This method calculates R_0_1 and R_0_2 from element's parameters:
        If resistance is expressed with non-linear equations,
        non linear value will overwrite the linear one.
        '''
        Element.InputParameters(self)
        if type(self.Resistance_0_1) is  str:
            evaluator.Evaluate(self.Resistance_0_1)
        else:    
            self.R_0_1 = self.Resistance_0_1
        if type(self.Resistance_0_2) is  str:
            evaluator.Evaluate(self.Resistance_0_2)
        else:    
            self.R_0_2 = self.Resistance_0_2
        self.Initialized = True
        
    def GetFlowProximal(self, info, timeIndex=0):
        '''
        This method returns volumetric flow rate calculated on the poiseuille resistance.(mL/min)
        '''
        self.Flow = self.Proximal.GetFlow(info, timeIndex)
        return self.Flow
    
    def GetFlowVein(self, info, timeIndex=0):
        '''
        This method returns volumetric flow rate calculated on the poiseuille resistance.(mL/min)
        '''           
        self.Flow = self.Vein.GetFlow(info, timeIndex)
        return self.Flow
    
    def GetFlowDistal(self, info, timeIndex=0):
        '''
        This method returns volumetric flow rate calculated on the poiseuille resistance.(mL/min)
        '''
        self.Flow = self.Distal.GetFlow(info, timeIndex)
        return self.Flow
    
    def GetRadiusProximal(self, info, timeIndex=0):
        '''
        This method returns vessel's Cross-Sectional Radius
        '''
        Radius = self.Proximal.GetRadius(info, timeIndex=0)
        RadiusProx = Radius[len(Radius)-1]
        return RadiusProx
    
    def GetRadiusVein(self, info, timeIndex=0):
        '''
        This method returns vessel's Cross-Sectional Radius
        '''
        if type(self.Vein.GetRadius(info, timeIndex=0)) is dict:
            Radius = self.Vein.GetRadius(info, timeIndex=0)[0.0]
        else:
            Radius = self.Vein.GetRadius(info, timeIndex=0)[0]
        return Radius
    
    def GetAreaProximal(self, info, timeIndex=0):
        '''
        This method returns vessel's Cross-Sectional Area
        '''
        Area = self.Proximal.GetArea(info, timeIndex)[1]
        return Area
    
    def GetAreaDistal(self, info, timeIndex=0):
        '''
        This method returns vessel's Cross-Sectional Area
        '''
        Area = self.Distal.GetArea(info, timeIndex)[0]
        return Area
    
    def GetAreaVein(self, info, timeIndex=0):
        '''
        This method returns vessel's Cross-Sectional Area
        '''
        Area = self.Vein.GetArea(info, timeIndex)[0]
        return Area
    
    def GetFlowRatio(self, info, timeIndex=0):
        '''
        This method returns the ratio between the volumetric flow rate
        through the vein and the volumetric flow rate through the proximal artery.(mL/min)
        '''
        self.FlowRatio = self.GetFlowVein(info, timeIndex)/self.GetFlowProximal(info, timeIndex)
        return self.FlowRatio
    
    def GetRadiusRatio(self, info, timeIndex=0):
        '''
        This method returns the ratio between the cross-sectional Radius of the vein and
        the cross-sectional Radius of the proximal artery (m)
        '''
        self.RadiusRatio = self.GetRadiusVein(info, timeIndex)/self.GetRadiusProximal(info, timeIndex)
        return self.RadiusRatio
        
    def GetAreaRatio(self, info, timeIndex=0):
        '''
        This method returns the ratio between the cross-sectional Area of the vein and
        the cross-sectional Area of the proximal artery (m2)
        '''
        self.AreaRatio = self.GetAreaVein(info, timeIndex)/self.GetAreaProximal(info, timeIndex)
        return self.AreaRatio
    
    def GetCircuitMatrix(self):
        '''
        This method builds element's circuit matrix
        Each Row is an edge, Node1 - Node2 - C - R - L
        '''        
        CircuitMatrix = array ([[self.dof[0], self.dof[1], 0, self.R_0_1, 0],         #Resistance_0_1
                                [self.dof[0], self.dof[2], 0, self.R_0_2, 0]])        #Resistance_0_2  
        return CircuitMatrix
        
    def GetLocalDof (self, NodeId):
        '''
        This method returns Local dof number corresponding to specific NodeId 
        '''
        if NodeId == self.NodeIds[0]:
            LocalDof = 0
        if NodeId == self.NodeIds[1]:
            LocalDof = 1
        if NodeId == self.NodeIds[2]:
            LocalDof = 2
        return LocalDof
    
    def GetExternalPressureLocalDofs(self):
        '''
        Setting Transmural pressure in the correct local dofs.
        '''
        return []
    
    def GetNodeLocalDofs(self):
        '''
        This method returns local dof number corresponding to its NodeId (if exist)
        '''
        for dofs in self.NodeIds:
            if dofs == self.NodeIds[0]:
                NodeDof1 = 0
            if dofs == self.NodeIds[1]:
                NodeDof2 = 1
            if dofs == self.NodeIds[2]:
                NodeDof3 = 2 
        NodeDofs = [NodeDof1, NodeDof2, NodeDof3]
        return NodeDofs
        
    def GetDofNodes(self):
        '''
        This method returns NodeId corresponding to local dof number (if exist)
        '''
        for dofs in self.dof:
            if dofs == 0:
                DofNodeId1 = self.NodeIds[0]
            if dofs == 1:
                DofNodeId2 = self.NodeIds[1]
            if dofs == 2:
                DofNodeId3 = self.NodeIds[2]
        DofNodes = [DofNodeId1, DofNodeId2, DofNodeId3]
        return DofNodes
    
RegisterElementType("Anastomosis",Anastomosis)


class ResistanceElement(Element):
    '''
    Resistance Element is an element realized with a non-linear resistance between two nodes.
    Resistance is marked by 2 nodes and unique Id.
    Resistance has 2 local dofs.
    Side: Arterial or Venous side (optional)
    Name: Resistance Name (optional)
    Each Resistance Element is modeled like a resistance (linear or not) between two nodes.(R)
    This class provides the following methods:
    SetResistance: a method for setting resistance.
    InputParameters: Setting Resistance value.
    GetPoiseuilleDofs: a method for getting poiseuille's resistance local dofs.
    GetCircuitMatrix: a method for building local circuit matrix.
    GetExternalPressureLocalDofs: a method for setting Transmural pressure in the correct local dofs.
    GetFlow : a method for calculating volumetric flow rate through the element
    GetDofNodes: a method for mapping local dof numbers in his NodeIds.
    GetLocalDof: a method for mapping element's NodeIds in local dof (if possible).
    '''
    def __init__(self, elId, nodeIds, elementParameters, side=None, name=None):
        '''
        Constructor
        '''
        Element.__init__(self)
        self.Type = "Resistance"
        self.Side = side
        self.Id = elId
        self.Name = name
        self.NodeIds = []
        self.NodeIds[:] = nodeIds[:]
        try:
            self.Resistance = elementParameters["resistance"]
        except KeyError:
            self.Resistance = None
        for name in elementParameters:
            if type(elementParameters[name]) is str:
                self.nonLinear = True
        self.R = None
        self.dof = [0,1]
        self.Flow = None
        self.Initialized = False
        
    def SetResistance(self, resistance, info=None, timeIndex=0):
        '''
        This method sets non linear resistance.
        '''
        self.R = resistance
        if info['history'] != []:
            self.LastIncrementNumber = info['history'][0]
            self.SetParameterInHistory('Resistance', self.R, info['incrementNumber'])   

    def InputParameters(self, evaluator=None):
        '''
        This method calculates R from element's parameters:
        If resistance is expressed with non-linear equations,
        non linear value will overwrite the linear one.
        '''
        Element.InputParameters(self)
        if type(self.Resistance) is  str:
            evaluator.Evaluate(self.Resistance)
        else:    
            self.R = self.Resistance
        self.Initialized = True
    
    def GetPoiseuilleDofs(self):
        '''
        This method return Poiseuille's resistance local dofs
        '''
        return [self.dof[0], self.dof[1]]
    
    def GetPressureGradient (self, info, timeIndex=0):
        '''
        This method returns pressure gradient on the resistance.
        If cycle is not specified, default cycle is the last one.
        '''
        # t=0, no pressure.
        if info['solution'] is None:
            self.Pressure = 1e-12
            return self.Pressure 
        try:
            self.Period = self.simulationContext.Context['period']
        except KeyError:
            print "Error, Please set period in Boundary Conditions XML File"
            raise
        try:
            self.Cycles = self.simulationContext.Context['cycles']
        except KeyError:
            print "Error, Please set cycles number in Boundary Conditions XML File"
            raise
        try:
            self.TimeStep = self.simulationContext.Context['timestep']
        except KeyError:
            print "Error, Please set timestep in Boundary Conditions XML File"
            raise  
        try:
            solution = info['solution'][timeIndex]
        except KeyError:
            print "Error, Please provide Solution"
            raise
        try:
            dofmap = info['dofmap']
        except KeyError:
            print "Error, Please provide Dofmap"
            raise
        try:
            Cycle = info['cycle']
        except KeyError:
            Cycle = self.Cycles
        
        dofs = self.GetPoiseuilleDofs()
        self.Pressure1 = (solution[(dofmap.DofMap[self.Id, dofs[0]]),:])
        self.Pressure2 = (solution[(dofmap.DofMap[self.Id, dofs[1]]),:])
        
        if len(self.Pressure1) != 1:
            self.Pressure1 = mean(self.Pressure1[(int(self.Period/self.TimeStep)*(Cycle-1)):(int(self.Period/self.TimeStep)*(Cycle))])
        else:
            self.Pressure1 = self.Pressure1[0] 
        if self.Pressure1 <= 0.0:
            self.Pressure1 = 1e-12
        if len(self.Pressure2) != 1:
            self.Pressure2 = mean(self.Pressure2[(int(self.Period/self.TimeStep)*(Cycle-1)):(int(self.Period/self.TimeStep)*(Cycle))])
        else:
            self.Pressure2 = self.Pressure2[0] 
        if self.Pressure2 <= 0.0:
            self.Pressure2 = 1e-12
            
        self.PressureGradient = self.Pressure2-self.Pressure1
        print self.PressureGradient
        return self.PressureGradient
    
    def GetCircuitMatrix(self):
        '''
        This method builds element's circuit matrix
        Each Row is an edge, Node1 - Node2 - C - R - L
        '''        
        CircuitMatrix = array ([[self.dof[0], self.dof[1], 0, self.R, 0]])        #Resistance
        return CircuitMatrix      
    
    def GetExternalPressureLocalDofs(self):
        '''
        Setting Transmural pressure in the correct local dofs.
        '''
        return []
    
    def GetFlow(self, info, timeIndex=0):
        '''
        This method returns volumetric flow rate calculated on the resistance.
        '''
        # t=0, no flow.
        if info['solution'] is None:
            self.Flow = 1.0e-25
            return self.Flow   
        try:
            solution = info['solution'][timeIndex]
        except KeyError:
            print "Error, Please provide Solution"
            raise
        try:
            dofmap = info['dofmap']
        except KeyError:
            print "Error, Please provide Dofmap"
            raise     
        self.Flow = (solution[(dofmap.DofMap[self.Id, self.dof[0]]),:] - solution[(dofmap.DofMap[self.Id, self.dof[1]]),:])/self.R
        return self.Flow
    
    def GetLocalDof (self, NodeId):
        '''
        This method returns Local dof number corresponding to specific NodeId 
        '''
        if NodeId == self.NodeIds[0]:
            LocalDof = 0
        if NodeId == self.NodeIds[1]:
            LocalDof = 1   
        return LocalDof
    
    def GetNodeLocalDofs(self):
        '''
        This method returns local dof number corresponding to its NodeId (if exist)
        '''
        for dofs in self.NodeIds:
            if dofs == self.NodeIds[0]:
                NodeDof1 = 0
            if dofs == self.NodeIds[1]:
                NodeDof2 = 1
        NodeDofs = [NodeDof1, NodeDof2]
        return NodeDofs
        
    def GetDofNodes(self):
        '''
        This method returns NodeId corresponding to local dof number (if exist)
        '''
        for dofs in self.dof:
            if dofs == 0:
                DofNodeId1 = self.NodeIds[0]
            if dofs == 1:
                DofNodeId2 = self.NodeIds[1]
        DofNodes = [DofNodeId1,DofNodeId2]
        return DofNodes
    
RegisterElementType("Resistance",ResistanceElement)

     
class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass