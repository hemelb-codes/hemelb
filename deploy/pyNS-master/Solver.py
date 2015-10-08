#!/usr/bin/env python

## Program:   PyNS
## Module:    Solver.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from Assembler import Assembler
from numpy.core.numeric import zeros, arange
from numpy.lib.function_base import delete
from numpy.lib.index_tricks import s_
from numpy.linalg.linalg import solve
from numpy.linalg import norm
from numpy.core.numeric import Inf, dot
from numpy.ma.core import ceil
import sys

class Solver(object):
    '''
    This is a general Solver Class. It doesn't provide any solver methods.
    This class provides only parameters setting methods
    '''
   
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkMesh = None
        self.SimulationContext = None
        self.BoundaryConditions = None
        self.Evaluator = None
        self.Solutions = None
        self.TimeStep = None
        self.SquareTimeStep = None
        self.Period = None
        self.CardiacFreq = None
        self.Cycles = None
        self.SimulationDays = []                            # Days' list for adaptation 
        self.NumberOfIncrements = None                     
        self.IncrementNumber = 1                            # increment number
        self.EndIncrementTime = 0.0                         # end increment time
        self.nltol = float(1e-3)                            # nonLinear convergence criterium
        self.convergence = float(1e-4)                      # steady convergence limit for steady pre runs    
        self.Flow = None
        self.PrescribedPressures = None
        self.LinearZeroOrderGlobalMatrix = None
        self.LinearFirstOrderGlobalMatrix = None
        self.LinearSecondOrderGlobalMatrix = None
        
    def SetNetworkMesh(self,networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetEvaluator(self, evaluator):
        '''
        Setting Evaluator
        '''
        self.Evaluator = evaluator
        
    def SetSimulationContext(self, simulationContext):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = simulationContext
        
    def SetBoundaryConditions(self, boundaryConditions):
        '''
        Setting BoundaryConditions
        '''
        self.BoundaryConditions = boundaryConditions
        
    def SetSteadyFlow(self):
        '''
        Setting Steady Flow
        '''
        self.steady = True
        
    def SetPulseFlow(self):
        '''
        Setting Pulsatile Flow
        '''
        self.steady = False
    
    def SetNonLinearTolerance(self, nltol):
        '''
        Setting Non Linear Tolerance Value
        '''
        self.nltol = float(nltol)
       
    def SetSteadyConvergenceLimit(self, convergence):
        '''
        Setting convergence limit for steady pre-run simulations
        '''
        self.convergence = convergence
 
class SolverFirstTrapezoid(Solver):
    '''
    This class provide a method to solving the system with "First Order Trapezium Method"
    '''

    def __init__(self):
        '''
        Constructor
        '''
        Solver.__init__(self)
        
    def Solve(self):
        '''
        This method builds System Matrix and gets Solution
        '''
        if self.SimulationContext.Id != self.NetworkMesh.Id:
            raise self.SimulationContext.XMLIdError()        
        try:
            self.TimeStep = self.SimulationContext.Context['timestep']
            self.SquareTimeStep = self.TimeStep*self.TimeStep
        except KeyError:
            print "Error, Please set timestep in Simulation Context XML File"
            raise
        try:
            self.Period = self.SimulationContext.Context['period']
            self.TimeStepFreq = int(self.Period/self.TimeStep)
        except KeyError:
            print "Error, Please set period in Simulation Context XML File"
            raise
        try:
            self.Cycles = self.SimulationContext.Context['cycles']
            self.NumberOfIncrements = (self.Cycles*self.TimeStepFreq)
        except KeyError:
            print "Error, Please set cycles number in Simulation Context XML File"
            raise
        
        history = []
        assembler = Assembler()
        assembler.SetNetworkMesh(self.NetworkMesh)
        assembler.SetBoundaryConditions(self.BoundaryConditions)
        info = {'dofmap':assembler.DofMap,'solution':None,'incrementNumber':self.IncrementNumber,'history':history}
        self.Evaluator.SetInfo(info)
        
        self.PrescribedPressures = assembler.AssembleBoundaryConditions(self.SimulationContext)
        
        self.LinearZeroOrderGlobalMatrix, self.LinearFirstOrderGlobalMatrix, self.LinearSecondOrderGlobalMatrix = \
        assembler.AssembleInit(self.SimulationContext, self.Evaluator)
        
        self.ZeroOrderGlobalMatrix = assembler.ZeroOrderGlobalMatrix
        self.FirstOrderGlobalMatrix = assembler.FirstOrderGlobalMatrix
        self.SecondOrderGlobalMatrix = assembler.SecondOrderGlobalMatrix  

        NumberOfGlobalDofs = assembler.GetNumberOfGlobalDofs()          # number of dofs                                             
        self.UnknownPressures = arange(0,NumberOfGlobalDofs).reshape(NumberOfGlobalDofs,1)          # unknown pressures        
        self.UnknownPressures = delete(self.UnknownPressures, s_[self.PrescribedPressures[:,0]], axis=0)
        PressuresMatrix = zeros((NumberOfGlobalDofs, self.NumberOfIncrements))                                  
        self.p = zeros((NumberOfGlobalDofs,1))
        self.pt = zeros((NumberOfGlobalDofs,1))
        self.ptt = zeros((NumberOfGlobalDofs,1))
        self.dp = zeros((NumberOfGlobalDofs,1))
        self.ddp = zeros((NumberOfGlobalDofs,1))
        self.dpt = zeros((NumberOfGlobalDofs,1))
        self.ddpt = zeros((NumberOfGlobalDofs,1))
        self.fe = zeros((NumberOfGlobalDofs,1))
        self.fet = zeros((NumberOfGlobalDofs,1))
        self.dfe = zeros((NumberOfGlobalDofs,1))
        self.dfet = zeros((NumberOfGlobalDofs,1))
        self.fi = zeros((NumberOfGlobalDofs,1))
        self.fit = zeros((NumberOfGlobalDofs,1))
        self.sumv = zeros((NumberOfGlobalDofs,1))
        sumvbk = zeros((NumberOfGlobalDofs,1))
        nonLinear = False
        for el in self.NetworkMesh.Elements:
            if el.IsNonLinear() == True:
                nonLinear = True
                break

        while self.IncrementNumber<=self.NumberOfIncrements:                              
            icc = (self.IncrementNumber%self.TimeStepFreq)
            if icc == 0:
                icc = self.TimeStepFreq
                
            #for flow in self.BoundaryConditions.elementFlow:
            for el in self.BoundaryConditions.elementFlow:
              if self.steady == True:
                  self.Flow = assembler.BoundaryConditions.GetSteadyFlow(el, self.TimeStep,icc*self.TimeStep)
              else:  
                  self.Flow = assembler.BoundaryConditions.GetTimeFlow(el, icc*self.TimeStep)
              self.fe[assembler.FlowDof[el.Id]]= self.Flow

            CoeffRelax = 0.9
            nltol = self.nltol
            self.pi = None
            pI = None
            sumvbk[:,:] = self.sumv[:,:]
            counter = 0
            while True:
                #Build the algebric equation system for the increment       
                SystemMatrix = (2.0/self.TimeStep)*self.SecondOrderGlobalMatrix + self.FirstOrderGlobalMatrix + (self.TimeStep/2.0)*self.ZeroOrderGlobalMatrix    #system matrix
                RightVector = self.fe + (2.0/self.TimeStep)*dot(self.SecondOrderGlobalMatrix,(self.pt)) + dot(self.SecondOrderGlobalMatrix,(self.dpt)) - dot(self.ZeroOrderGlobalMatrix,(self.sumv))-(self.TimeStep/2.0)*dot(self.ZeroOrderGlobalMatrix,(self.pt)) # right hand side vector                
                #The reduced (partioned) system of equations is generated.    
                RightVector[:,:] = RightVector[:,:] - dot(SystemMatrix[:,self.PrescribedPressures[:,0]],self.PrescribedPressures[:,1:])
                SystemMatrix = SystemMatrix[:,s_[self.UnknownPressures[:,0]]]
                if SystemMatrix.shape[0]> 0.0:     
                    SystemMatrix = SystemMatrix[s_[self.UnknownPressures[:,0]],:]
                RightVector = RightVector[s_[self.UnknownPressures[:,0]],:]
                #Unknown nodal point values are solved from this system.
                #  Prescribed nodal values are inserted in the solution vector.
                Solution = solve(SystemMatrix,RightVector) # solutions, unknown pressures            
                self.p[self.UnknownPressures,0] = Solution[:,:]          
                self.p[self.PrescribedPressures[:,0],0] = self.PrescribedPressures[:,1]
                #Calculating derivatives.
                #Calculating internal nodal flow values.
                self.dp = dot((2.0/self.TimeStep),(self.p-self.pt)) - self.dpt
                self.ddp = dot((4.0/self.SquareTimeStep),(self.p-self.pt)) - dot((4.0/self.TimeStep),self.dpt) -self.ddpt
                self.sumv = sumvbk + dot((self.TimeStep/2.0),(self.pt+self.p))
                self.fi = dot(self.SecondOrderGlobalMatrix,(self.dp)) + dot(self.FirstOrderGlobalMatrix,(self.p)) + dot(self.ZeroOrderGlobalMatrix,(self.sumv))             
                if not nonLinear :
                    break
                
                if self.pi == None:
                    self.pi = zeros((NumberOfGlobalDofs,1))
                    self.pi[:,:] = self.pt[:,:]
                pI = CoeffRelax * self.p + self.pi * (1.0-CoeffRelax)
                self.p[:,:] = pI[:,:]
                den = norm(self.pi,Inf)
                if den < 1e-12:
                    den = 1.0
                nlerr = norm(self.p-self.pi,Inf) / den
              
                info = {'dofmap':assembler.DofMap,'solution':[self.p, self.pt, self.ptt],'incrementNumber':self.IncrementNumber,'history':history}
                self.Evaluator.SetInfo(info)

                assembler.Assemble(self.SimulationContext, self.Evaluator, self.LinearZeroOrderGlobalMatrix, self.LinearFirstOrderGlobalMatrix, self.LinearSecondOrderGlobalMatrix)
                self.ZeroOrderGlobalMatrix = assembler.ZeroOrderGlobalMatrix
                self.FirstOrderGlobalMatrix = assembler.FirstOrderGlobalMatrix
                self.SecondOrderGlobalMatrix = assembler.SecondOrderGlobalMatrix        
                
                #Dynamic nonlinear relaxing coefficient 
                if counter == 100:
                    print "relaxing..."
                    print nlerr, nltol, CoeffRelax
                    counter = 0
                    self.pi[:,:] = None
                    self.sumv[:,:] = sumvbk[:,:]
                    CoeffRelax *= 0.6
                    nltol *= 0.95
                if nlerr < nltol:
                    nltol = self.nltol
                    counter = 0 
                    break
                counter+=1
                self.pi[:,:] = self.p[:,:]
                                 
            self.ptt[:,:] = self.pt[:,:]
            self.pt[:,:] = self.p[:,:]
            self.dpt[:,:] = self.dp[:,:]
            self.ddpt[:,:] = self.ddp[:,:]
            self.fet[:,:] = self.fe[:,:]
            self.fit[:,:] = self.fi[:,:]                
            PressuresMatrix[:,(self.IncrementNumber-1)] = self.p[:,0]  
            history.insert(0,self.IncrementNumber)
            history = history[:3]
            
            if self.steady == True:
                self.MinimumIncrementNumber = 0.01* self.NumberOfIncrements
                if norm(self.fi-self.fe,Inf)<self.convergence and self.IncrementNumber > self.MinimumIncrementNumber:
                    self.IncrementNumber = self.NumberOfIncrements 
                else:
                    pass

            if self.IncrementNumber==ceil(0.05*self.NumberOfIncrements):
                print "->5%"
            if self.IncrementNumber==ceil(0.25*self.NumberOfIncrements):
                print "->25%"
            if self.IncrementNumber==ceil(0.5*self.NumberOfIncrements):
                print "->50%"   
            if self.IncrementNumber==ceil(0.70*self.NumberOfIncrements):
                print "->70%" 
            if self.IncrementNumber==ceil(0.90*self.NumberOfIncrements):
                print "->90%"  
            if self.IncrementNumber==ceil(0.99*self.NumberOfIncrements):
                print "->99%"   
            
            self.IncrementNumber = self.IncrementNumber+1
            self.EndIncrementTime = self.EndIncrementTime + self.TimeStep    # increment 
        info = {'dofmap':assembler.DofMap,'solution':[self.p, self.pt, self.ptt],'incrementNumber':self.IncrementNumber,'history':history,'allSolution':PressuresMatrix}      
        self.Evaluator.SetInfo(info)
        self.Solutions = PressuresMatrix
        return PressuresMatrix