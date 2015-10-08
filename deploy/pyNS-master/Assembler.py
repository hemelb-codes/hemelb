#!/usr/bin/env python

## Program:   PyNS
## Module:    Assembler.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from DofMap import DofMap
from numpy.core.numeric import zeros
from numpy.numarray.numerictypes import Int32

class Assembler(object):
    '''
    Assembler Class generates Global Matrices from Local Matrices and dofmap.
    This class provides the following methods:
    SetNetworkMesh: a method for setting NetworkMesh input.
    SetBoundaryConditions: a method for setting BoundaryConditions input.
    GetNumberOfGlobalDofs: a method for calculating number of global degrees of freedom.
    Assemble: a method for building boundary conditions vectors and assembling
    each local Zero, First and Second Order matrix into global Zero, First and Second Order matrix.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkMesh = None
        self.BoundaryConditions = None
        self.PrescribedPressures = None
        self.Flow = None
        self.FlowDof = {}
        self.DofMap = None
        self.LinearZeroOrderGlobalMatrix = None
        self.LinearFirstOrderGlobalMatrix = None
        self.LinearSecondOrderGlobalMatrix = None
        self.ZeroOrderGlobalMatrix = None
        self.FirstOrderGlobalMatrix = None
        self.SecondOrderGlobalMatrix = None
        self.Evaluator = None
        self.Initialized = False
        
    def SetNetworkMesh(self, networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetBoundaryConditions(self, boundaryConditions):
        '''
        Setting BoundaryConditions
        '''
        self.BoundaryConditions = boundaryConditions

    def GetNumberOfGlobalDofs(self):
        '''
        This method returns number of global degrees of freedom
        '''
        return self.DofMap.NumberOfGlobalDofs
    
    def AssembleBoundaryConditions(self, simulationContext):
        '''
        This method assembles arrays of prescribed pressures
        '''
        self.DofMap = DofMap()
        self.DofMap.SetNetworkMesh(self.NetworkMesh)
        self.DofMap.Build()
             
        # Searching for prescribed pressure output/input
        numberOfElements = 0
        for element in self.NetworkMesh.Elements:  
            for dof in element.GetExternalPressureLocalDofs():
                    numberOfElements+=1
        if self.BoundaryConditions.OutP is not None:
            numberOfElements+=1
        if self.BoundaryConditions.InP is not None:
            numberOfElements+=1
      
            
        # Setting Transmural Pressures for windkessel elements and wave propagation elements
        PrescribedPressures = zeros((numberOfElements,2))
        done = 0
        i = 0
        for element in self.NetworkMesh.Elements:
            for dof in element.GetExternalPressureLocalDofs():    
                if self.BoundaryConditions.OutP is not None:
                    if element.Id == self.BoundaryConditions.elementOut.Id:
                        if done == 0:
                            value1 = self.DofMap.DofMap[(self.BoundaryConditions.elementOut.Id,self.BoundaryConditions.elementOut.GetLocalDof(int(self.BoundaryConditions.NodeOut)))]
                            value2 = self.BoundaryConditions.OutP
                            PrescribedPressures[i,0] = value1
                            PrescribedPressures[i,1] = value2
                            done = 1
                            i+=1
                if self.BoundaryConditions.InP is not None:
                    if element.Id == self.BoundaryConditions.elementIn.Id:
                        if done == 0:
                            value1 = self.DofMap.DofMap[(self.BoundaryConditions.elementIn.Id,self.BoundaryConditions.elementIn.GetLocalDof(int(self.BoundaryConditions.NodeIn)))]
                            value2 = self.BoundaryConditions.InP
                            PrescribedPressures[i,0] = value1
                            PrescribedPressures[i,1] = value2
                            done = 1
                            i+=1
                            
                value1 = self.DofMap.DofMap[(element.Id,dof)]
                value2 = self.BoundaryConditions.PressureValues[element.Id]
                PrescribedPressures[i,0] = value1
                PrescribedPressures[i,1] = value2     
                i+=1
        
        self.BoundaryConditions.SetSimulationContext(simulationContext)
        for el in self.BoundaryConditions.elementFlow:
            self.FlowDof[el.Id] = self.DofMap.DofMap[(el.Id,el.GetLocalDof(int(self.BoundaryConditions.NodeFlow[el.Id])))]
        self.PrescribedPressures = PrescribedPressures.astype(Int32)
        return self.PrescribedPressures

    def AssembleInit(self, simulationContext, evaluator):
        '''
        This method calculates:
        Pressures GlobalDofs and Input Flow Vector from BoundaryConditions (Prescribed Pressures matrix, each row ---> Dof, Value).
        Zero, First and Second Order Global Matrix from Local Matrix of each element.
        '''
         
        #Assembling global matrices from local matrices.
        self.LinearZeroOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
        self.LinearFirstOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
        self.LinearSecondOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
        self.ZeroOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
        self.FirstOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
        self.SecondOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
        
        nonLinear = False
        #Building Global Linear Matrices      
        for element in self.DofMap.NetworkMesh.Elements:
            if element.Initialized == False:             
                element.Initialize(simulationContext)
            if element.nonLinear == False: 
                element.InputParameters(evaluator)
                zeroOrderMatrix = element.GetZeroOrderMatrix()
                firstOrderMatrix = element.GetFirstOrderMatrix()
                secondOrderMatrix = element.GetSecondOrderMatrix()
                for rowLocalDof in element.dof:
                    rowGlobalDof = self.DofMap.GetDof(element.Id,rowLocalDof)
                    for columnLocalDof in element.dof:
                        columnGlobalDof = self.DofMap.GetDof(element.Id,columnLocalDof)
                        self.LinearZeroOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += zeroOrderMatrix[rowLocalDof,columnLocalDof]   
                        self.LinearFirstOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += firstOrderMatrix[rowLocalDof,columnLocalDof]
                        self.LinearSecondOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += secondOrderMatrix[rowLocalDof,columnLocalDof]
            else:
                nonLinear = True
            
        if nonLinear == True:        
            self.ZeroOrderGlobalMatrix,self.FirstOrderGlobalMatrix,self.SecondOrderGlobalMatrix = \
            self.Assemble(simulationContext, evaluator, self.LinearZeroOrderGlobalMatrix, self.LinearFirstOrderGlobalMatrix, self.LinearSecondOrderGlobalMatrix)
        
        else:
            self.ZeroOrderGlobalMatrix[:,:] = self.LinearZeroOrderGlobalMatrix[:,:]
            self.FirstOrderGlobalMatrix[:,:] = self.LinearFirstOrderGlobalMatrix[:,:]
            self.SecondOrderGlobalMatrix[:,:] = self.LinearSecondOrderGlobalMatrix[:,:]
        
        return self.LinearZeroOrderGlobalMatrix, self.LinearFirstOrderGlobalMatrix, self.LinearSecondOrderGlobalMatrix
        
        
    def Assemble(self, simulationContext, evaluator, linearZeroOrder, linearFirstOrder, linearSecondOrder):  
        '''
        This method builds non-linear elements of the global matrices from linear matrices.
        '''
        
        #Building non linear matrices
        self.ZeroOrderGlobalMatrix[:,:] = linearZeroOrder[:,:]  
        self.FirstOrderGlobalMatrix[:,:] = linearFirstOrder[:,:]
        self.SecondOrderGlobalMatrix[:,:] = linearSecondOrder[:,:]  
              
        for element in self.DofMap.NetworkMesh.Elements:
            if element.nonLinear == True:
                element.Initialize(simulationContext)
                element.InputParameters(evaluator)
                zeroOrderMatrix = element.GetZeroOrderMatrix()
                firstOrderMatrix = element.GetFirstOrderMatrix()
                secondOrderMatrix = element.GetSecondOrderMatrix()
                for rowLocalDof in element.dof:
                    rowGlobalDof = self.DofMap.GetDof(element.Id,rowLocalDof)
                    for columnLocalDof in element.dof:
                        columnGlobalDof = self.DofMap.GetDof(element.Id,columnLocalDof)
                        self.ZeroOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += zeroOrderMatrix[rowLocalDof,columnLocalDof]   
                        self.FirstOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += firstOrderMatrix[rowLocalDof,columnLocalDof]
                        self.SecondOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += secondOrderMatrix[rowLocalDof,columnLocalDof]
        return self.ZeroOrderGlobalMatrix, self.FirstOrderGlobalMatrix, self.SecondOrderGlobalMatrix     