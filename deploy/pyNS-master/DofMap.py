#!/usr/bin/env python

## Program:   PyNS
## Module:    DofMap.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

class DofMap(object):
    '''
    This class maps local degrees of freedom (dof) into global dof.
    Each element knows only its global nodeIds and its local dofs.
    DofMap must consider also various connections between each element.
    This class provides the following methods:
    SetNetworkMesh: a method for setting NetworkMesh input.
    GetDof: a method which returns global dof from (element.Id, localdof).
    Build: a method for building DofMap.
    DofMapOutput: a method for printing DofMap Output ordered by element's Id.
    '''
   
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkMesh = None
        self.DofMap = {}
        self.NumberOfGlobalDofs = 0

    def SetNetworkMesh(self,networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh

    def GetDof(self,elementId,localDof):
        '''
        This method calculates Global dof from (element.Id, localdof)
        '''
        return self.DofMap[(elementId,localDof)]

    def Build(self):
        '''
        Building DofMap, mapping local dofs of each element into global network dofs
        '''
        nodesToElements = self.NetworkMesh.BuildNodesToElements()
        newDofId = 0
        for element in self.NetworkMesh.Elements:
            dofNodes = element.GetDofNodes()
            matchedLocalDofs = []
            for dofNode in dofNodes:
                nodeElementIds = nodesToElements[dofNode]
                aNeighborElementId = None
                for nodeElementId in nodeElementIds:
                    if element.Id == nodeElementId: 
                        break
                    else:
                        aNeighborElementId = nodeElementId
                    for el in self.NetworkMesh.Elements:
                        if el.Id == aNeighborElementId:  
                            neighborElement = el
                            localDof = element.GetLocalDof(dofNode) 
                            matchedLocalDofs.append(localDof)
                            aNeighborLocalDof = neighborElement.GetLocalDof(dofNode)
                            self.DofMap[(element.Id,localDof)] = self.DofMap[aNeighborElementId,aNeighborLocalDof]        
            for localDof in element.dof:
                if localDof in matchedLocalDofs:
                    continue
                self.DofMap[(element.Id,localDof)] = newDofId
                newDofId += 1
        self.NumberOfGlobalDofs = newDofId        
        return DofMap
    
    def DofMapOutput(self):
        '''
        This method returns DofMap Output ordered by element's IDs.
        '''
        for i in sorted(self.DofMap.items(),
        key=lambda x:int(x[1])):
            print "%-10s %3s" % i