#!/usr/bin/env python

## Program:   PyNS
## Module:    MeshGenerator.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

import Elements
from NetworkMesh import Entity
from numpy.ma.core import ceil

class MeshGenerator(object):
    '''
    This Class builds NetworkMesh from NetworkGraph,
    providing two different type of meshing (Maximum length
    or Radius Tolerance).
    This Class provides the following methods:
    SetMeshType: A method for setting specific MeshType for each edge.
    SetNetworkGraph: a method for setting NetworkGraph input.
    SetNetworkMesh: a method for setting NetworkMesh input.
    SetTolerance: a method for setting Radius Tolerance % value.
    SetMaxLength: a method for setting maximum length for a single mesh.
    ParsingEdgeProperties: a method for parsing properties for each edge of the NetworkGraph.
    GenerateMesh: a method for generating NetworkMesh from NetworkGraph. If MaxLength is specified,
    GenerateMesh uses MaxLengthMeshing method, otherwise ToleranceMeshing method is used.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.MeshType = {} # edgeId:elementType   (default 0D)
        self.NetworkGraph = None
        self.NetworkMesh = None
        self.MaxLength = None
        self.Tolerance = None

    def SetMeshType(self, edgeId, elementType):
        '''
        Setting meshType for each edge
        '''
        self.MeshType[edgeId] = elementType
        
    def SetNetworkGraph(self,networkGraph):
        '''
        Setting NetworkGraph
        '''
        self.NetworkGraph = networkGraph
        
    def SetNetworkMesh(self, networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
    
    def SetTolerance(self, radiusTolerance):
        '''
        Setting Tolerance %
        '''
        self.Tolerance = {'radius':radiusTolerance*1e-2}
    
    def SetMaxLength(self, maxLength):
        '''
        Setting Maximum segment length
        '''
        self.MaxLength = maxLength
    
    def ParsingEdgeProperties(self, edge):
        '''
        Parsing NetworkGraph.Edge properties into 
        a dictionary edgeproperty {name:value}
        '''
        edgeProperty = {} #dictionary property:{s:value}
        
        if edge.Type is not None:
            self.SetMeshType(edge.Id, edge.Type)
    
        if edge.Length.has_key('array') == False:
            length = edge.Length['value']
            edgeProperty['length'] = length

        if edge.Radius.has_key('array') == False:                      
            try:
                radius = edge.Radius['value']
                meshradius = {}
                meshradius[0.0] = radius
                meshradius[1.0] = radius
                edgeProperty['radius'] = meshradius
            except KeyError:
                edgeProperty['radius'] = None
                
            if edge.xRadius.has_key('value') == True:
                xradius = edge.xRadius['value']
            else:
                xradius = None
            edgeProperty['xradius'] = {0.0:xradius, 1.0:xradius}
            if edge.yRadius.has_key('value') == True:
                yradius = edge.yRadius['value']
            else:
                yradius = None
            edgeProperty['yradius'] = {0.0:yradius, 1.0:yradius}  
          
        if edge.Radius.has_key('array') == True: 
            meshradius = {}
            s_list = []
            for s in edge.Radius['array'].iterkeys():
                s_list.append(s)			
                s_list.sort()
                if meshradius.has_key(s):
                    meshradius[s].append(edge.Radius['array'][s])    
                else:
                    meshradius[s]=edge.Radius['array'][s]
                edgeProperty['radius'] = meshradius
        
        if edge.xRadius.has_key('array') == True:
            meshradius = {}
            s_list = []
            for s in edge.xRadius['array'].iterkeys():
                s_list.append(s)
                s_list.sort()
                if meshradius.has_key(s):
                    meshradius[s].append(edge.xRadius['array'][s])    
                else:
                    meshradius[s]=edge.xRadius['array'][s]
                edgeProperty['xradius'] = meshradius
                
        elif edge.xRadius.has_key('value') == True:
            xradius = edge.xRadius['value']
            edgeProperty['xradius'] = {0.0:xradius, 1.0:xradius}  
        else:
            xradius = None
            edgeProperty['xradius'] = {0.0:xradius, 1.0:xradius}     
        
        if edge.yRadius.has_key('array') == True:
            meshradius = {}
            s_list = []
            for s in edge.yRadius['array'].iterkeys():
                s_list.append(s)
                s_list.sort()
                if meshradius.has_key(s):
                    meshradius[s].append(edge.yRadius['array'][s])    
                else:
                    meshradius[s]=edge.yRadius['array'][s]
                edgeProperty['yradius'] = meshradius
                
        elif edge.yRadius.has_key('value') == True:
            yradius = edge.yRadius['value']
            edgeProperty['yradius'] = {0.0:yradius, 1.0:yradius}
        else:
            yradius = None
            edgeProperty['yradius'] = {0.0:yradius, 1.0:yradius}     
                      
        if edge.WallThickness.has_key('value') == True:
            wallThickness = edge.WallThickness['value'] 
            meshwt = {}
            meshwt[0.0] = wallThickness
            meshwt[1.0] = wallThickness
            edgeProperty['wallthickness'] = meshwt          
                  
        if edge.WallThickness.has_key('array') == True:
            meshwt = {}
            s_list = []
            for s in edge.WallThickness['array'].iterkeys():
                s_list.append(s)
                s_list.sort()
                if meshwt.has_key(s):
                    meshwt[s].append(edge.WallThickness['array'][s])     
                else:
                    meshwt[s]=edge.WallThickness['array'][s]
            edgeProperty['wallthickness'] = meshwt 
            
        if edge.WallThickness.has_key('expression') == True:
            edgeProperty['wallthickness'] = edge.WallThickness['expression']
            
        if edge.YoungModulus.has_key('array') == False:
            try:
                youngModulus = edge.YoungModulus['value']
            except:
                youngModulus = 0.0
            meshym = {}
            meshym[0.0] = youngModulus
            meshym[1.0] = youngModulus
            edgeProperty['youngmodulus'] = meshym
        
        if edge.ScalarRadius != {}:
            edgeProperty['ScalarRadius'] = edge.ScalarRadius
        else:
            edgeProperty['ScalarRadius'] = {}
        
        if edgeProperty['radius'] is not None:
            for rad in edgeProperty['radius'].itervalues():
                if type(rad) is not str:
                    edgeProperty['ScalarRadius'] = edgeProperty['radius']
        
        return edgeProperty

    def GenerateMesh(self):
        '''
        This method generates Network Mesh from Network Graph
        '''
        self.NetworkMesh.Id = self.NetworkGraph.Id    
        if self.Tolerance is not None:   
            self.ToleranceMeshing(self.Tolerance)
        else:
            if self.MaxLength is not None:
                self.MaxLengthMeshing(self.MaxLength)
            else:
                raise MeshingError   
    
    def MaxLengthMeshing(self, maxLength):
        '''
        This method build NetworkMesh from NetworkGraph using Maximum Length meshing criterium.
        If edge's length is greater than maxLength parameter, a new mesh is created.
        '''
        meshId = 1 #First mesh Id index
        self.NetworkMesh.meshToEdges = {} # dict edgenode:meshnode
        self.NetworkMesh.s_mesh = {}  # dict (s,edge):meshnode
        self.NetworkMesh.MeshToGraph = {} #  {meshId : edgeId}
        self.NetworkMesh.MeshToS = {}  #  {meshId : [s]}
        self.NetworkMesh.Entities = {} # {entity:[mesh]}       
        self.NetworkMesh.GraphEdgeToMesh = {} #{edge:[meshId:[s0,s1]]}
        self.NetworkMesh.ElementIdsToElements = {} # {mesh.Id:mesh}
        self.NetworkMesh.GraphNodeToMesh = {} #{edge:[meshId:[s0,s1]]}
        nLeak = 0
        meshNode1 = None
        meshNode2 = None
        nodes_list = []  #nodes list
        nodes_list_last = len(nodes_list)  
        meshNodeAnast0 = None
        meshNodeAnast1 = None
        meshNodeAnast2 = None
        done = 0  #if = 1, AVF resistances meshed
        for edgeId, edge in self.NetworkGraph.Edges.iteritems():
            if self.NetworkMesh.GraphEdgeToMesh.has_key(edge):
                pass
            else:
                edgeProperties = self.ParsingEdgeProperties(edge)
                # Searching for anastomosis element (arterial and venous side)
                for nodeId, node in self.NetworkGraph.Nodes.iteritems():
                    if node.Type == 'anastomosis':
                        if int(nodeId) in edge.NodeIds:                          
                            if edge.Side == 'venous':  
                                if done == 0:  
                                    try:
                                        meshNode0 = self.NetworkMesh.meshToEdges[edge.NodeIds[0]]
                                    except KeyError:
                                        meshNode0 = nodes_list_last+1
                                    meshNodeAnast0 = meshNode0  
                                    if meshNode0 not in nodes_list:
                                        nodes_list.append(meshNode0)
                                        nodes_list_last = len(nodes_list)   
                                    meshNode1 = nodes_list_last+1
                                    if meshNode1 not in nodes_list:
                                        nodes_list.append(meshNode1)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast1 = meshNode1
                                    meshNode2 = nodes_list_last+1
                                    if meshNode2 not in nodes_list:
                                        nodes_list.append(meshNode2)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast2 = meshNode2   
                                else:
                                    meshNode0 = meshNodeAnast0
                                    meshNode1 = meshNodeAnast1
                                    meshNode2 = meshNodeAnast2                                
                                self.NetworkMesh.meshToEdges[edge.NodeIds[0]] = meshNodeAnast2        
                            if edge.Side == 'arterial' and  int(nodeId) == edge.NodeIds[1]:
                                if done == 0:
                                    try:
                                        meshNode0 = self.NetworkMesh.meshToEdges[edge.NodeIds[1]]
                                    except KeyError:
                                        meshNode0 = nodes_list_last+1                                               
                                    meshNodeAnast0 = meshNode0
                                    if meshNode0 not in nodes_list:
                                        nodes_list.append(meshNode0)
                                        nodes_list_last = len(nodes_list)
                                    meshNode1 = nodes_list_last+1
                                    if meshNode1 not in nodes_list:
                                        nodes_list.append(meshNode1)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast1 = meshNode1
                                    meshNode2 = nodes_list_last+1
                                    if meshNode2 not in nodes_list:
                                        nodes_list.append(meshNode2)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast2 = meshNode2 
                                else:
                                    meshNode0 = meshNodeAnast0
                                    meshNode1 = meshNodeAnast1
                                    meshNode2 = meshNodeAnast2
                                self.NetworkMesh.meshToEdges[edge.NodeIds[1]] = meshNodeAnast0       
                            if edge.Side == 'arterial' and  int(nodeId) == edge.NodeIds[0]:
                                if done == 0:                                   
                                    try:
                                        meshNode1 = self.NetworkMesh.meshToEdges[edge.NodeIds[0]]
                                    except KeyError:
                                        meshNode1 = nodes_list_last+1                                               
                                    meshNodeAnast1 = meshNode1       
                                    if meshNode1 not in nodes_list:
                                        nodes_list.append(meshNode1)
                                        nodes_list_last = len(nodes_list)
                                    meshNode0 = nodes_list_last+1
                                    if meshNode0 not in nodes_list:
                                        nodes_list.append(meshNode0)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast0 = meshNode0
                                    meshNode2 = nodes_list_last+1
                                    if meshNode2 not in nodes_list:
                                        nodes_list.append(meshNode2)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast2 = meshNode2 
                                else:                                   
                                    meshNode0 = meshNodeAnast0
                                    meshNode1 = meshNodeAnast1
                                    meshNode2 = meshNodeAnast2                                 
                                self.NetworkMesh.meshToEdges[edge.NodeIds[0]] = meshNodeAnast1                                                               
                            #NEW ELEMENT
                            if done == 0:  
                                anast_res_0_1 = node.Properties['arterial_resistance']
                                anast_res_0_2 = node.Properties['venous_resistance']
                                elementParameters = {}
                                elementParameters['resistance_0_1'] = anast_res_0_1
                                elementParameters['resistance_0_2'] = anast_res_0_2
                                newElement = Elements.Anastomosis(str(meshId), [meshNode0,meshNode1,meshNode2], elementParameters, name = 'anastomosis')
                                self.NetworkMesh.Elements.append(newElement)
                                self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement
                                self.NetworkMesh.GraphNodeToMesh[node] = newElement.Id
                                meshId+=1
                                entity = Entity()
                                entity.SetId(node.Name)   
                                if self.NetworkMesh.Entities.has_key(entity):
                                    self.NetworkMesh.Entities[entity].append(newElement)
                                else:
                                    self.NetworkMesh.Entities[entity] = [newElement]
                                done = 1                                          
                #MESHNODE1########
                if meshNode1 == None:
                    meshNode1 = 1
                else:
                    try:
                        meshNode1 = self.NetworkMesh.meshToEdges[edge.NodeIds[0]]
                    except KeyError:
                        meshNode1 = nodes_list_last+1      
                ##################      
                if edgeProperties['length'] <= self.MaxLength:
                    elementParameters = {}
                    self.NetworkMesh.meshToEdges[edge.NodeIds[0]] = meshNode1
                    if meshNode1 not in nodes_list:
                            nodes_list.append(meshNode1)
                            nodes_list_last = len(nodes_list)
                    #MESHNODE2#######
                    try:
                        meshNode2 = self.NetworkMesh.meshToEdges[edge.NodeIds[1]]      
                    except KeyError:
                        meshNode2 = nodes_list_last+1                   
                    self.NetworkMesh.meshToEdges[edge.NodeIds[1]] = meshNode2
                    if meshNode2 not in nodes_list:
                            nodes_list.append(meshNode2)
                            nodes_list_last = len(nodes_list)
                    ###############                    
                    #parameters####
                    
                    elementParameters['s1'] = 0.0 
                    elementParameters['s2'] = 1.0
                    elementParameters['length'] = edgeProperties['length']
                    if self.MeshType.get(edgeId,None) == "Resistance":
                        pass
                    else:
                        if edgeProperties['radius'] is not None:
                            if type(edgeProperties['radius'][0]) is str or type(edgeProperties['radius'][1]) is str:
                                elementParameters['radiusAtRest'] = edgeProperties['ScalarRadius']
                                elementParameters['radius'] = edgeProperties['radius']
                            else:
                                elementParameters['radius'] = edgeProperties['radius']
                        else:
                            elementParameters['radius'] = {elementParameters['s1']:edgeProperties['xradius'][elementParameters['s1']]*edgeProperties['yradius'][elementParameters['s1']]**0.5, elementParameters['s2']:edgeProperties['xradius'][elementParameters['s2']]*edgeProperties['yradius'][elementParameters['s2']]**0.5}
                        elementParameters['wall_thickness'] = edgeProperties['wallthickness']
                        elementParameters['young_modulus'] = edgeProperties['youngmodulus']                             
                        if edgeProperties['xradius'][0.0] is not None:     
                            elementParameters['xradius'] = edgeProperties['xradius']
                        else:
                            elementParameters['xradius'] = None
                        if edgeProperties['yradius'][0.0] is not None:     
                            elementParameters['yradius'] = edgeProperties['yradius']
                        else:
                            elementParameters['yradius'] = None  
                        
                    
                    elementParameters["resistance"] = edge.Resistance

                    if edge.NlCompliance != {}:
                        elementParameters["compliance"] = edge.NlCompliance['expression']
                    else:
                        elementParameters["compliance"] = edge.Compliance
                    elementParameters["leakage"] = edge.QLeakage
                    if edge.QLeakage:
                        nLeak+=1
                    ##############################################
                    # MESHING AN EDGE ACCORDING TO ITS MESH TYPE #
                    ##############################################
                    
                    newElement = Elements.NewElement(self.MeshType.get(edgeId,None), str(meshId), [meshNode1,meshNode2], elementParameters, edge.Side, edge.Name)                                       
                                 
                    self.NetworkMesh.Elements.append(newElement)                
                    #Dicts###
                    self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement    
                    self.NetworkMesh.s_mesh[(0.0,edge)] = meshNode1
                    self.NetworkMesh.s_mesh[(1.0,edge)] = meshNode2
                    self.NetworkMesh.MeshToGraph[str(meshId)] = edge.Id
                    self.NetworkMesh.MeshToS[meshId] = [0.0,1.0]
                    if self.NetworkMesh.GraphEdgeToMesh.has_key(edge):                        
                        meshtoSappend = {meshId:self.NetworkMesh.MeshToS[meshId]} 
                        self.NetworkMesh.GraphEdgeToMesh[edge].append(meshtoSappend)
                    else:
                        meshtoSappend = {meshId:self.NetworkMesh.MeshToS[meshId]}
                        self.NetworkMesh.GraphEdgeToMesh[edge] = [meshtoSappend]   
                    ##############    
                    meshId+=1         
                else:
                    elementParameters = {}
                    numEl = int(ceil(edgeProperties['length'] / self.MaxLength))         
                    i = 0.0
                    nameId = 1                   
                    self.NetworkMesh.meshToEdges[edge.NodeIds[0]] = meshNode1                   
                    if meshNode1 not in nodes_list:
                            nodes_list.append(meshNode1)
                            nodes_list_last = len(nodes_list)                            
                    while i < numEl:
                        s1 = i / numEl
                        s2 = (i+1) / numEl 
                        elLength = edgeProperties['length']/numEl                                                      
                        if s2 >= 1.0:
                            s2 = 1.0
                            try:
                                meshNode2 = self.NetworkMesh.meshToEdges[edge.NodeIds[1]]
                            except KeyError:
                                meshNode2 = nodes_list_last+1
                            self.NetworkMesh.meshToEdges[edge.NodeIds[1]] = meshNode2                 
                        else:
                            meshNode2 = nodes_list_last+1    
                                           
                        elementParameters['s1'] = s1 
                        elementParameters['s2'] = s2 
                        if edgeProperties['radius'] is not None:
                            if type(edgeProperties['radius'][0]) is str or type(edgeProperties['radius'][1]) is str:
                                elementParameters['radiusAtRest'] = {s1:edgeProperties['ScalarRadius'][0]-((edgeProperties['ScalarRadius'][0]-edgeProperties['ScalarRadius'][1])*s1), s2:edgeProperties['ScalarRadius'][0]-((edgeProperties['ScalarRadius'][0]-edgeProperties['ScalarRadius'][1])*s2)}
                                elementParameters['radius'] = {s1:edgeProperties['radius'][0], s2:edgeProperties['radius'][1]}
                            else: 
                                elementParameters['radius'] = {s1:edgeProperties['radius'][0]-((edgeProperties['radius'][0]-edgeProperties['radius'][1])*s1), s2:edgeProperties['radius'][0]-((edgeProperties['radius'][0]-edgeProperties['radius'][1])*s2)}
                        else:
                            rs1 = ((edgeProperties['xradius'][0]-((edgeProperties['xradius'][0]-edgeProperties['xradius'][1])*s1))*(edgeProperties['yradius'][0]-((edgeProperties['yradius'][0]-edgeProperties['yradius'][1])*s1)))**0.5
                            rs2 = ((edgeProperties['xradius'][0]-((edgeProperties['xradius'][0]-edgeProperties['xradius'][1])*s2))*(edgeProperties['yradius'][0]-((edgeProperties['yradius'][0]-edgeProperties['yradius'][1])*s2)))**0.5
                            elementParameters['radius'] = {s1:rs1, s2:rs2}   
                        if edgeProperties['xradius'][0.0] is not None:
                            elementParameters['xradius'] = {s1:edgeProperties['xradius'][0]-((edgeProperties['xradius'][0]-edgeProperties['xradius'][1])*s1), s2:edgeProperties['xradius'][0]-((edgeProperties['xradius'][0]-edgeProperties['xradius'][1])*s2)}   
                        else:
                            elementParameters['xradius'] = None
                        if edgeProperties['yradius'][0.0] is not None:
                            elementParameters['yradius'] = {s1:edgeProperties['yradius'][0]-((edgeProperties['yradius'][0]-edgeProperties['yradius'][1])*s1), s2:edgeProperties['yradius'][0]-((edgeProperties['yradius'][0]-edgeProperties['yradius'][1])*s2)}
                        else:
                            elementParameters['yradius'] = None                            
                        if type(edgeProperties['wallthickness']) == dict:
                            elementParameters['wall_thickness'] = {s1:edgeProperties['wallthickness'][0]-((edgeProperties['wallthickness'][0]-edgeProperties['wallthickness'][1])*s1), s2:edgeProperties['wallthickness'][0]-((edgeProperties['wallthickness'][0]-edgeProperties['wallthickness'][1])*s2)}
                        else:                         
                            elementParameters['wall_thickness']=edgeProperties['wallthickness']                         
                        elementParameters['young_modulus'] = {s1:edgeProperties['youngmodulus'][0]-((edgeProperties['youngmodulus'][0]-edgeProperties['youngmodulus'][1])*s1), s2:edgeProperties['youngmodulus'][0]-((edgeProperties['youngmodulus'][0]-edgeProperties['youngmodulus'][1])*s2)}
                        
                        elementParameters['length'] = elLength
                        elementParameters["resistance"] = edge.Resistance

                        
                        if edge.NlCompliance != {}:
                            elementParameters["compliance"] = edge.NlCompliance['expression']
                        else:
                            elementParameters["compliance"] = edge.Compliance
                        
                        elementParameters["leakage"] = edge.QLeakage
                        if edge.QLeakage:
                            nLeak+=2                      
                        name = edge.Name + "_" + str(nameId)                                        
                        #############                        
                         
                        newElement= Elements.NewElement(self.MeshType.get(edgeId,None),str(meshId), [meshNode1,meshNode2], elementParameters, edge.Side, name) 
                                                 
                        self.NetworkMesh.Elements.append(newElement)
                        self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement                                           
                        self.NetworkMesh.s_mesh[(s1,edge)] = meshNode1
                        self.NetworkMesh.s_mesh[(s2,edge)] = meshNode2                                                       
                        self.NetworkMesh.MeshToGraph[str(meshId)] = edge.Id
                        self.NetworkMesh.MeshToS[meshId] = [s1,s2]                        
                        if self.NetworkMesh.GraphEdgeToMesh.has_key(edge):                        
                            meshtoSappend = {meshId:self.NetworkMesh.MeshToS[meshId]} 
                            self.NetworkMesh.GraphEdgeToMesh[edge].append(meshtoSappend)
                        else:
                            meshtoSappend = {meshId:self.NetworkMesh.MeshToS[meshId]}
                            self.NetworkMesh.GraphEdgeToMesh[edge] = [meshtoSappend]                                               
                        meshId+=1                        
                        if meshNode2 not in nodes_list:
                            nodes_list.append(meshNode2)
                            nodes_list_last = len(nodes_list)               
                        meshNode1 = meshNode2
                        nameId+=1
                        i+=1 
        #SETTING ANASTOMOSIS CONNECTIONS DICTIONARIES
        for nodeId, node in self.NetworkGraph.Nodes.iteritems():                              
            if node.Type == 'anastomosis':
                edge = self.NetworkGraph.Edges[node.Properties['proximal']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][1] == 1.0:
                            proximal = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetProximal(proximal)
                try:
                    edge = self.NetworkGraph.Edges[node.Properties['distal']]
                    for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                        for meshId in meshlist:
                            if self.NetworkMesh.MeshToS[meshId][0] == 0.0:
                                distal = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                                self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetDistal(distal)
                except KeyError:
                    pass
                edge = self.NetworkGraph.Edges[node.Properties['vein']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][0] == 0.0:
                            vein = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetVein(vein)    
        #SETTING LEAKAGES
        for el in self.NetworkMesh.Elements:
            if el.Type == "WavePropagation":
                el.Leakages = nLeak
                           
                       
        #ENTITIES DICTIONARY
        for sedge in self.NetworkGraph.SuperEdges.values():
            entity = Entity()
            entity.SetId(sedge.Name)
            for edge in sedge.Edges.values():              
                for key, value in self.NetworkMesh.MeshToGraph.iteritems():                    
                    if value == edge.Id:
                        for el in self.NetworkMesh.Elements:
                            if el.Type != 'Resistance':  
                                if el.Id == key:                             
                                    if self.NetworkMesh.Entities.has_key(entity):
                                        self.NetworkMesh.Entities[entity].append(el)
                                    else:
                                        self.NetworkMesh.Entities[entity] = [el]          
        # MESHING END SEGMENTS (WINDKESSEL ELEMENTS) 
        entity = Entity()
        entity.SetId('end_segment')        
        nodeslist = []
        for nodeId, node in self.NetworkGraph.Nodes.iteritems():
            nodeslist.append(int(node.Id))
            nodeslist.sort()
        nodes_list_last = len(nodes_list)
        meshNode2 = nodes_list_last+1
        meshEndId = 1 
        for nd in nodeslist:
            for nodeId, node in self.NetworkGraph.Nodes.iteritems():
                if nd == int(node.Id):
                    if node.Type == 'downstream network':      
                        meshNode1 = self.NetworkMesh.meshToEdges[int(node.Id)]                      
                        for el in self.NetworkMesh.Elements:
                            if el.NodeIds[1] == self.NetworkMesh.meshToEdges[int(node.Id)]:                                      
                                newElement = Elements.WindkesselElement(meshEndId, [meshNode1,meshNode2], node.Name, el.Side)
                                newId = newElement.Id
                                newElement.SetLastElement(el)
                                newElement.RelExpression = node.Properties['windkessel']                         
                                self.NetworkMesh.Elements.append(newElement)     
                        if self.NetworkMesh.Entities.has_key(entity):
                            self.NetworkMesh.Entities[entity].append(newElement)
                        else:
                            self.NetworkMesh.Entities[entity] = [newElement]                                      
                        self.NetworkMesh.s_mesh[(newId, 1.0)] = meshNode2
                        self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement
                        self.NetworkMesh.GraphNodeToMesh[node] = newElement.Id                 
                        meshNode2+=1
                        meshEndId+=1
     
    def ToleranceMeshing(self, Tolerance):
        '''
        This method build NetworkMesh from NetworkGraph using Radius Tolerance meshing criterium.
        If radius along the curvilinear abscissa increases more than Tolerance, a new mesh is created.
        '''
        meshId = 1 #First mesh Id index
        self.NetworkMesh.meshToEdges = {} # dict edgenode:meshnode
        self.NetworkMesh.s_mesh = {}  # dict (s,edge):meshnode
        self.NetworkMesh.MeshToGraph = {} #  {meshId : edgeId}
        self.NetworkMesh.MeshToS = {}  #  {meshId : [s]}
        self.NetworkMesh.Entities = {} # {entity:[mesh]}       
        self.NetworkMesh.GraphEdgeToMesh = {} #{edge:[meshId:[s0,s1]]}
        self.NetworkMesh.ElementIdsToElements = {} # {mesh.Id:mesh}
        self.NetworkMesh.GraphNodeToMesh = {} #{edge:[meshId:[s0,s1]]}
        nLeak = 0     
        meshNode1 = None
        meshNode2 = None        
        nodes_list = []  #nodes list
        nodes_list_last = len(nodes_list)                 
        meshNodeAnast0 = None
        meshNodeAnast1 = None
        meshNodeAnast2 = None
        done = 0  #if = 1, AVF resistances meshed
        for edgeId, edge in self.NetworkGraph.Edges.iteritems():
            if self.NetworkMesh.GraphEdgeToMesh.has_key(edge):
                pass
            else:
                edgeProperties = self.ParsingEdgeProperties(edge)
                # Searching for anastomosis element (arterial and venous side)
                for nodeId, node in self.NetworkGraph.Nodes.iteritems():
                    if node.Type == 'anastomosis':
                        if int(nodeId) in edge.NodeIds:                          
                            if edge.Side == 'venous':  
                                if done == 0:  
                                    try:
                                        meshNode0 = self.NetworkMesh.meshToEdges[edge.NodeIds[0]]
                                    except KeyError:
                                        meshNode0 = nodes_list_last+1
                                    meshNodeAnast0 = meshNode0  
                                    if meshNode0 not in nodes_list:
                                        nodes_list.append(meshNode0)
                                        nodes_list_last = len(nodes_list)   
                                    meshNode1 = nodes_list_last+1
                                    if meshNode1 not in nodes_list:
                                        nodes_list.append(meshNode1)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast1 = meshNode1
                                    meshNode2 = nodes_list_last+1
                                    if meshNode2 not in nodes_list:
                                        nodes_list.append(meshNode2)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast2 = meshNode2   
                                else:
                                    meshNode0 = meshNodeAnast0
                                    meshNode1 = meshNodeAnast1
                                    meshNode2 = meshNodeAnast2                                
                                self.NetworkMesh.meshToEdges[edge.NodeIds[0]] = meshNodeAnast2        
                            if edge.Side == 'arterial' and  int(nodeId) == edge.NodeIds[1]:
                                if done == 0:
                                    try:
                                        meshNode0 = self.NetworkMesh.meshToEdges[edge.NodeIds[1]]
                                    except KeyError:
                                        meshNode0 = nodes_list_last+1                                               
                                    meshNodeAnast0 = meshNode0
                                    if meshNode0 not in nodes_list:
                                        nodes_list.append(meshNode0)
                                        nodes_list_last = len(nodes_list)
                                    meshNode1 = nodes_list_last+1
                                    if meshNode1 not in nodes_list:
                                        nodes_list.append(meshNode1)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast1 = meshNode1
                                    meshNode2 = nodes_list_last+1
                                    if meshNode2 not in nodes_list:
                                        nodes_list.append(meshNode2)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast2 = meshNode2 
                                else:
                                    meshNode0 = meshNodeAnast0
                                    meshNode1 = meshNodeAnast1
                                    meshNode2 = meshNodeAnast2
                                self.NetworkMesh.meshToEdges[edge.NodeIds[1]] = meshNodeAnast0       
                            if edge.Side == 'arterial' and  int(nodeId) == edge.NodeIds[0]:
                                if done == 0:                                   
                                    try:
                                        meshNode1 = self.NetworkMesh.meshToEdges[edge.NodeIds[0]]
                                    except KeyError:
                                        meshNode1 = nodes_list_last+1                                               
                                    meshNodeAnast1 = meshNode1       
                                    if meshNode1 not in nodes_list:
                                        nodes_list.append(meshNode1)
                                        nodes_list_last = len(nodes_list)
                                    meshNode0 = nodes_list_last+1
                                    if meshNode0 not in nodes_list:
                                        nodes_list.append(meshNode0)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast0 = meshNode0
                                    meshNode2 = nodes_list_last+1
                                    if meshNode2 not in nodes_list:
                                        nodes_list.append(meshNode2)
                                        nodes_list_last = len(nodes_list)
                                    meshNodeAnast2 = meshNode2 
                                else:                                   
                                    meshNode0 = meshNodeAnast0
                                    meshNode1 = meshNodeAnast1
                                    meshNode2 = meshNodeAnast2                                 
                                self.NetworkMesh.meshToEdges[edge.NodeIds[0]] = meshNodeAnast1                                                               
                            #NEW ELEMENT
                            if done == 0:  
                                anast_res_0_1 = node.Properties['arterial_resistance']
                                anast_res_0_2 = node.Properties['venous_resistance']
                                elementParameters = {}
                                elementParameters['resistance_0_1'] = anast_res_0_1
                                elementParameters['resistance_0_2'] = anast_res_0_2
                                newElement = Elements.Anastomosis(str(meshId), [meshNode0,meshNode1,meshNode2], elementParameters, name = 'anastomosis')
                                self.NetworkMesh.Elements.append(newElement)
                                self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement
                                self.NetworkMesh.GraphNodeToMesh[node] = newElement.Id
                                meshId+=1
                                entity = Entity()
                                entity.SetId(node.Name)   
                                if self.NetworkMesh.Entities.has_key(entity):
                                    self.NetworkMesh.Entities[entity].append(newElement)
                                else:
                                    self.NetworkMesh.Entities[entity] = [newElement]
                                done = 1               
                #MESHNODE1########
                if meshNode1 == None:
                    meshNode1 = 1
                else:
                    try:
                        meshNode1 = self.NetworkMesh.meshToEdges[edge.NodeIds[0]]
                    except KeyError:
                        meshNode1 = nodes_list_last+1
                ##################
                elementParameters = {}             
                nameId = 1               
                self.NetworkMesh.meshToEdges[edge.NodeIds[0]] = meshNode1               
                if meshNode1 not in nodes_list:
                        nodes_list.append(meshNode1)
                        nodes_list_last = len(nodes_list)               
                i = 0
                s_list=[]              
                for s in edgeProperties['radius'].iterkeys():
                    s_list.append(s)
                    s_list.sort()     
                    
                print edge.Name, s_list        
                while i < (len(s_list)-1):
                    radius_1 = edgeProperties['radius'][s_list[i]]
                    maxVar = radius_1*self.Tolerance['radius']
                    radius_2 = edgeProperties['radius'][s_list[i+1]]                  
                    if edgeProperties['xradius'][0.0] is not None:
                        elementParameters['xradius'] = {s_list[i]:edgeProperties['xradius'][s_list[i]], s_list[i+1]:edgeProperties['xradius'][s_list[i+1]]}   
                    else:
                        elementParameters['xradius'] = None
                    if edgeProperties['yradius'][0.0] is not None:
                        elementParameters['yradius'] = {s_list[i]:edgeProperties['yradius'][s_list[i]], s_list[i+1]:edgeProperties['yradius'][s_list[i+1]]} 
                    else:
                        elementParameters['yradius'] = None
                                                
                    if abs(radius_1-radius_2) <= maxVar:
                        try:
                            elementParameters['radius'].update({s_list[i]:radius_1})
                        except KeyError:
                            elementParameters['radius']={s_list[i]:radius_1}
                        if s_list[i+1] == 1.0:
                            elementParameters['radius'].update({s_list[i+1]:radius_2})               
                        if edgeProperties['xradius'][0.0] is not None:
                            elementParameters['xradius'] = {s_list[i]:edgeProperties['xradius'][s_list[i]], s_list[i+1]:edgeProperties['xradius'][s_list[i+1]]}   
                        else:
                            elementParameters['xradius'] = None
                        if edgeProperties['yradius'][0.0] is not None:
                            elementParameters['yradius'] = {s_list[i]:edgeProperties['yradius'][s_list[i]], s_list[i+1]:edgeProperties['yradius'][s_list[i+1]]} 
                        else:
                            elementParameters['yradius'] = None    
                        for s in elementParameters['radius'].iterkeys():
                            
                            if s not in s_list:
                                s_list.append(s)
                            try:                               
                                elementParameters['wall_thickness'].update({s:edgeProperties['wallthickness'][s]})
                            except KeyError:                             
                                elementParameters['wall_thickness']={s:edgeProperties['wallthickness'][s]}
                            try: 
                                elementParameters['young_modulus'].update({s:edgeProperties['youngmodulus'][s]})
                            except KeyError:
                                elementParameters['young_modulus']={s:edgeProperties['youngmodulus'][s]}      
                        s_start = min(s_list)
                        s_end = max(s_list)
                        elementParameters['s1'] = s_start
                        elementParameters['s2'] = s_end
                        elementParameters['length'] = (s_end-s_start)*edgeProperties['length']
                        elementParameters["resistance"] = edge.Resistance
                        elementParameters["compliance"] = edge.Compliance
                        elementParameters["leakage"] = edge.QLeakage
                        if edge.QLeakage:
                            nLeak+=1                       
                        try:
                            meshNode2 = self.NetworkMesh.meshToEdges[edge.NodeIds[1]]      
                        except KeyError:
                            meshNode2 = nodes_list_last+1               
                        self.NetworkMesh.meshToEdges[edge.NodeIds[1]] = meshNode2
                        if meshNode2 not in nodes_list:
                            nodes_list.append(meshNode2)
                            nodes_list_last = len(nodes_list)                       
                        name = edge.Name + "_" + str(nameId)  
                        
                        newElement = Elements.NewElement(self.MeshType.get(edgeId,None),str(meshId), [meshNode1,meshNode2], elementParameters, edge.Side, name) 
                        
                        self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement 
                        self.NetworkMesh.Elements.append(newElement)                        
                        self.NetworkMesh.s_mesh[(s_start,edge)] = meshNode1
                        self.NetworkMesh.s_mesh[(s_end,edge)] = meshNode2                                                       
                        self.NetworkMesh.MeshToGraph[str(meshId)] = edge.Id
                        self.NetworkMesh.MeshToS[meshId] = [s_start,s_end]                       
                        if self.NetworkMesh.GraphEdgeToMesh.has_key(edge):
                            self.NetworkMesh.GraphEdgeToMesh[edge].append((meshId, self.NetworkMesh.MeshToS[meshId]))
                        else:
                            self.NetworkMesh.GraphEdgeToMesh[edge] = [(meshId, self.NetworkMesh.MeshToS[meshId])]                                            
                        meshId+=1                       
                        if meshNode2 not in nodes_list:
                            nodes_list.append(meshNode2)
                            nodes_list_last = len(nodes_list)                
                        meshNode1 = meshNode2
                        nameId+=1                               
                        i+=1                       
                    if abs(radius_1-radius_2) > maxVar:                       
                        s1 = s_list[i]
                        s2 = s_list[i+1]
                        try:
                            elementParameters['radius'].update({s1:radius_1, s2:radius_2})
                        except KeyError:
                            elementParameters['radius']={s1:radius_1, s2:radius_2}                     
                        if edgeProperties['xradius'][0.0] is not None:
                            elementParameters['xradius'] = {s1:edgeProperties['xradius'][s1], s2:edgeProperties['xradius'][s2]}   
                        else:
                            elementParameters['xradius'] = None
                        if edgeProperties['yradius'][0.0] is not None:
                            elementParameters['yradius'] = {s1:edgeProperties['yradius'][s1], s2:edgeProperties['yradius'][s2]} 
                        else:
                            elementParameters['yradius'] = None                          
                        for s in elementParameters['radius'].iterkeys():
                            if s not in s_list:
                                s_list.append(s)
                                s_list.append(s)
                                try:     
                                    elementParameters['wall_thickness'].update({s:edgeProperties['wallthickness'][s]})
                                except KeyError:                           
                                    elementParameters['wall_thickness']={s:edgeProperties['wallthickness'][s]}
                                try:    
                                    elementParameters['young_modulus'].update({s:edgeProperties['youngmodulus'][s]})
                                except KeyError: 
                                    elementParameters['young_modulus']={s:edgeProperties['youngmodulus'][s]}                      
                        s_start = min(s_list)
                        s_end = max(s_list)
                        elementParameters['s1'] = s1
                        elementParameters['s2'] = s2
                        elementParameters['length'] = (s_end-s_start)*edgeProperties['length']
                        elementParameters["resistance"] = edge.Resistance
                        elementParameters["compliance"] = edge.Compliance
                        elementParameters["leakage"] = edge.QLeakage
                        if edge.QLeakage:
                            nLeak+=1                    
                        try:
                            meshNode2 = self.NetworkMesh.meshToEdges[edge.NodeIds[1]]      
                        except KeyError:
                            meshNode2 = nodes_list_last+1                   
                        self.NetworkMesh.meshToEdges[edge.NodeIds[1]] = meshNode2
                        if meshNode2 not in nodes_list:
                            nodes_list.append(meshNode2)
                            nodes_list_last = len(nodes_list)                       
                        name = edge.Name + "_" + str(nameId) 
                               
                        newElement= newElement.NewElement(self.MeshType.get(edgeId,None), str(meshId), [meshNode1,meshNode2], elementParameters, edge.Side, name) 
                        
                        self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement 
                        self.NetworkMesh.Elements.append(newElement)                      
                        self.NetworkMesh.s_mesh[(s1,edge)] = meshNode1
                        self.NetworkMesh.s_mesh[(s2,edge)] = meshNode2                                                      
                        self.NetworkMesh.MeshToGraph[str(meshId)] = edge.Id
                        self.NetworkMesh.MeshToS[meshId] = [s1,s2]                     
                        if self.NetworkMesh.GraphEdgeToMesh.has_key(edge):
                            self.NetworkMesh.GraphEdgeToMesh[edge].append((meshId, self.NetworkMesh.MeshToS[meshId]))
                        else:
                            self.NetworkMesh.GraphEdgeToMesh[edge] = [(meshId, self.NetworkMesh.MeshToS[meshId])]                          
                        meshId+=1      
                        if meshNode2 not in nodes_list:
                            nodes_list.append(meshNode2)
                            nodes_list_last = len(nodes_list)
                        meshNode1 = meshNode2
                        nameId+=1
                    i+=1 
        #SETTING ANASTOMOSIS CONNECTIONS DICTIONARIES
        for nodeId, node in self.NetworkGraph.Nodes.iteritems():                              
            if node.Type == 'anastomosis':
                edge = self.NetworkGraph.Edges[node.Properties['proximal']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][1] == 1.0:
                            proximal = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetProximal(proximal)
                edge = self.NetworkGraph.Edges[node.Properties['distal']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][0] == 0.0:
                            distal = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetDistal(distal)
                edge = self.NetworkGraph.Edges[node.Properties['vein']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][0] == 0.0:
                            vein = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetVein(vein) 
        #SETTING LEAKAGES
        for el in self.NetworkMesh.Elements:
            if el.Type == "WavePropagation":
                el.Leakages = nLeak
        #ENTITIES DICTIONARY                
        for sedge in self.NetworkGraph.SuperEdges.values():
            entity = Entity()
            entity.SetId(sedge.Name)
            for edge in sedge.Edges.values(): 
                for key, value in self.NetworkMesh.MeshToGraph.iteritems():                    
                    if value == edge.Id:
                        for el in self.NetworkMesh.Elements:
                            if el.Type != 'Resistance':  
                                if el.Id == key:                             
                                    if self.NetworkMesh.Entities.has_key(entity):
                                        self.NetworkMesh.Entities[entity].append(el)
                                    else:
                                        self.NetworkMesh.Entities[entity] = [el]
                                        
        #SETTING ANASTOMOSIS CONNECTIONS DICTIONARIES
        for nodeId, node in self.NetworkGraph.Nodes.iteritems():                              
            if node.Type == 'anastomosis':
                edge = self.NetworkGraph.Edges[node.Properties['proximal']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][1] == 1.0:
                            proximal = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetProximal(proximal)
                edge = self.NetworkGraph.Edges[node.Properties['distal']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][0] == 0.0:
                            distal = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetDistal(distal)
                edge = self.NetworkGraph.Edges[node.Properties['vein']]
                for meshlist in self.NetworkMesh.GraphEdgeToMesh[edge]:
                    for meshId in meshlist:
                        if self.NetworkMesh.MeshToS[meshId][0] == 0.0:
                            vein = self.NetworkMesh.ElementIdsToElements[str(meshId)]
                            self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]].SetVein(vein)                               
        # MESHING END SEGMENTS (WINDKESSEL ELEMENTS) 
        entity = Entity()
        entity.SetId('end_segment')        
        nodeslist = []
        for nodeId, node in self.NetworkGraph.Nodes.iteritems():
            nodeslist.append(int(node.Id))
            nodeslist.sort()
        nodes_list_last = len(nodes_list)
        meshNode2 = nodes_list_last+1
        meshEndId = 1 
        for nd in nodeslist:
            for nodeId, node in self.NetworkGraph.Nodes.iteritems():
                if nd == int(node.Id):
                    if node.Type == 'downstream network':                     
                        meshNode1 = self.NetworkMesh.meshToEdges[int(node.Id)]                       
                        for el in self.NetworkMesh.Elements:                           
                            if el.NodeIds[1] == self.NetworkMesh.meshToEdges[int(node.Id)]:                               
                                newElement = Elements.WindkesselElement(meshEndId, [meshNode1,meshNode2], node.Name, el.Side)
                                newId = newElement.Id
                                newElement.SetLastElement(el)
                                newElement.RelExpression = node.Properties['windkessel']
                                self.NetworkMesh.Elements.append(newElement)
                                self.NetworkMesh.ElementIdsToElements[newElement.Id] = newElement
                        self.NetworkMesh.GraphNodeToMesh[node] = newElement.Id
                        if self.NetworkMesh.Entities.has_key(entity):
                            self.NetworkMesh.Entities[entity].append(newElement)
                        else:
                            self.NetworkMesh.Entities[entity] = [newElement]  
                        self.NetworkMesh.s_mesh[(newId, 1.0)] = meshNode2                
                        meshNode2+=1
                        meshEndId+=1          
                  
class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass

class MeshingError(Error):
    '''
    Exception raised for Meshing Failure
    '''
    def __init__(self):
        print "Error, Please set MaxLength or Tolerance before meshing the network."