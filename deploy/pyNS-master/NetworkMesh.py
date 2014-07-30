#!/usr/bin/env python

## Program:   PyNS
## Module:    NetworkMesh.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from numpy.core.numeric import arange
from xml.etree import ElementTree as etree

class NetworkMesh(object):
    '''
    NetworkMesh is a list of elements and entities.
    Each entities is a group of elements.
    Each element is a circuit.
    NodesToElement is a dictionary: {nodeId:[Element]}. Each node (dictionary's key) is associated with its own element(s) (list).
    This class provide the following methods:
    WriteToXML writes network mesh XML file according to its schema.
    BuildNodesToElement builds nodes'dictionary.
    CheckLinearConsistence checks and fixes the correct proportion between the meshes of each edge.
    '''
   
    def __init__(self):
        '''
        Class Constructor
        '''
        self.Id = None
        self.Elements = []
        self.NodesToElement = {}
        self.ElementIdsToElements = {} # {mesh.Id:mesh}
        self.MeshToGraph = {} #  {meshId : edgeId}
        self.MeshToS = {}  #  {meshId : [s]}
        self.s_mesh = {}  # {(s,edge):meshnode}
        self.meshToEdges = {}  # {edgenode:meshnode}
        self.Entities = {} # {entity:[mesh]}       
        self.GraphEdgeToMesh = {} #{edge:[meshId:[s0,s1]]}
        self.GraphNodeToMesh = {} #{node:meshId}
      
    def BuildNodesToElements(self):
        '''
        Building nodes' dictionary, each node has associated its element(s)
        '''
        nodesToElement = {}
        for element in self.Elements:
            for nodes in element.NodeIds:
                if nodesToElement.has_key(nodes):
                    nodesToElement[nodes].append(element.Id)
                else:
                    nodesToElement[nodes] = [element.Id]
        return nodesToElement
    
    def WriteToXML(self, xmlmeshpath):
        '''
        This method writes network_mesh XML file.
        '''
        print "Writing xml Mesh file..."
        node_list = []
        root = etree.Element("NetworkMesh", id=self.Id, version="2.0")
        xmlmesh = etree.ElementTree(root)
        meshnodes = etree.SubElement(root, "meshnodes")
        entities = etree.SubElement(root, "entities")
        elements = etree.SubElement(root, "elements")
        
        #MESHNODES
        for meshnode in self.s_mesh.itervalues():
            node_list.append(meshnode)
            node_list = list(set(node_list))
            node_list.sort()
        for n in node_list:
            etree.SubElement(meshnodes, "meshnode", id = str(n))
            
        #ENTITIES 
        for ent, el_list in self.Entities.iteritems():
            entity = etree.SubElement(entities, "entity", id = ent.Id)
            for el in el_list:       
                etree.SubElement(entity, "mesh", id = el.Id)
      
        #ELEMENTS    
        for mesh in self.Elements:
            if self.MeshToGraph.has_key(mesh.Id):
                element = etree.SubElement(elements, "element", id = str(mesh.Id), nodeIds = str(mesh.NodeIds).strip('[]'), type = mesh.Type)
                edgeId = self.MeshToGraph[mesh.Id]            
                mesh_coordinates = etree.SubElement(element, "pcoord", edgeId=edgeId)
                s1 = etree.SubElement(mesh_coordinates,"s1")
                s2 = etree.SubElement(mesh_coordinates,"s2")
                s_coord1 = self.MeshToS[int(mesh.Id)][0]
                s_coord2 = self.MeshToS[int(mesh.Id)][1]
                s1.text = str(s_coord1)
                s2.text = str(s_coord2)  
                mesh_parameters = etree.SubElement(element, "parameters")               
                if mesh.Type == 'Resistance':
                    mesh_R = etree.SubElement(mesh_parameters, "Resistance")
                    R_value = etree.SubElement(mesh_R, "scalar")
                    R_value.text = str(mesh.R)
                else:
                    mesh_length = etree.SubElement(mesh_parameters, "length", unit="m")
                    length_value = etree.SubElement(mesh_length, "scalar")
                    length_value.text = str(mesh.Length)
                    mesh_radius = etree.SubElement(mesh_parameters, "radius", unit="m")
                    radius_v1 = etree.SubElement(mesh_radius, "value", s=str(s_coord1))
                    radius_value1 = etree.SubElement(radius_v1, "scalar") 
                    if type(mesh.Radius) is dict:
                        radius_value1.text = str(mesh.Radius[s_coord1])
                    else:
                        radius_value1.text = str(mesh.Radius[0])
                    radius_v2 = etree.SubElement(mesh_radius, "value", s=str(s_coord2))
                    radius_value2 = etree.SubElement(radius_v2, "scalar")
                    if type(mesh.Radius) is dict:
                        radius_value2.text = str(mesh.Radius[s_coord2])
                    else:
                        radius_value2.text = str(mesh.Radius[len(mesh.Radius)-1])
                    wall_thickness = etree.SubElement(mesh_parameters, "wall_thickness", unit="m")
                    if type(mesh.WallThickness) is str:
                        wall_thickness_v = etree.SubElement(wall_thickness, "value")
                        wall_thickness_value = etree.SubElement(wall_thickness_v, "expression")
                        wall_thickness_value.text = str(mesh.WallThickness)
                    if type(mesh.WallThickness) is not str:
                        wall_thickness_v1 = etree.SubElement(wall_thickness, "value", s=str(s_coord1))
                        wall_thickness_value1 = etree.SubElement(wall_thickness_v1, "scalar")  
                        wall_thickness_value1.text = str(mesh.WallThickness[s_coord1])
                        wall_thickness_v2 = etree.SubElement(wall_thickness, "value", s=str(s_coord2))
                        wall_thickness_value2 = etree.SubElement(wall_thickness_v2, "scalar")
                        wall_thickness_value2.text = str(mesh.WallThickness[s_coord2]) 
                    young_modulus = etree.SubElement(mesh_parameters, "young_modulus", unit="Pa")
                    young_modulus_v1 = etree.SubElement(young_modulus, "value", s=str(s_coord1))
                    young_modulus_value1 = etree.SubElement(young_modulus_v1, "scalar")
                    young_modulus_value1.text = str(mesh.YoungModulus[s_coord1])
                    young_modulus_v2 = etree.SubElement(young_modulus, "value", s=str(s_coord2))
                    young_modulus_value2 = etree.SubElement(young_modulus_v2, "scalar")
                    young_modulus_value2.text = str(mesh.YoungModulus[s_coord2])              
                    mesh_R = etree.SubElement(mesh_parameters, "Resistance")
                    R_value = etree.SubElement(mesh_R, "scalar")
                    R_value.text = str(mesh.R)                   
                    mesh_C = etree.SubElement(mesh_parameters, "Compliance")
                    C_value = etree.SubElement(mesh_C, "scalar")
                    C_value.text = str(mesh.C)
        for mesh in self.Elements:       
            if self.MeshToGraph.has_key(mesh.Id) == False:               
                if mesh.Type == 'Anastomosis':
                    element = etree.SubElement(elements, "element", id = str(mesh.Id), nodeIds = str(mesh.NodeIds).strip('[]'), type = mesh.Type)
                    for node in self.GraphNodeToMesh:
                        if node.Type == 'anastomosis':
                            nodeid = node.Id
                    mesh_coordinates = etree.SubElement(element, "pcoord", nodeId=nodeid)
                    mesh_parameters = etree.SubElement(element, "parameters")
                    mesh_Resistance01 = etree.SubElement(mesh_parameters, "Resistance_0_1")
                    Resistance01_value = etree.SubElement(mesh_Resistance01, "scalar")
                    Resistance01_value.text = str(mesh.R_0_1)
                    mesh_Resistance02 = etree.SubElement(mesh_parameters, "Resistance_0_2")
                    Resistance_value02 = etree.SubElement(mesh_Resistance02, "scalar")
                    Resistance_value02.text = str(mesh.R_0_2)
                else:
                    element = etree.SubElement(elements, "element", id = str(mesh.Id), nodeIds = str(mesh.NodeIds).strip('[]'), type = mesh.Type)
                    for node in self.GraphNodeToMesh:
                        if node.Type == 'downstream network':
                            if self.GraphNodeToMesh[node] == mesh.Id:
                                nodeid = node.Id
                                mesh_coordinates = etree.SubElement(element, "pcoord", nodeId=nodeid)
                                mesh_parameters = etree.SubElement(element, "parameters")
                                mesh_R1 = etree.SubElement(mesh_parameters, "Wave_Impedance")
                                R1_value = etree.SubElement(mesh_R1, "scalar")
                                R1_value.text = str(mesh.R1)
                                mesh_R2 = etree.SubElement(mesh_parameters, "Peripheral_Resistance")
                                R2_value = etree.SubElement(mesh_R2, "scalar")
                                R2_value.text = str(mesh.R2)     
                                mesh_C = etree.SubElement(mesh_parameters, "Compliance")
                                C_value = etree.SubElement(mesh_C, "scalar")
                                C_value.text = str(mesh.C)
        indent(root)                
        xmlmesh.write (xmlmeshpath, encoding='iso-8859-1')       

    def checkLinearConsistence(self):
        '''
        This method checks and fixes the correct proportion between the meshes of each edge.
        '''
        for edge in self.GraphEdgeToMesh.iterkeys():
            if edge.Side =='venous':
                meshes = self.GraphEdgeToMesh[edge]      
                startingMesh = meshes[0]
                startingRadius = self.ElementIdsToElements[str(startingMesh.keys()[0])].Radius[0]
                endingMesh = meshes[len(meshes)-1]
                endingRadius = self.ElementIdsToElements[str(endingMesh.keys()[0])].Radius[len(self.ElementIdsToElements[str(endingMesh.keys()[0])].Radius)-1]
                elLen = edge.Length['value']/len(meshes)
                self.dz = elLen/1.0e5
                z = arange(0.0,elLen,self.dz)
                dr = (endingRadius-startingRadius)/(len(meshes))
                for mesh in meshes:
                    if mesh == startingMesh:
                        r1 = startingRadius
                        r2 = r1+dr
                    elif mesh == endingMesh:     
                        r1 = r2
                        r2 = endingRadius
                    else: 
                        r1 = r2
                        r2+=dr
                    r_z = r1+((dr/elLen)*z)
                    self.ElementIdsToElements[str(mesh.keys()[0])].Radius = r_z
                    
         
class Entity(object):
    '''
    Each entities is a group of elements, 
    identified by an ID (according to superedges).
    '''
   
    def __init__(self):
        '''
        Class Constructor
        '''
        self.Id = None
        self.Leakage = False
        
    def SetId(self, id):
        '''
        This method sets unique Id
        '''
        self.Id = id
        
    def SetLeakage(self):
        '''
        This method sets Leakage True or False
        '''
        self.Leakage = True    
      
class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass

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