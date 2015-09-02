#!/usr/bin/env python

## Program:   PyNS
## Module:    NetworkGraph.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390


from numpy.core.numeric import array
from numpy.core.fromnumeric import mean
from math import pi
import csv, sys

class NetworkGraph(object):
    '''
    NetworkGraph represents the vascular network, composed by Nodes, SuperEdges and Edges.
    Nodes, SuperEdges and Edges are objects themselves.
    This class provides the following methods:
    SetNodes, SetSuperEdges, SetEdges, SetId for setting elements' class.
    GetEdge: a method for returning the corresponding edge from provided nodes.
    GetNodeEdges: a method for returning all edges connected to provided node.
    GetSuperEdge: a method for returning the corresponding superedge(s) from provided edge.
    ReadFromXML: a method for reading XML NetworkGraph file.
    WriteCsv: a method for writing network graph properties (length, radius and compliance) in a .csv file.
    '''
   
    def __init__(self):
        '''
        Class Constructor
        '''
        self.xmlgraphpath = None
        self.Id = None
        self.PatientId = None
        self.Visit = None
        self.Nodes = {} #node_id:Node
        self.NodeNamesToIds = {} # {nodename:nodeId}
        self.SuperEdges = {} #superEdge_id:superEdge
        self.SuperEdgeNamesToIds = {} # {superedgename:superedgeId}
        self.Edges = {} #edge_id:Edge
        self.EdgeNamesToIds = {} # {edgename:edgeId}
        
    def SetNodes(self, node):
        '''
        This method sets Nodes dictionary.
        Node_Id:Node
        '''
        self.Nodes[node.Id] = node
        if node.Name is not None:
            self.NodeNamesToIds[node.Name] = node.Id
        
    def SetSuperEdges(self, superEdge):
        '''
        This method sets SuperEdges dictionary.
        SuperEdge_Id:SuperEdge
        '''
        self.SuperEdges[superEdge.Id] = superEdge
        self.SuperEdgeNamesToIds[superEdge.Name] = superEdge.Id
        
    def SetEdges(self, edge):
        '''
        This method sets Edges dictionary.
        Edge_Id:Edge
        '''
        self.Edges[edge.Id] = edge
        self.EdgeNamesToIds[edge.Name] = edge.Id
        
    def SetId(self, id):
        '''
        This method sets unique ID
        '''
        self.Id = id
        
    def GetEdge(self, nodes):
        '''
        This method returns the corresponding edge from provided nodes
        '''
        nodes.sort()
        for edgeId, edge in self.Edges.iteritems():
            edge.NodeIds.sort()
            if edge.NodeIds == nodes:
                return edgeId
               
    def GetNodeEdges(self, node):
        '''
        This method returns all edges connected to provided node.
        '''
        nodeEdges = []
        for edge in self.Edges.itervalues():
            if node in edge.NodeIds:
                nodeEdges.append(edge)
        return nodeEdges
                
    def GetSuperEdge(self, edge):
        '''
        This method returns the corresponding superedge(s) from provided edge
        '''
        superEdge = []
        for superedge in self.SuperEdges.itervalues():
            if superedge.Edges.has_key(edge):
                superEdge.append(superedge)
        return superEdge
        
    def ReadFromXML(self, xmlgraphpath, xsdgraphpath=None):
        '''
        This method reads Network Graph XML File.
        If XML schema is given and lxml package is installed,
        XML file is validated first.
        '''
        self.xmlgraphpath = xmlgraphpath
        
        try:
            from lxml import etree
            lxml = True
        except:
            LXMLError()
            lxml = False
            from xml.etree import ElementTree as etree
        
        if lxml:
            if not xsdgraphpath:
                NoXSDWarning()
            while True:
                try:
                    schemagraphfile = open(xsdgraphpath)
                except:
                    WrongXSDPathError()
                    break
                try:
                    xmlschema_doc = etree.parse(schemagraphfile)
                    xmlschema = etree.XMLSchema(xmlschema_doc)
                    docgraphfile = open(xmlgraphpath)
                    docgraph = etree.parse(docgraphfile)
                    xmlschema.assert_(docgraph)
                    print "Network Graph Xml File has been validated."
                    break
                except AssertionError:   
                    XMLValidationError(xmlschema)
                            
        docgraphfile = open(xmlgraphpath)
        graphtree=etree.parse(docgraphfile)
        graph=graphtree.getroot()
        graph_dict = graph.attrib
        self.SetId(graph_dict['id'])  # Setting Graph ID
        #CASE
        for case in graph.findall(".//case"):
            for patid in case.findall(".//patient_id"):
                self.PatientId = patid.text
            for visit in case.findall(".//visit"):
                self.Visit = visit.text
        # NODES
        for nodeg in graph.findall(".//node"):
            node_dict = nodeg.attrib
            node = Node()  # New Node
            node.SetId(node_dict['id'])  # Setting Node ID
            try:
                node.SetType(node_dict['type'])  # Setting Node Type
                node.SetName(node_dict['name'])  # Setting Node Name
            except KeyError:
                pass          
            for prop in nodeg.findall(".//properties"):
                for data in prop:              
                    resistance = {}
                    compliance = {}
                    connections = {}                                 
                    if data.tag == 'arterial_resistance':
                        for a_resistance in data.findall(".//expression"):
                            resistance['arterial']=a_resistance.text
                            node.SetProperties(resistance)
                        for a_resistance in data.findall(".//scalar"):
                            resistance['arterial']=float(a_resistance.text)                    
                            node.SetProperties(resistance)                    
                    if data.tag == 'venous_resistance':
                        for v_resistance in data.findall(".//expression"):
                            resistance['venous']=v_resistance.text
                            node.SetProperties(resistance)
                        for v_resistance in data.findall(".//scalar"):
                            resistance['venous']=float(v_resistance.text)
                            node.SetProperties(resistance)
                    if data.tag == 'connections':
                        for prox in data.findall(".//proximal_artery"):
                            prox_dict = prox.attrib
                            connections['proximal'] = prox_dict['edge_id']
                            node.SetProperties(connections)
                        for prox in data.findall(".//distal_artery"):
                            prox_dict = prox.attrib
                            connections['distal'] = prox_dict['edge_id']
                            node.SetProperties(connections)
                        for prox in data.findall(".//proximal_vein"):
                            prox_dict = prox.attrib
                            connections['vein'] = prox_dict['edge_id']
                            node.SetProperties(connections)                      
                    if data.tag == 'windkessel':
                        for windk in data.findall(".//expression"):
                            resistance['windkessel']=windk.text
                            node.SetProperties(resistance)                                                      
            self.SetNodes(node)  # Adding New Node in Network Graph            
        # EDGES  
        for edgeg in graph.findall(".//edge"):
            edge_dict = edgeg.attrib     
            edge = Edge()  #New Edge
            edge.SetId(edge_dict['id'])  #Setting Edge Id
            edge.SetNodes([int(edge_dict['node1_id']), int(edge_dict['node2_id'])])  # Setting Edge Nodes
            edge.SetSide(edge_dict['side'])  # Setting Edge Side
            edge.SetName(edge_dict['name'])  # Setting Edge Name
            try:
                edge.SetType(edge_dict['type'])  # Setting Edge Type        
            except:
                pass
            for geometry in edgeg.findall(".//geometry"):
                for data in geometry:
                    if data.tag == "length":                       
                        length_value = []                          
                        for length in data.findall(".//vector"):
                            length_value.append(float(length.text))   
                            edge.SetLength(length_value)  # Setting Edge Length (vector type)      
                        for length in data.findall(".//scalar"):
                            length_value.append(float(length.text)) 
                            edge.SetLength(length_value)  # Setting Edge Length (scalar type)      
                        for length in data.findall(".//expression"):
                            length_value.append(length.text)
                            edge.SetLength(length_value)  # Setting Edge Length(expression type)  
                                        
                    if data.tag == "coordinates_array":                       
                        coord_value = {}                        
                        for coord in geometry.findall(".//coordinates"):
                            coord_dict = coord.attrib
                            for point in coord:    
                                coord_value[coord_dict['s']].append(point.text)                               
                            edge.SetCoordinates(coord_value)  # Setting Edge Coordinates                        
            for properties in edgeg.findall(".//properties"):
                for data in properties:                            
                    if data.tag == "radius_array":                        
                        radius_value = {}
                        radius_dict = {}                       
                        for curv in data.findall(".//value"):
                            curv_dict = curv.attrib
                            s = float(curv_dict['s'])
                            for radius in curv.findall(".//scalar"):
                                radius_v = float(radius.text)                              
                                radius_dict[s] = radius_v
                                edge.SetScalarRadius(s, radius_v)
                            for radius in curv.findall(".//expression"):
                                radius_v = radius.text
                                radius_dict[s] = radius_v   # Setting Edge Radius(expression type)                                                   
                        radius_value['array'] = radius_dict              
                        edge.SetRadius(radius_value)  # Setting Edge Radius (array type) 
                        
                    if data.tag == "radius":                        
                        radius_value = {}                        
                        for radius in data.findall(".//scalar"):
                            radius_value['value'] = float(radius.text)                            
                            edge.SetRadius(radius_value)  # Setting Edge Radius (scalar type)
                        for radius in data.findall(".//expression"):
                            radius_value['expression'] = radius.text
                            edge.SetRadius(radius_value)  # Setting Edge Radius(expression type)
                                                                         
                    if data.tag == "radius_a":                       
                        radius_valueA = {}                                       
                        for radiusA in data.findall(".//scalar"):
                            radius_valueA['value'] = float(radiusA.text)                       
                            edge.SetRadiusxAxis(radius_valueA)  # Setting Edge Radius X axis (scalar type)
                        for radiusA in data.findall(".//expression"):
                            radius_valueA['expression'] = radiusA.text
                            edge.SetRadiusxAxis(radius_valueA)  # Setting Edge Radius X axis (expression type) 
                                                                  
                    if data.tag == "radius_b":                        
                        radius_valueB = {}  
                        radius_valueAB = {}
                        for radiusB in data.findall(".//expression"):
                            radius_valueB['expression'] = radiusB.text
                            edge.SetRadiusyAxis(radius_valueB)  # Setting Edge Radius Y axis (expression type)
                                                      
                        for radiusB in data.findall(".//scalar"):
                            radius_valueB['value'] = float(radiusB.text)                             
                            edge.SetRadiusyAxis(radius_valueB)  # Setting Edge Radius Y axis (scalar type)                            
                        radius_valueAB['value'] = (radius_valueA['value']*radius_valueB['value'])**0.5                                                                  
                        edge.SetRadius(radius_valueAB)  # Setting Edge Radius (scalar type)
                                               
                    if data.tag == "radius_a_array":
                        radius_valueA = {}
                        radius_dict = {}                        
                        for curv in data.findall(".//value"):
                            curv_dict = curv.attrib
                            s = float(curv_dict['s'])
                            for radius in curv.findall(".//scalar"):
                                radius_v_A = float(radius.text) 
                                radius_dict[s] = radius_v_A                               
                        radius_valueA['array'] = radius_dict                        
                        edge.SetRadiusxAxis(radius_valueA)  # Setting Edge Radius X axis (array type)                      
                    if data.tag == "radius_b_array":                        
                        radius_valueB = {}
                        radius_dict = {}                        
                        for curv in data.findall(".//value"):
                            curv_dict = curv.attrib
                            s = float(curv_dict['s'])
                            for radius in curv.findall(".//scalar"):
                                radius_v_B = float(radius.text)                              
                                radius_dict[s] = radius_v_B                                
                        radius_valueB['array'] = radius_dict                       
                        edge.SetRadiusyAxis(radius_valueB)  # Setting Edge Radius X axis (array type) 
                                                
                    if data.tag == "distensibility":                       
                        distensibility_value = {}                        
                        for distensibility in data.findall(".//scalar"):
                            distensibility_value['value'] = float(distensibility.text)                            
                            edge.SetDistensibility(distensibility_value)  # Setting Edge Distensibility (scalar type)                           
                    if data.tag == "distensibility_array":   
                        distensibility_value = {}
                        distensibility_dict = {}  
                        for curv in data.findall(".//value"):
                            curv_dict = curv.attrib
                            s = float(curv_dict['s'])
                            for distensibility in curv.findall(".//scalar"):
                                distensibility_v = float(distensibility.text)
                                distensibility_dict[s] = distensibility_v
                        distensibility_value['array'] = distensibility_dict
                        edge.SetDistensibility(distensibility_value)  # Setting Edge Distensibility (array type)
                    if data.tag == "wall_thickness":
                        wall_thickness_value = {}  
                        for wall_thickness in data.findall(".//scalar"):
                            wall_thickness_value['value'] = float(wall_thickness.text)
                            edge.SetWallThickness(wall_thickness_value)  # Setting Edge Wall Thickness (scalar type)                            
                        for wall_thickness in data.findall(".//expression"):                            
                            wall_thickness_value['expression'] = wall_thickness.text                            
                            edge.SetWallThickness(wall_thickness_value)  # Setting Edge Wall Thickness (scalar type)     
                    if data.tag == "wall_thickness_array":                        
                        wall_thickness_value = {}
                        wall_thickness_dict = {}                        
                        for curv in data.findall(".//value"):
                            curv_dict = curv.attrib
                            s = float(curv_dict['s'])
                            for wall_thickness in curv.findall(".//scalar"):
                                wall_thickness_v = float(wall_thickness.text)                               
                                wall_thickness_dict[s] = wall_thickness_v                                
                        wall_thickness_value['array'] = wall_thickness_dict
                        edge.SetWallThickness(wall_thickness_value)  # Setting Edge Wall Thickness (array type)                                           
                    if data.tag == "young_modulus":                        
                        young_modulus_value = {}                        
                        for young_modulus in data.findall(".//scalar"):
                            young_modulus_value['value'] = float(young_modulus.text)                           
                            edge.SetYoungModulus(young_modulus_value)  # Setting Edge Young Modulus (scalar type)
                        for young_modulus in data.findall(".//expression"):
                            young_modulus_value['expression'] = young_modulus.text
                            edge.SetYoungModulus(young_modulus_value)  # Setting Edge Young Modulus(expression type)                            
                    if data.tag == "young_modulus_array":                        
                        young_modulus_value = {}
                        young_modulus_dict = {}                        
                        for curv in data.findall(".//value"):
                            curv_dict = curv.attrib
                            s = float(curv_dict['s'])
                            for young_modulus in curv.findall(".//scalar"):
                                young_modulus_v = float(young_modulus.text)                               
                                young_modulus_dict[s] = young_modulus_v                                
                        young_modulus_value['array'] = young_modulus_dict
                        edge.SetYoungModulus(young_modulus_value)  # Setting Edge Young Modulus (array type)
                    
                    
                    if data.tag == "leakage":                        
                        for leakage in data.findall(".//expression"):
                            edge.SetQLeakage(leakage.text)                                 
                    if data.tag == "resistance":                        
                        for resistance in data.findall(".//expression"):                       
                            edge.SetResistance(resistance.text)     
                        for resistance in data.findall(".//scalar"):
                            edge.SetResistance(float(resistance.text) )                                                  
                    if data.tag == "compliance":                        
                        for compliance in data.findall(".//expression"):                       
                            edge.SetCompliance(compliance.text)     
                        for compliance in data.findall(".//scalar"):
                            edge.SetCompliance(float(compliance.text))
                    if data.tag == "nl_compliance":                        
                        for nlcompliance in data.findall(".//expression"):                       
                            edge.SetNlCompliance(nlcompliance.text)
                                                
            for features in edgeg.findall(".//features"):
                for data in features:                          
                    if data.tag == "stenosis":
                        resistance_ste = None
                        compliance_ste = None
                        data_dict = data.attrib
                        s_ste = float(data_dict['s'])
                        for radius in data.findall(".//radius"):
                            for rad in radius.findall(".//scalar"):
                                rad_ste = float(rad.text)
                        for length in data.findall(".//length"):
                            for len in length.findall(".//scalar"):
                                len_ste = float(len.text)
                        for resistance in data.findall(".//resistance"):
                            for resist in data.findall(".//expression"):
                                resistance_ste = resist.text 
                        for compliance in data.findall(".//compliance"):
                            for compl in data.findall(".//expression"):
                                compliance_ste = compl.text                        
                        edge.SetStenosis(s_ste, rad_ste, len_ste, resistance_ste, compliance_ste)                        
                    if data.tag == "kink":
                        resistance_k = None
                        compliance_k = None
                        data_dict = data.attrib
                        s_k = float(data_dict['s'])
                        for radius in data.findall(".//radius"):
                            for rad in radius.findall(".//scalar"):
                                rad_k = float(rad.text)
                        for curvature in data.findall(".//curvature"):
                            for curv in curvature.findall(".//scalar"):
                                curv_k = float(curv.text)
                        for resistance in data.findall(".//resistance"):
                            for resist in data.findall(".//expression"):
                                resistance_k = resist.text 
                        for compliance in data.findall(".//compliance"):
                            for compl in data.findall(".//expression"):
                                compliance_k = compl.text                       
                        edge.SetKink(s_k, rad_k, curv_k, resistance_k, compliance_k)                                                               
            self.SetEdges(edge) # Setting New Edge in Network Graph        
                                   
        # SUPEREDGES
        for superedgeg in graph.findall(".//superedge"):
            superedge_dict = superedgeg.attrib
            superEdge = SuperEdge() # New SuperEdge
            superEdge.SetId(superedge_dict['id']) # Setting SuperEdge Id
            superEdge.SetName(superedge_dict['name']) # Setting SuperEdge Name
            self.SetSuperEdges(superEdge)  # Setting New SuperEdge in Network Graph
            for children in superedgeg.findall('.//superedges'):
                for superedgeg_child in children:
                    superedgeg_child_dict = superedgeg_child.attrib               
                    superEdge_child = SuperEdge() # New SuperEdge
                    superEdge_child.SetId(superedgeg_child_dict['id']) # Setting SuperEdge Id
                    superEdge_child.SetName(superedgeg_child_dict['name']) # Setting SuperEdge Name
                    superEdge.SetSuperEdges(superEdge_child) # Setting SuperEdges dictionary
                    self.SetSuperEdges(superEdge_child)  # Setting New SuperEdge in Network Graph  
            for edgeg in superedgeg.findall(".//edgeIds"):
                edgeg_dict = edgeg.attrib           
                if self.Edges.has_key(edgeg_dict['edge_id']):
                    superEdge.SetEdges(self.Edges[edgeg_dict['edge_id']])  # Setting Edges dictionary
                    
    def WriteCsv(self):
        '''
        This method writes network graph arterial properties
        (length, radius and compliance) in a csv file.
        '''
        
        edges_list = []
        for edge in self.Edges.iterkeys():
            edges_list.append(int(edge))
        edges_list.sort()
        
        path = self.xmlgraphpath+'.csv' 
        ofile  = open(path, "wb")
        csv_writer = csv.writer(ofile, delimiter=",", quoting=csv.QUOTE_ALL)
        
        for edg in edges_list:    
            for e in self.Edges.itervalues():
                if e.xRadius or e.yRadius:
                    ellipticGeometry = True
                else:
                    ellipticGeometry = False
        
        
        if ellipticGeometry == True:
            csv_writer.writerow(["Name","Side", "Length", "Radius s=0", "Radius s=1","xRadius s=0", "xRadius s=1","yRadius s=0", "yRadius s=1", "Compliance", "YoungModulus"])
            csv_writer.writerow(["","", "cm", "mm", "mm","mm", "mm","mm", "mm", "mm2/kPa", "Pa"])
        if ellipticGeometry == False:
            csv_writer.writerow(["Name","Side", "Length", "Radius s=0", "Radius s=1", "Compliance", "YoungModulus"])
            csv_writer.writerow(["","", "cm", "mm", "mm", "mm2/kPa", "Pa"])
            
        for edg in edges_list:    
            for e in self.Edges.itervalues():
                if e.Id == str(edg):
                    try:
                        if 'value' in e.Radius:
                            e.Radius_0 = e.Radius['value']
                            e.Radius_1 = e.Radius['value']
                        else:
                            if type(e.Radius['array'][0.0]) is str:
                                e.Radius_0 = e.ScalarRadius[0.0]
                            else:
                                e.Radius_0 = e.Radius['array'][0.0]
                            if type(e.Radius['array'][1.0]) is str:
                                e.Radius_1 = e.ScalarRadius[1.0]
                            else:
                                e.Radius_1 = e.Radius['array'][1.0]   
                        e.xRadius_0 = e.yRadius_0 = e.xRadius_1 = e.yRadius_1 = 0.0
                    except KeyError:
                        if 'value' in e.xRadius:
                            e.xRadius_0 = e.xRadius['value']
                            e.xRadius_1 = e.xRadius['value']
                        else:
                            try:
                                e.xRadius_0 = e.xRadius['array'][0.0]
                                e.xRadius_1 = e.xRadius['array'][1.0]
                            except:
                                e.xRadius_0 = 0
                                e.xRadius_1 = 0
                        if 'value' in e.yRadius:
                            e.yRadius_0 = e.yRadius['value']
                            e.yRadius_1 = e.yRadius['value']
                        else:
                            try:
                                e.yRadius_0 = e.yRadius['array'][0.0]
                                e.yRadius_1 = e.yRadius['array'][1.0]
                            except:
                                e.yRadius_0 = 0
                                e.yRadius_1 = 0
                        e.Radius_0 = e.Radius_1 = 0.0
                        
                    if e.Compliance is not None:
                        C = e.Compliance*1e9
                    else:
                        C = ''
                    if 'value' in e.YoungModulus:
                        ym = e.YoungModulus['value']
                        rm = ((e.Radius_0+e.Radius_1)/2)*1e3
                        wt = rm * 0.2
                        C = (((2.0*pi*rm**2)*(((2.0*rm**2*(1.0-3.e-3**2))/(wt**2))+((1.0+3.e-3)*(((2.0*rm)/wt)+1.0))))/(ym*(((2.0*rm)/wt)+1.0)))*1e3
                    else:
                        ym = ''
                    if ellipticGeometry == True:
                        csv_writer.writerow([e.Name, e.Side, e.Length['value']*1e2, e.Radius_0*1e3, e.Radius_1*1e3,e.xRadius_0*1e3, e.xRadius_1*1e3,e.yRadius_0*1e3, e.yRadius_1*1e3, C, ym])
                    if ellipticGeometry == False:
                        csv_writer.writerow([e.Name, e.Side, e.Length['value']*1e2, e.Radius_0*1e3, e.Radius_1*1e3, C, ym])
        csv_writer.writerow([])
        csv_writer.writerow([])
        csv_writer.writerow(["idpat", "gender", "age", "arm", "fistula type", "height", "weight", "bsa", "pressure", "cardiac output", "cardiac frequency", "brachial flow", "radial flow", "ulnar flow", "hematocrit", "plasma concentration","dynamic_viscosity", "blood_density","hypertension", "diabetes"])
        csv_writer.writerow(["", "", "" , "", "", "cm", "kg", "m2", "mmHg", "mL/min", "Hz", "mL/min", "mL/min", "mL/min", "%", "g/dL", "cP", "Kg/m3", "", ""])
        
      
class Node(object):
    '''
    Node is a component of the vascular network.
    Each node is identified by an unique id. If useful, we can classify a node,
    defining a type and a content with a specific name referred to its location.
    If node is classified (inside XML file) as capillary or downstream network, we must provide additional
    information about its properties, such as windkessel parameters (Wave Impedance, Peripheral Resistance and Compliance).
    This class has the following methods:
    SetId, SetType, SetName and SetProperties for setting each node.
    '''
    
    def __init__(self):
        '''
        Class Constructor
        '''
        self.Id = None
        self.Type = None
        self.Name = None
        self.Properties = {}
        
    def SetId(self, id):
        '''
        This method sets unique ID
        '''
        self.Id = id
    
    def SetType(self, type):
        '''
        This method sets Type
        '''
        self.Type = type
        
    def SetName(self, name):
        '''
        This method sets Type
        '''
        self.Name = name
    
    def SetProperties(self, property = None):
        '''
        This method sets Properties (if needed)
        '''
        try:
            self.Properties['arterial_resistance'] = property['arterial']
        except KeyError:
            try:
                self.Properties['venous_resistance'] = property['venous']
            except KeyError:
                pass
        try:
            self.Properties['venous_resistance'] = property['venous']
        except KeyError:
            try:
                self.Properties['arterial_resistance'] = property['arterial']
            except KeyError:
                pass  
        try:
            self.Properties['proximal'] = property['proximal']
        except KeyError:
            pass
        try:
            self.Properties['distal'] = property['distal']
        except KeyError:
            pass
        try:
            self.Properties['vein'] = property['vein']
        except KeyError:
            pass     
        try:
            self.Properties['windkessel'] = property['windkessel']
        except KeyError:
            pass
               
class SuperEdge(object):
    '''
    SuperEdge is a component of the vascular network.
    Each of these is identified by an unique id and
    we can classify a SuperEdge defining a name.
    Superedge is a sequence of edges (superedge is formed by one or more edges) 
    or can be a sequence of others superedges
    This class has the following methods:
    SetId and SetName for setting each superedge.
    SetSuperEdges builds dictionary of superedges contained in each superedge.
    SetEdges builds dictionary of edges contained in each superedge.
    '''
    
    def __init__(self):
        '''
        Class Constructor
        '''
        self.Id = None
        self.Name = None
        self.SuperEdges = {} # Dictionary of superedges contained in current superedge. superEdge_id:superEdge 
        self.Edges = {} # Dictionary of Edges contained in current superedge. edge_id:edge
     
    def SetId(self, id):
        '''
        This method sets unique Id
        '''
        self.Id = id 
        
    def SetName(self, name):
        '''
        This method sets name
        '''
        self.Name = name
        
    def SetSuperEdges(self, superEdge):
        '''
        This method sets superedge dictionary, sequence of superedges inside a superedge.
        '''
        self.SuperEdges[superEdge.Id] = superEdge
        
    def SetEdges(self, edge):
        '''
        This method sets edges dictionary, sequence of edges inside a superedge
        '''
        self.Edges[edge.Id] = edge
        
    def GetRadius(self,info):
        '''
        This method returns superedge's radius
        '''
        pass
                
class Edge(object):
    '''
    Edge represents a vessel of vascular network layout.
    Edge is marked by n nodes, name, id, side.
    Edge has geometries and properties parameters.
    This class has the following methods:
    SetId, SetName, SetNodes, SetSide, SetLength, SetCoordinates, SetTrasformationId, SetRadius, SetRadiusxAxis,
    SetRadiusyAxis, SetDistensibility, SetWallThickness and SetYoungModulus for setting each edge's property.
    SetQLeakage for setting leakage resistance expression.
    SetResistance for setting non linear resistance.
    SetCompliance and SetNlCompliance for setting linear and non linear compliance.
    SetStenosis for setting stenosis feature.
    SetKink for setting kink feature.
    GetRadius, GetLength and GetYoungModulus for returning edge' radius, length or Young modulus respectively.
    '''
    
    def __init__(self):
        '''
        Class Constructor
        '''
        self.NodeIds = []
        self.Name = None
        self.Id = None
        self.Side = None
        self.Type = None
        self.TransformationId = None
        self.Length = {}
        self.Coordinates = {}
        self.Radius = {}
        self.ScalarRadius = {}
        self.xRadius = {}
        self.yRadius = {}
        self.Distensibility = {}
        self.WallThickness = {}
        self.YoungModulus = {}
        self.QLeakage = None
        self.Resistance = None
        self.Compliance = None
        self.NlCompliance = {}
        self.Stenosis = None
        self.Kink = None
        self.edgeAbscissa = None
        
    def SetNodes(self, NodeIds):
        '''
        This method sets nodes
        '''
        self.NodeIds[:] = NodeIds[:]
        
    def SetName(self, name):
        '''
        This method sets name 
        '''
        self.Name = name
        
    def SetId(self, id):
        '''
        This method sets Id
        '''
        self.Id = id
        
    def SetSide(self, side):
        '''
        This method sets side (arterial or venous)
        '''
        self.Side = side
    
    def SetType(self, type):
        '''
        This method sets type. Elements belonging to 
        this edge will be meshed according to this specific type.
        '''
        self.Type = type
        
    def SetTransformationId(self, transformationId):
        '''
        This method sets transformation Id
        '''
        self.TransformationId = transformationId
    
    def SetLength(self, length):
        '''
        This method sets length.
        Length can be a single value or an array of values (x, y, z)
        '''
        if type(length) == float:
            self.Length['value'] = length
        else:    
            for x in length:
                if type(x) == str:
                    self.Length['expression'] = length[0]
            if len(length) == 1:
                self.Length['value'] = length[0]           
            if len(length) == 3:
                self.Length['array'] = length
        
    def SetCoordinates(self, coordinates):
        '''
        This method sets coordinates (s:[x, y ,z])
        '''
        self.Coordinates.update(coordinates)
        
    def SetRadius(self, radius):
        '''
        This method sets Radius
        Radius can be a single value or an array of values.
        s:value        
        '''
        if type(radius) is not dict:
            if self.edgeAbscissa is not None:
                self.Radius['array'][self.edgeAbscissa] = radius    
            else:
                self.Radius.update({'value':radius})       
        else:
            self.Radius.update(radius)
            
    def SetScalarRadius(self, s, radius):
        '''
        This method sets Radius
        Radius can be a single value or an array of values.
        s:value    
        '''
        self.ScalarRadius[s] = radius
        
        
    def SetRadiusxAxis(self, radius):
        '''
        This method set X Axis Radius (Venous has an elliptic geometry)
        Radius can be a single value or an array of values.
        s:value
        '''
        self.xRadius.update(radius)
    
    def SetRadiusyAxis(self, radius):
        '''
        This method set Y Axis Radius (Venous has an elliptic geometry)
        Radius can be a single value or an array of values.
        s:value
        '''
        self.yRadius.update(radius)
        
    def SetDistensibility(self, distensibility):
        '''
        This method sets Distensibility
        Distensibility can be a single value or an array of values.
        s:value
        '''
        self.Distensibility.update(distensibility)
        
    def SetWallThickness(self, wallThickness):
        '''
        This method sets Wall Thickness
        Wall Thickness can be a single value or an array of values.
        s:value
        '''
        self.WallThickness.update(wallThickness)
    
    def SetYoungModulus(self, youngModulus):
        '''
        This method sets Young Modulus
        Young Modulus can be a single value or an array of values.
        s:value
        '''
        if type(youngModulus) == float:
            self.YoungModulus['value']=youngModulus
        else:
            self.YoungModulus.update(youngModulus) 
    
    def SetQLeakage(self, qleakage):
        '''
        This method sets q_leakage expression.
        '''
        self.QLeakage = qleakage
    
    def SetResistance(self, resistance):
        '''
        This method sets non linear resistance.
        '''
        self.Resistance = resistance
    
    def SetCompliance(self, compliance):
        '''
        This method sets compliance.
        '''
        if type(compliance) is str:
            self.Compliance = {}
            self.Compliance['expression'] = compliance
        else:
            self.Compliance = compliance
            
    def SetNlCompliance(self, compliance):
        '''
        This method sets non linear compliance.
        '''
        self.NlCompliance['expression'] = compliance
        
    def SetStenosis(self, s, radius, length, resistance, compliance):
        '''
        This method sets stenosis feature.
        '''
        Stenosis = {}       
        Stenosis[s]= [radius]
        Stenosis[s].append(length)
        Stenosis[s].append(resistance)
        Stenosis[s].append(compliance)        
        self.Stenosis = Stenosis
    
    def SetKink(self, s, radius, curvature, resistance, compliance):
        '''
        This method sets kink feature.
        '''
        Kink = {}        
        Kink[s]= [radius]
        Kink[s].append(curvature)
        Kink[s].append(resistance)
        Kink[s].append(compliance)       
        self.Kink = Kink
    
    def GetRadius(self, abscissa):
        '''
        This method returns edge's radius
        '''
        if 'value' in self.Radius:
            return self.Radius['value']
        if 'array' in self.Radius:
            if abscissa is not None:
                return self.Radius['array'][abscissa]
            else:
                if self.edgeAbscissa is None: 
                    return mean(array((self.Radius['array'].values())))     
                else:
                    return self.Radius['array'][self.edgeAbscissa]    
    
    def GetLength(self,info):
        '''
        This method returns edge's length
        '''
        if 'value' in self.Length:
            return self.Length['value']
    
    def GetYoungModulus(self,info):
        '''
        This method returns edge's young's modulus
        '''
        if 'value' in self.YoungModulus:
            return self.YoungModulus['value']
            
    
class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass

class XMLValidationError(Error):
    '''
    Exception raised for XML validation failure
    '''
    def __init__(self,xmlschema):
        print "Error, Invalid Network Graph Xml File."
        print xmlschema.error_log
        sys.exit()
        
class WrongXSDPathError(Error):
    '''
    Exception raised if a wrong xsd path is provided.
    '''
    def __init__(self):
        print "Warning, Xml schema file not found.\nNetwork Graph Xml file can not be validated."
        
class NoXSDWarning(Error):
    '''
    Exception raised if no xsd file is provided.
    '''
    def __init__(self):
        print "Warning, XML schema file was not provided.\nNetwork Graph Xml file can not be validated."
        
class LXMLError(Error):
    '''
    Exception raised if lxml package is not installed.
    '''
    def __init__(self):
        print "Warning, Lxml package was not provided.\nNetwork Graph Xml file can not be validated." 