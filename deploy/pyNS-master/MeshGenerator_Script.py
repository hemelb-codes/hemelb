#!/usr/bin/env python

## Program:   PyNS
## Module:    MeshGenerator_Script.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

'''
THIS SCRIPT IS USED BY ARCH NETWORK EDITOR (GUI)
IF YOU ARE USING PYNS WITHOUT THE GUI YOU CAN DELETE THIS FILE
'''

from NetworkGraph import NetworkGraph
from NetworkMesh import NetworkMesh
from MeshGenerator import MeshGenerator
import sys, getopt, os

'''Default Values'''
wdir = 'XML/'   #Working Directory (-w or --wdir)
xdir = 'XML/XSD/' #XSD schema files Working Directory (-x or --xdir)
xsdNet =  'vascular_network_v3.2.xsd' #Vascular Network Graph XSD Schema  (-t or --xsdNet)
xsdMesh = 'vascular_mesh_v2.0.xsd' #Vascular Network Mesh XSD Schema  (-h or --xsdMesh)
ToleranceValue = float(5e-2) #default value for tolerance

try:                                
    opts, args = getopt.getopt(sys.argv[1:], "x:w:i:t:o:h:m:v:", ["xdir=", "wdir=", "xmlNet=", "xsdNet=", "xmlMesh=", "xsdMesh=","method=","tolValue="]) 
except getopt.GetoptError: 
    print "Wrong parameters, please use -shortname parameter or --longname=parameter"                                  
    sys.exit(2)  

for opt, arg in opts:
    if opt in ("-w", "--wdir"):
        wdir = arg
    if opt in ("-x", "--xdir"):
        xdir = arg 
    if opt in ("-i", "--xmlNet"):
        xmlNet = arg
    if opt in ("-t", "--xsdNet"):
        xsdNet = arg
    if opt in ("-o", "--xmlMesh"):
        xmlMesh = arg 
    if opt in ("-h", "--xsdMesh"):
        xsdMesh = arg
    if opt in ("-m", "--method"):
        method = arg
    if opt in ("-v", "--tolValue"):
        ToleranceValue = float(arg)
      
xmlnetpath = os.path.join(wdir, xmlNet)   
xsdnetpath = os.path.join(xdir, xsdNet)
xmlmeshpath = os.path.join(wdir, xmlMesh)
xsdmeshpath = os.path.join(xdir, xsdMesh)

'''Creating NetworkGraph Object From its XML'''
networkGraph = NetworkGraph()
networkGraph = NetworkGraph()
networkGraph.ReadFromXML(xmlnetpath, xsdnetpath)

'''Mesh generation, XML Network Graph is needed for creating XML Network Mesh.'''
meshGenerator = MeshGenerator()
meshGenerator.SetNetworkGraph(networkGraph)
networkMesh = NetworkMesh()
meshGenerator.SetNetworkMesh(networkMesh)
meshGenerator.SetMaxLength(ToleranceValue)
meshGenerator.GenerateMesh()
networkMesh.WriteToXML(xmlmeshpath)