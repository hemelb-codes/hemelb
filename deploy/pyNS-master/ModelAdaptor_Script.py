#!/usr/bin/env python

## Program:   PyNS
## Module:    ModelAdaptor_Script.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

'''
THIS SCRIPT IS USED BY ARCH NETWORK EDITOR (GUI)
IF YOU ARE USING PYNS WITHOUT THE GUI YOU CAN DELETE THIS FILE
'''

from NetworkGraph import NetworkGraph
from SimulationContext import SimulationContext
from ModelAdaptor import ModelAdaptor
from Evaluator import Evaluator
import sys, getopt, os

'''Default Values'''
wdir = 'XML/'   #Working Directory (-w or --wdir)
xdir = 'XML/XSD/' #XSD schema files Working Directory (-x or --xdir)
xsdNet =  'vascular_network_v3.2.xsd' #Vascular Network Graph XSD Schema  (-t or --xsdNet)
xsdBound = 'boundary_conditions_v3.1.xsd' #Boundary Condition XSD Schema  (-h or --xsdBound)

try:                                
    opts, args = getopt.getopt(sys.argv[1:], "x:w:i:t:o:h:", ["xdir=", "wdir=", "xmlNet=", "xsdNet=", "xmlBound=", "xsdBound="]) 
except getopt.GetoptError: 
    print "Wrong parameters, please use -shortname parameter or --longname=parameter"                                  
    sys.exit(2)  

for opt, arg in opts:
    if opt in ("-w", "--wdir"):
        wdir = arg
    if opt in ("-x", "--xdir"):
        xdir = arg 
    if opt in ("-i", "--xmlNet"):
        xmlNetGeneric = arg
    if opt in ("-t", "--xsdNet"):
        xsdNet = arg
    if opt in ("-o", "--xmlBound"):
        xmlBoundGeneric = arg 
    if opt in ("-h", "--xsdBound"):
        xsdBound = arg


'''Setting Simulation Context Parameters for Simulation'''
simulationContext = SimulationContext()
evaluator = Evaluator()
evaluator.SetSimulationContext(simulationContext)
simulationContext.SetEvaluator(evaluator)
    
'''Parameters Model Adaptor'''
modelAdaptor = ModelAdaptor()
modelAdaptor.SetSimulationContext(simulationContext)
modelAdaptor.SetEvaluator(evaluator)

xmlnetpathGeneric = os.path.join(wdir, xmlNetGeneric)
xmlboundpathGeneric = os.path.join(wdir, xmlBoundGeneric)
        
xsdnetpath = os.path.join(xdir, xsdNet)
xsdboundpath = os.path.join(xdir, xsdBound)
print wdir
simulationContext.ReadFromXML(xmlboundpathGeneric, xsdboundpath)

'''Creating NetworkGraph Object From its XML'''
networkGraph = NetworkGraph()
networkGraph.ReadFromXML(xmlnetpathGeneric, xsdnetpath)
modelAdaptor.Idpat=networkGraph.PatientId
netPost = modelAdaptor.Idpat+'_'+xmlNetGeneric
boundPost = modelAdaptor.Idpat+'_BC_'+xmlNetGeneric
xmlnetpath = os.path.join(wdir, netPost)
xmlboundpath = os.path.join(wdir, boundPost)

#starting customization
modelAdaptor.AdaptingParameters(xmlboundpathGeneric,xmlboundpath)

'''NetworkGraph Model Adaptor'''
modelAdaptor.SetNetworkGraph(networkGraph)
evaluator.SetNetworkGraph(networkGraph)
modelAdaptor.AdaptingModel(xmlnetpathGeneric,xmlnetpath)