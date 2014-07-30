#!/usr/bin/env python

## Program:   PyNS
## Module:    Solver_Script.py
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

from Evaluator import Evaluator
from SimulationContext import SimulationContext
from Solver import SolverFirstTrapezoid
from NetworkMesh import NetworkMesh
from NetworkGraph import NetworkGraph
from NetworkSolutions import NetworkSolutions
from BoundaryConditions import BoundaryConditions
from MeshGenerator import MeshGenerator
import sys, getopt, os

'''Default Values'''
wdir = 'XML/'   #Working Directory (-w or --wdir)
xdir = 'XML/XSD/' #XSD schema files Working Directory (-x or --xdir)
xsdNet =  'vascular_network_v3.2.xsd' #Vascular Network Graph XSD Schema  (-t or --xsdNet)
xsdMesh = 'vascular_mesh_v2.0.xsd' #Vascular Network Mesh XSD Schema  (-h or --xsdMesh)
xsdBound = 'boundary_conditions_v3.1.xsd' #Boundary Conditions XSD Schema 
ToleranceValue = float(5e-2) #default value for tolerance
xmlNet = 'v.xml'
xmlBound = 'b.xml'
xmlOut = 'o.xml'
try:                                
    opts, args = getopt.getopt(sys.argv[1:], "x:w:i:t:o:h:m:v:b:d:u:", ["xdir=", "wdir=", "xmlNet=", "xsdNet=", "xmlMesh=", "xsdMesh=","method=","tolValue=","xmlBound=","xsdBound=","xmlOut="]) 
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
    if opt in ("-b", "--xmlBound"):
        xmlBound = arg
    if opt in ("-d", "--xsdBound"):
        xmlBound = arg
    if opt in ("-u", "--xmlOut"):
        xmlOut = arg
  
xmlnetpath = os.path.join(wdir, xmlNet)   
xsdnetpath = os.path.join(xdir, xsdNet)
xsdmeshpath = os.path.join(xdir, xsdMesh)
xmlboundpath = os.path.join(wdir, xmlBound)
xsdboundpath = os.path.join(xdir, xsdBound)
xmloutpath = os.path.join(wdir, xmlOut)
images = wdir+'Images/'
if not os.path.exists (images):
    os.mkdir(images)

'''Creating NetworkGraph Object From its XML'''
networkGraph = NetworkGraph()
networkGraph.ReadFromXML(xmlnetpath, xsdnetpath)

'''Mesh generation, XML Network Graph is needed for creating XML Network Mesh.'''
meshGenerator = MeshGenerator()
meshGenerator.SetNetworkGraph(networkGraph)
networkMesh = NetworkMesh()
meshGenerator.SetNetworkMesh(networkMesh)
meshGenerator.SetMaxLength(ToleranceValue)
meshGenerator.GenerateMesh()

'''Setting Boundary Conditions Mesh input and reading XML Boundary Conditions File'''
simulationContext = SimulationContext()
simulationContext.ReadFromXML(xmlboundpath, xsdboundpath)
evaluator = Evaluator()
evaluator.SetSimulationContext(simulationContext)
simulationContext.SetEvaluator(evaluator)
boundaryConditions = BoundaryConditions()
boundaryConditions.SetSimulationContext(simulationContext)
boundaryConditions.SetNetworkMesh(networkMesh)
boundaryConditions.ReadFromXML(xmlboundpath, xsdboundpath)

'''Setting Evaluator'''
evaluator.SetNetworkGraph(networkGraph)
evaluator.SetNetworkMesh(networkMesh)
preRun = False
for el in networkMesh.Elements:
    if el.Type == 'WavePropagation' and el.nonLinear is True:
        preRun = True
        break

'''Pre-run'''
if preRun is True:
    ''' Setting Solver Class'''
    solver = SolverFirstTrapezoid()  
    solver.SetNetworkMesh(networkMesh)
    solver.SetBoundaryConditions(boundaryConditions)
    solver.SetSimulationContext(simulationContext)
    solver.SetEvaluator(evaluator)
    solver.SetSteadyFlow()
    print "Steady Pre-Run, setting non-linear parameters"
    solver.Solve() 
    parameters = ["Radius","Compliance"]
    for el in networkMesh.Elements:
        el.SetLinearValues(parameters)

'''Run'''
evaluator.ExpressionCache = {}
solver = SolverFirstTrapezoid() 
solver.SetNetworkMesh(networkMesh)
solver.SetBoundaryConditions(boundaryConditions)
solver.SetSimulationContext(simulationContext)
solver.SetEvaluator(evaluator) 
solver.SetPulseFlow()
print "Solving System"
solver.Solve()

'''Post Processing: Setting Solutions input and plotting some information and/or writing solutions to XML Solutions File'''
networkSolutions = NetworkSolutions()
networkSolutions.SetNetworkMesh(networkMesh)
networkSolutions.SetNetworkGraph(networkGraph)
networkSolutions.SetSimulationContext(simulationContext)
networkSolutions.SetSolutions(solver.Solutions)
networkSolutions.SetImagesPath({'im':images})
for element in networkMesh.Elements:
    if element.Type == 'WavePropagation':
        networkSolutions.PlotFlow(element.Id)
        networkSolutions.PlotPressure(element.Id)
networkSolutions.WriteToXML(xmloutpath)