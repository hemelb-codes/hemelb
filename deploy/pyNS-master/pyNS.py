#!/usr/bin/env python

## Program:   PyNS
## Module:    pyNS.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from ModelAdaptor import ModelAdaptor
from NetworkGraph import NetworkGraph
from NetworkMesh import NetworkMesh
from MeshGenerator import MeshGenerator
from BoundaryConditions import BoundaryConditions
from Solver import SolverFirstTrapezoid
from NetworkSolutions import NetworkSolutions
from SimulationContext import SimulationContext
from Evaluator import Evaluator
from Adaptation import Adaptation, linspace
from Export import export as exporting
import os, sys, shutil, SimpleHTTPServer, SocketServer, webbrowser, time, argparse

def mylistdir(directory):
    '''
    A specialized version of os.listdir() that ignores files that
    start with a leading period.
    '''
    filelist = os.listdir(directory)
    return [x for x in filelist
            if not (x.startswith('.'))]

def runSimulation(simType, defaultNet, wdir, odir, images, xsd, net, mesh, xmlout, bound, netSchema, boundSchema, template, parameters, diameters, days, xmlSol, xmlMesh, writeCsv, plotImages, plotPressure, plotFlow, plotWss, plotReynolds, writePressure, writeFlow, writeWss, writeReynolds, velocityProfile, results, excludeWss, export, automaticResults, inputGnuid):
    
    '''Welcome and instructions messages.'''
    
    print "##########################################"
    print "############ Welcome to pyNS #############"
    print "## ./pyNS -h or --help for instructions ##"
    print "##########################################\n"
            
    '''Exporting results into txt files'''   
    if export is not False:
        if not os.path.exists ('Results/%s/exportedSolutions' % export):
            os.mkdir('Results/%s/exportedSolutions' % export)
        for f in mylistdir('Results/%s/json' % export):
            if f == 'info.json':
                pass
            else:
                print "exporting Results/%s/json/" % export + f
                exporting('Results/%s/json/' % export + f)
                new_file = f.split('.')[0]+'.txt'
                shutil.move('Results/%s/json/' % export + new_file, 'Results/%s/exportedSolutions/' % export + new_file)
        sys.exit('All %s solutions exported successfully in Results/%s/exportedSolutions/ folder' % (export,export))
    
    if not results:
        if defaultNet is True:
            simType = 'specific'
            net = 'vascular_network_arterial_right_arm.xml'
            bound = 'boundary_conditions_arterial_right_arm.xml'
        elif template == 'willis':
            simType = 'specific'
            wdir = 'XML/Models/WillisCircle'
            net = 'vascular_network_willis.xml'
            bound = 'boundary_conditions_willis.xml'
        elif simType == 'specific':
            if net is None and bound is not None:
                sys.exit("Please provide a network graph XML input file or choose a generic simulation type.")
            elif net is not None and bound is None:
                sys.exit("Please provide a boundary conditions XML input file or choose a generic simulation type.")
            elif net is None and bound is None:
                sys.exit("Please provide either a network graph XML input file and a boundary conditions XML input file or choose a generic simulation type.")
    
    '''Checking matplotlib module for optional plotting methods.'''
    if plotImages or plotFlow or plotPressure or plotWss or plotReynolds or velocityProfile is True:
        try:
            import matplotlib
        except ImportError:
            sys.exit('Matplotlib package is required for plotting solutions in .png files or computing velocityProfile videos.\nPlease download matplotlib from matplotlib.sourceforge.net.')
            
    '''Loading previous specific results.'''
    if results is not False:
        while True:
            print "Starting webServer for post-processing results. Close it with CTRL-C."
            Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
            try:
                port = 8000
                httpd = SocketServer.TCPServer(("localhost", port), Handler)
            except:
                try:
                    pid = None
                    for line in os.popen("lsof -i:8000"):
                        fields = line.split()
                        pid = fields[1]
                    if pid:
                        os.system("kill %s" %pid)
                        time.sleep(5)
                    httpd = SocketServer.TCPServer(("localhost", port), Handler)
                except:
                    connected = False
                    startPort = 8000
                    while not connected:
                        try:
                            httpd = SocketServer.TCPServer(("localhost", startPort), Handler)
                            connected = True
                            port = startPort
                        except:
                            startPort+=1
                    
            if results == 'last':
                ip = "http://localhost:%s" %port
                webbrowser.open_new_tab(ip+'/Results/results.html')
            else:
                if os.path.exists('Results/'+results):
                    ip = "http://localhost:%s" %port
                    webbrowser.open_new_tab(ip+"/Results/"+results+"/results.html")
                else:
                    sys.exit('Error: '+results+' directory does not exist.')
            httpd.serve_forever()
        
    '''Checking for webserver instance'''
    if automaticResults:
        try:
            ip = "http://localhost:8000"
            pid = None
            for line in os.popen("lsof -i:8000"):
                fields = line.split()
                pid = fields[1]
            if pid:
                os.system("kill %s" %pid)
            Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
            httpd = SocketServer.TCPServer(("localhost", 8000), Handler)
        except:
            connected = False
            startPort = 8000
            while not connected:
                try:
                    Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
                    httpd = SocketServer.TCPServer(("localhost", startPort), Handler)
                    connected = True
                    port = startPort
                    ip = "http://localhost:%s" %port
                except:
                    startPort+=1
    
    '''SIMULATION'''
    
    '''Create XML and image directories'''
    if not os.path.exists (wdir):
        os.mkdir(wdir)
    if not os.path.exists (xsd):
        os.mkdir(xsd)

    '''If needed, creating output directory(s).'''
    if xmlSol is True or xmlMesh is True or writeFlow is True or writePressure is True or  writeWss is True or writeReynolds is True:
        if not os.path.exists (odir):
            os.mkdir(odir)
    if writeFlow is True:
        ofdir = os.path.join(odir, 'Flow/')
        if not os.path.exists (ofdir):
            os.mkdir(ofdir)
    if writePressure is True:
        opdir = os.path.join(odir, 'Pressure/')
        if not os.path.exists (opdir):
            os.mkdir(opdir)
    if writeWss is True:
        owdir = os.path.join(odir, 'Wss/')
        if not os.path.exists (owdir):
            os.mkdir(owdir)
    if writeReynolds is True:
        oodir = os.path.join(odir, 'Other/')
        if not os.path.exists (oodir):
            os.mkdir(oodir)

    '''If needed, creating images directory.'''
    if plotImages is True:
        f_images = os.path.join(images, 'Flow/')
        p_images = os.path.join(images, 'Pressure/')
        w_images = os.path.join(images, 'Wss/')
        o_images = os.path.join(images, 'Other/')
        if not os.path.exists (images):
            os.mkdir(images)
            os.mkdir(f_images)
            os.mkdir(p_images)
            os.mkdir(w_images)
            os.mkdir(o_images)

    '''Setting variables.'''
    testTube = 'XML/TEST/CircularStraightTube/'
    netTube = 'vascular_network_v3.0_TUBE.xml'
    boundTube = 'boundary_conditions_v2.0_TUBE.xml'
    testTape = 'XML/TEST/CircularTaperedTube/'
    netTape = 'vascular_network_v3.0_TAPE.xml'
    boundTape = 'boundary_conditions_v2.0_TAPE.xml'
    testSimple = 'XML/TEST/SimpleNetwork/'
    netSimple = 'vascular_network_simple.xml'
    boundSimple = 'boundary_conditions_simple.xml'
    testing = 'XML/TEST/Testing/'
    testingNetwork = 'vascular_network_test.xml'
    testingBoundary = 'boundary_conditions_test.xml'

    if simType == 'specific':
        xmlnetpath = os.path.join(wdir, net)
        xmlboundpath = os.path.join(wdir, bound)
        preRun = True
    if simType == 'tube':
        xmlnetpath = os.path.join(testTube,netTube)
        xmlboundpath = os.path.join(testTube, boundTube)
        preRun = False
    if simType == 'tape':
        xmlnetpath = os.path.join(testTape,netTape)
        xmlboundpath = os.path.join(testTape, boundTape)
        preRun = False
    if simType == 'simple':
        xmlnetpath = os.path.join(testSimple,netSimple)
        xmlboundpath = os.path.join(testSimple, boundSimple)
        preRun = False
    if simType == 'testing':
        xmlnetpath = os.path.join(testing,testingNetwork)
        xmlboundpath = os.path.join(testing, testingBoundary)
        preRun = False
  
    xmlmeshpath = os.path.join(wdir, mesh)
    xmloutpath = os.path.join(odir, xmlout)
    xsdnetpath = os.path.join(xsd, netSchema)
    xsdboundpath = os.path.join(xsd, boundSchema)

    '''Setting adaptation and simulation days'''
    adaptation = Adaptation()
    daysList = map(int,list(linspace(-1,days,days+2)))
    if excludeWss is True and days > 0:
        sys.exit("Error: You can't exclude Wss computing for adaptation algorithm")
 
    '''Setting Simulation Context Parameters for Simulation'''
    simulationContext = SimulationContext()
    evaluator = Evaluator()
    evaluator.SetSimulationContext(simulationContext)
    simulationContext.SetEvaluator(evaluator)

    for day in daysList:
        if day <= 0:
            '''Parameters Model Adaptor'''
            if simType == 'generic':
                modelAdaptor = ModelAdaptor()
                modelAdaptor.SetSimulationContext(simulationContext)
                modelAdaptor.SetEvaluator(evaluator)
                modelAdaptor.ChoosingTemplate(parameters)
                if template == 'arm':
                    if day == -1:
                        modelAdaptor.ftype = 7
                    if modelAdaptor.arm == 0:
                        if modelAdaptor.ftype == 0:
                            wdir = 'XML/Models/Left_Arm/#0.Lower_RC_EE'
                            preRun = True
                        if modelAdaptor.ftype == 1:
                            wdir = 'XML/Models/Left_Arm/#1.Lower_RC_ES'
                            preRun = True
                        if modelAdaptor.ftype == 2:
                            pass
                        if modelAdaptor.ftype == 3:
                            wdir = 'XML/Models/Left_Arm/#3.Upper_BC_ES'
                            preRun = True
                        if modelAdaptor.ftype == 4:
                            pass
                        if modelAdaptor.ftype == 5:
                            wdir = 'XML/Models/Left_Arm/#5.Upper_BB_ES'
                            preRun = True
                        if modelAdaptor.ftype == 6:
                            pass
                        if modelAdaptor.ftype == 7:
                            wdir = 'XML/Models/Left_Arm/PRE'
                            preRun = False
                    if modelAdaptor.arm == 1:
                        if modelAdaptor.ftype == 0:
                            wdir = 'XML/Models/Right_Arm/#0.Lower_RC_EE'
                            preRun = True
                        if modelAdaptor.ftype == 1:
                            wdir = 'XML/Models/Right_Arm/#1.Lower_RC_ES'
                            preRun = True
                        if modelAdaptor.ftype == 2:
                            pass
                        if modelAdaptor.ftype == 3:
                            wdir = 'XML/Models/Right_Arm/#3.Upper_BC_ES'
                            preRun = True
                        if modelAdaptor.ftype == 4:
                            pass
                        if modelAdaptor.ftype == 5:
                            wdir = 'XML/Models/Right_Arm/#5.Upper_BB_ES'
                            preRun = True
                        if modelAdaptor.ftype == 6:
                            pass
                        if modelAdaptor.ftype == 7:
                            wdir = 'XML/Models/Right_Arm/PRE'
                            preRun = False
                    
                netPostGeneric = 'vascular_network.xml'
                boundPostGeneric = 'boundary_conditions.xml'
                netPost = modelAdaptor.Idpat+'_vascular_network.xml'
                boundPost = modelAdaptor.Idpat+'_boundary_conditions.xml'
                xmlnetpathGeneric = os.path.join(wdir, netPostGeneric)
                xmlboundpathGeneric = os.path.join(wdir, boundPostGeneric)
                xmlnetpath = os.path.join(wdir, netPost)
                xmlboundpath = os.path.join(wdir, boundPost)
                simulationContext.ReadFromXML(xmlboundpathGeneric, xsdboundpath)
            else:  
                simulationContext.ReadFromXML(xmlboundpath, xsdboundpath)
            
            if simType == 'generic':  
                modelAdaptor.SettingParameters(parameters)
                modelAdaptor.AdaptingParameters(xmlboundpathGeneric,xmlboundpath)
            
            '''Creating NetworkGraph Object From its XML'''
            networkGraph = NetworkGraph()
            if simType == 'generic':
                networkGraph.ReadFromXML(xmlnetpathGeneric, xsdnetpath)
            else:
                networkGraph.ReadFromXML(xmlnetpath, xsdnetpath)
            
            '''NetworkGraph Model Adaptor'''
            if simType == 'generic':
                modelAdaptor.SetNetworkGraph(networkGraph)
                evaluator.SetNetworkGraph(networkGraph)
                if diameters is False:
                    csvfilepath = modelAdaptor.AdaptingModel(xmlnetpathGeneric,xmlnetpath)
                else:
                    csvfilepath = modelAdaptor.AdaptingModel(xmlnetpathGeneric,xmlnetpath,diameters)   
                       
            '''Setting results directory based on PatientID in networkGraph XML file'''
            
            if plotImages is False:
                try:
                    shutil.rmtree('Results/json')
                except:
                    pass
                try:
                    os.mkdir('Results/json')
                except:
                    pass
                if simType == 'generic':
                    idPat = modelAdaptor.Idpat
                elif template == 'willis':
                    idPat = template
                else:
                    idPat = simType
                if os.path.exists('Results/%s' % idPat):
                    pass
                else:
                    os.mkdir('Results/%s' % idPat)
                    os.mkdir('Results/%s/json' % idPat)
                    shutil.copytree('Results/css','Results/%s/css'  % idPat)
                    shutil.copytree('Results/js','Results/%s/js'  % idPat)
                    shutil.copy('Results/results.html','Results/%s/results.html'  % idPat)

            '''Mesh generation, XML Network Graph is needed for creating XML Network Mesh.'''
            meshGenerator = MeshGenerator()
            meshGenerator.SetNetworkGraph(networkGraph)
            networkMesh = NetworkMesh()
            meshGenerator.SetNetworkMesh(networkMesh)
            meshGenerator.SetMaxLength(5.0e-2)
            meshGenerator.GenerateMesh()
            
        '''Setting Boundary Conditions Mesh input and reading XML Boundary Conditions File'''
        boundaryConditions = BoundaryConditions()
        boundaryConditions.SetSimulationContext(simulationContext)
        boundaryConditions.SetNetworkMesh(networkMesh)
        boundaryConditions.ReadFromXML(xmlboundpath, xsdboundpath)
        boundaryConditions.SetSpecificCardiacOutput()
        
        
        '''In case of a generic simulation, patient-specific generated files will be moved to Results folder.'''
        if simType == 'generic' and day < 0:
            shutil.move(os.path.abspath(xmlnetpath),('Results/%s/%s_pre_vascular_network.xml' % (idPat,idPat)))
            shutil.move(os.path.abspath(xmlboundpath),('Results/%s/%s_pre_boundary_conditions.xml' % (idPat,idPat)))
            shutil.move(os.path.abspath(csvfilepath),('Results/%s/%s_pre_patient_specific.csv' % (idPat,idPat)))
        if simType == 'generic' and day == 0:
            shutil.copy(os.path.abspath(xmlnetpath),('Results/%s/%s_post_vascular_network.xml' % (idPat,idPat)))
            shutil.copy(os.path.abspath(xmlboundpath),('Results/%s/%s_post_boundary_conditions.xml' % (idPat,idPat)))
            shutil.copy(os.path.abspath(csvfilepath),('Results/%s/%s_post_patient_specific.csv' % (idPat,idPat)))
        if simType == 'generic' and day > 0 and day == days:
            shutil.move(os.path.abspath(xmlnetpath),('Results/%s/%s_adapted_vascular_network.xml' % (idPat,idPat)))
            shutil.move(os.path.abspath(xmlboundpath),('Results/%s/%s_adapted_boundary_conditions.xml' % (idPat,idPat)))
            shutil.move(os.path.abspath(csvfilepath),('Results/%s/%s_adapted_patient_specific.csv' % (idPat,idPat)))
        
        '''Setting Evaluator'''
        evaluator.SetNetworkGraph(networkGraph)
        evaluator.SetNetworkMesh(networkMesh)

        '''Adaptation Model'''
        adaptation.SetBoundaryConditions(boundaryConditions)
        adaptation.SetSimulationContext(simulationContext)
        preRun = adaptation.Adapt(day)
        if len(daysList)==1:
            pass
        else:
            print "Day %d " %(day*10)  	#1 step represent 10 days

        ''' Setting Solver Class'''
        solver = SolverFirstTrapezoid()  
        solver.SetNetworkMesh(networkMesh)
        solver.SetBoundaryConditions(boundaryConditions)
        solver.SetSimulationContext(simulationContext)
        solver.SetEvaluator(evaluator)
    
        '''Pre-run'''
        if preRun is True:
            solver.SetSteadyFlow()
            print "Steady Pre-Run, setting non-linear parameters"
            solver.Solve()
            parametersToLinear = ["Radius","Compliance"]
            for el in networkMesh.Elements:
                el.SetLinearValues(parametersToLinear)
            networkMesh.checkLinearConsistence()
    
        '''Run'''
        evaluator.ExpressionCache = {}
        solver = SolverFirstTrapezoid()
        solver.SetNetworkMesh(networkMesh)
        solver.SetBoundaryConditions(boundaryConditions)
        solver.SetSimulationContext(simulationContext)
        solver.SetEvaluator(evaluator) 
        solver.SetPulseFlow()
        print "Solving system"
        solver.Solve()

        '''Post Processing: Setting Solutions input and plotting some information and/or writing solutions to XML Solutions File'''
        '''User can choose two different post processing strategies. Saving images using matplotlib or visualize results in its browser'''

        '''If needed, pyNS writes xml mesh file'''
        if xmlMesh is True:
            meshdirpath = os.path.join(odir,str(day))
            if not os.path.exists(meshdirpath):
                os.mkdir(meshdirpath)
            xmlmeshpath = os.path.join(meshdirpath,mesh)
            outdirpath = os.path.join(odir,str(day))
            if not os.path.exists(outdirpath):
                os.mkdir(outdirpath)
            xmloutpath = os.path.join(outdirpath,xmlout)
            networkMesh.WriteToXML(xmlmeshpath)
    
        '''Setting NetworkSolutions'''
        print "->100%, Running post-processing"
        networkSolutions = NetworkSolutions()
        networkSolutions.SetNetworkMesh(networkMesh)
        networkSolutions.SetNetworkGraph(networkGraph)
        networkSolutions.SetSimulationContext(simulationContext)
        networkSolutions.SetSolutions(solver.Solutions) 
        networkSolutions.WriteJsonInfo(days,networkMesh.Elements,idPat)
        adaptation.SetSolutions(day, networkSolutions)
        adaptation.SetRefValues(day, networkMesh)
    
        '''If needed, pyNS creates images subdirectory(s) for each adaptation step.'''
        if plotImages is True:
            daystr = str(day)+'/'
            f_dayImages = os.path.join(f_images,daystr)   
            p_dayImages = os.path.join(p_images,daystr)
            w_dayImages = os.path.join(w_images,daystr)
            o_dayImages = os.path.join(o_images,daystr)
            if not os.path.exists(images):
                os.mkdir(images)
            if not os.path.exists(f_dayImages):
                os.mkdir(f_dayImages)
            if not os.path.exists(p_dayImages):
                os.mkdir(p_dayImages)
            if not os.path.exists(w_dayImages):
                os.mkdir(w_dayImages)
            if not os.path.exists(o_dayImages):
                os.mkdir(o_dayImages)
            networkSolutions.SetImagesPath({'im':images,'f':f_dayImages,'p':p_dayImages,'w':w_dayImages,'o':o_dayImages})    
        
        '''If needed, pyNS creates output subdirectory(s) for each adaptation step.'''       
        if writeFlow is True:
            if day == -1:
                daystr = 'pre/'
            else:
                daystr = str(day)+'/'
            f_dayOutput = os.path.join(ofdir,daystr) 
            if not os.path.exists(f_dayOutput):
                os.mkdir(f_dayOutput)
        if writePressure is True:
            if day == -1:
                daystr = 'pre/'
            else:
                daystr = str(day)+'/'
            p_dayOutput = os.path.join(opdir,daystr) 
            if not os.path.exists(p_dayOutput):
                os.mkdir(p_dayOutput)
        if writeWss is True:
            if day == -1:
                daystr = 'pre/'
            else:
                daystr = str(day)+'/'
            w_dayOutput = os.path.join(owdir,daystr)
            if not os.path.exists(w_dayOutput):
                os.mkdir(w_dayOutput)
        if writeReynolds is True:
            if day == -1:
                daystr = 'pre/'
            else:
                daystr = str(day)+'/'
            o_dayOutput = os.path.join(oodir,daystr) 
            if not os.path.exists(o_dayOutput):
                os.mkdir(o_dayOutput)
        
        '''If needed, pyNS writes xml Solution file.'''
        if xmlSol is True:
            networkSolutions.WriteToXML(xmloutpath)
    
        '''Post process solution for each element of the network'''  
        for element in networkMesh.Elements:  
            if element.Type == 'WavePropagation' or element.Type == 'Resistance':
                networkSolutions.WriteJson(element.Id, day, excludeWss, idPat)
                if velocityProfile is True:
                    networkSolutions.SaveVelocityProfile(element,str(day))
                if plotFlow is True:
                    networkSolutions.PlotFlow(element.Id)
                if plotPressure is True:
                    networkSolutions.PlotPressure(element.Id)
                if plotWss is True:
                    networkSolutions.PlotWSS(element)
                if plotReynolds is True:
                    networkSolutions.PlotReynolds(element.Id)
                if writeFlow is True:
                    networkSolutions.WriteFlowOutput(element.Id,f_dayOutput+'Flow_'+element.Name+'.txt')
                if writePressure is True:
                    networkSolutions.WritePressureInput(element.Id,p_dayOutput+'/p_in_'+element.Name+'.txt')
                    networkSolutions.WritePressureOutput(element.Id,p_dayOutput+'/p_out_'+element.Name+'.txt')
                    networkSolutions.WritePressureDrop(element.Id,p_dayOutput+'/p_drop_'+element.Name+'.txt')
                if writeWss is True:
                    networkSolutions.WriteWSSOutput(element.Id,w_dayOutput+'WSS_'+element.Name+'.txt')
                if writeReynolds is True:
                    networkSolutions.WriteReynolds(element.Id,o_dayOutput+'Reynolds'+element.Name+'.txt')
                
    '''Adaptation data'''
    if days > 0:
        networkSolutions.WriteJsonAdapt(adaptation, idPat)
        if writeCsv is True:
            networkSolutions.WriteToCsv(adaptation, 'Diameter')
            networkSolutions.WriteToCsv(adaptation, 'Pressure')
            networkSolutions.WriteToCsv(adaptation, 'Flow')
            networkSolutions.WriteToCsv(adaptation, 'Wss')
    
    '''Export GNUID'''
    if inputGnuid:
        networkSolutions.GetGnuidInformation(idPat, inputGnuid)
     
    print "\nJOB FINISHED"
    if automaticResults:
        try:
            shutil.copytree('Results/%s/json' % idPat,'Results/json',symlinks=True)
        except OSError:
            shutil.rmtree('Results/json')
            shutil.copytree('Results/%s/json' % idPat,'Results/json',symlinks=True)
        print "Starting webServer for post-processing results. Close it with CTRL-C."
        webbrowser.open_new_tab(ip+'/Results/results.html')
        httpd.serve_forever()
        

if __name__ == "__main__":
        
    '''Command-line arguments.'''
   
    parser = argparse.ArgumentParser(description='pyNS, The python Network Solver')
    
    parser.add_argument("-s", "--simType", action="store",dest='simType', default="specific",
					          help="Simulation type, 'generic': fromGenericTemplate. 'specific':from specific xml file. 'tube':circular straight tube simulation. 'tape':circular tapered tube simulation. 'simple': simple network simulation. 'testing' : simple network of arteries vein and resistances.")
    parser.add_argument("--default", action="store_true", dest='defaultNet', default=False,
                    help = "Turn on this flag for simulating an example network representing an arterial network of a right arm.")
    parser.add_argument("-w", "--workingDir", action="store", dest='wdir', default='XML/',
	                  help = "Working directory path for xml input files. By default is located in 'XML/' pyNS subfolder.")
    parser.add_argument("-o", "--outputDir", action="store", dest='odir', default='Output/',
	                  help = "Output directory for subfolders and output files. By default is located in 'Output/' pyNS subfolder.")
    parser.add_argument("-i", "--imagesDir", action="store", dest='images', default='Images/',
					          help = "Images directory for subfolders and output images. By default is located in 'Images/' pyNS subfolder.")
    parser.add_argument("-x", "--xsdDir", action="store", dest='xsd', default = 'XML/XSD/',
                    help="XML schema files directory. By default is located in XML/XSD/ pyNS subfolder.")
    parser.add_argument("-n", "--net", action="store", dest='net', default = None,
	                  help="PreOperative vascular network xml file. By default a right arm case arterial network is loaded.")
    parser.add_argument("-m", "--mesh", action="store", dest='mesh', default = 'vascular_mesh_v1.1.xml',
                    help="Vascular network xml mesh file name. By default is specified as 'vascular_mesh_v1.1.xml'.")
    parser.add_argument("-l", "--xmlOut", action="store", dest="xmlout", default = 'vascular_output.xml',
			              help="Vascular network xml output solutions file name. By default is specified as 'vascular_output.xml'.")
    parser.add_argument("-b", "--bound", action="store", dest='bound', default = None,
			              help="Boundary conditions xml file for a preOperative simulation. By default a standard preOperative boundary condition file associated to default right arm case arterial network is loaded.")
    parser.add_argument("-c", "--netSchema", action="store", dest='netSchema', default = 'vascular_network_v3.2.xsd',
	                  help="Vascular network xml schema xsd file. By default is defined as 'vascular_network_v3.2.xsd' and located in the XML schema files directory.")
    parser.add_argument("-f", "--boundSchema", action="store", dest='boundSchema', default = 'boundary_conditions_v3.1.xsd',
	                  help="Boundary conditions xml schema xsd file. By default is defined as 'boundary_conditions_v3.1.xsd' and located in the XML schema files directory.")
    parser.add_argument("-g", "--template", action="store", dest='template', default = 'arm',
                    help="Specify a template network by choosing between currently implemented models: 'arm', 'willis'")
    parser.add_argument("-k", "--parameters", action="store", dest='parameters', default = 'XML/parameters.csv',
	                  help="Additional .csv file for patient-specific parameters. This allows the generation of a patient-specific network from a generic template. By default is located in 'XML/' pyNS subfolder.")
    parser.add_argument("-d", "--diameters", action="store", dest='diameters', default = False,
	                  help="Additional .csv file for patient-specific measured diameters. This enhance the patient-specific network generated from a generic template. By default does not exist.")
    parser.add_argument("-a", "--adaptation", action="store", dest='adaptation', type = int, default = -1,
	                  help="Turn on adaptation algorithm by setting the number of simulated steps. 1 step represents 10 days. By default simulation is performed for preoperative(-1day)")
    parser.add_argument("--xmlSolution", action="store_true", dest='xmlSol', default = False,
	                  help="Network Graph solution XML file will be saved in the Output directory if this feature is active. By default this feature is inactive.")
    parser.add_argument("--xmlMesh", action="store_true", dest='xmlMesh', default = False,
	                  help="Network Mesh XML file will be saved in the Output directory if this feature is active. By default this feature is inactive.")
    parser.add_argument("--writeCsv", action="store_true", dest='writeCsv', default = False,
	                  help="Adaptation results (flow rate, diameter, pressure and wss) will be saved in separates csv files if this feature is active. By default this feature is inactive.")
    parser.add_argument("-p", "--plotImages", action="store_true", dest='plotImages', default = False,
	                  help="Plot images using matplotlib library instead of using results.html. By default this feature is inactive.")
    parser.add_argument("--plotPressure", action="store_true", dest='plotPressure', default = False,
	                  help="Plot pressure solution for each element of the network. By default this feature is inactive.")
    parser.add_argument("--plotFlow", action="store_true", dest='plotFlow', default = False,
	                  help="Plot flow volume solution for each element of the network. By default this feature is inactive.")
    parser.add_argument("--plotReynolds", action="store_true", dest='plotReynolds', default = False,
	                  help="Plot Reynolds number solution for each element of the network. By default this feature is inactive.")
    parser.add_argument("--plotWss", action="store_true", dest='plotWss', default = False,
	                  help="Plot wall shear stress solution for each element of the network. By default this feature is inactive.")
    parser.add_argument("--writePressure", action="store_true", dest='writePressure', default = False,
	                  help="Write pressure solution for each element of the network in a .txt file. By default this feature is inactive.")
    parser.add_argument("--writeFlow", action="store_true", dest='writeFlow', default = False,
	                  help="Write flow volume solution for each element of the network in a .txt file. By default this feature is inactive.")
    parser.add_argument("--writeReynolds", action="store_true", dest='writeReynolds', default = False,
	                  help="Write Reynolds number solution for each element of the network in a .txt file. By default this feature is inactive.")
    parser.add_argument("--writeWss", action="store_true", dest='writeWss', default = False,
	                  help="Write wall shear stress solution for each element of the network in a .txt file. By default this feature is inactive.")
    parser.add_argument("--velocityProfile", action="store_true", dest='velocityProfile', default = False,
	                  help="Save velocity profile in a .avi file. By default this feature is inactive.")
    parser.add_argument("--results", action="store", dest='results', default = False,
                    help="If active pyNS will be launched in post-processing mode for inspecting existing results. If you want to load a specific result, please specify the name ")
    parser.add_argument("--excludeWss", action="store_true", dest='excludeWss', default = False,
                    help="If active pyNS will not compute wall shear stress improving computational time. For vascular adaptation algorithm excluding wss calculation is not admitted.")
    parser.add_argument("--export", action="store", dest='export', default = False,
                    help="If active pyNS will export to a .txt file the solution relative to the choosen mesh. Please specify a mesh name.")
    parser.add_argument("--noAutomaticResults", action="store_false", dest='automaticResults', default = True,
                    help="In case of multiple concurrent simulations it is reccomended to use this flag for avoiding the automatic load of pyNS local web server for postProcessing. User can manually inspect results by open results.html file located in the corresponding folder related to the simulation's patient_id.")
    parser.add_argument("--inputGnuid", action="store", dest='inputGnuid', default = None,
                    help="If active pyNS will compute information for Gnuid solver in a txt file for properly setting boundaryConditions for 3D CFD simulations using gnuid. Please specify the radius of the 3D mesh [cm]. By default this feature is Inactive.")
                      
    args = parser.parse_args()
        
    try:
        runSimulation(args.simType, args.defaultNet, args.wdir, args.odir, args.images, args.xsd, args.net, args.mesh, args.xmlout, args.bound, args.netSchema, args.boundSchema, args.template, args.parameters, args.diameters, args.adaptation, args.xmlSol, args.xmlMesh, args.writeCsv, args.plotImages, args.plotPressure, args.plotFlow, args.plotWss, args.plotReynolds, args.writePressure, args.writeFlow, args.writeWss, args.writeReynolds, args.velocityProfile, args.results, args.excludeWss, args.export, args.automaticResults, args.inputGnuid)
    except KeyboardInterrupt:
        print "\nLocal web server for post processing was shutdown successfully. pyNS is ready for next simulation."
        print "If you want to inspect last simulation results, type ./pyNS.py --results last"
        print "If you want to inspect a specific simulation results, type ./pyNS.py --results folderName"
