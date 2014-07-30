#!/usr/bin/env python

## Program:   PyNS
## Module:    SimulationContext.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from xml.etree import ElementTree as etree
import shutil

class SimulationContext(object):
    '''
    This class provides a dictionary for simulation parameters.
    This dictionary is created from XML BoundaryConditions File.
    This class provides the following methods:
    SetEvaluator: a method for setting evaluator class.
    ReadFromXML: a method for reading XML BoundaryConditions File.
    UpdateXML: a method for updating Boundary Conditions XML File Parameters from ModelAdaptor.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.Id = None
        self.Context = {} # dictionary of simulation parameters. {parameter:value}
        self.Defaults = {} # dictionary of defaults simulation parameters. {parameter:value}
        self.Evaluator = None
        
    def SetEvaluator(self, evaluator):
        '''
        This method sets Evaluator input
        '''
        self.Evaluator = evaluator
        
    def ReadFromXML(self, xmlcontextpath, xsdcontextpath = None):
        '''
        This method reads Boundary Conditions XML File.
        '''
        self.xmlcontextpath = xmlcontextpath
        doccontextfile = open(xmlcontextpath)
        contexttree = etree.parse(doccontextfile)
        contextgraph = contexttree.getroot()
        contextgraph_dict = contextgraph.attrib
        self.Id = contextgraph_dict['id']        
        for dataset in contextgraph.findall(".//specificdatasets"):
            for pdata in dataset.findall(".//patient_data"):
                for data in pdata:
                    if data.tag == "dos":
                        self.Context['dos'] = data.text
                    if data.tag == "dob":
                        self.Context['dob'] = data.text
                    if data.tag == "age":
                        for value in data.findall(".//scalar"):
                            self.Context['age'] = float(value.text)
                        for exp in data.findall(".//expression"):
                            self.Context['age'] = exp.text
                    if data.tag == "gender":
                        if data.text == "male":
                            self.Context['gender'] = int(1)
                        if data.text == "female":
                            self.Context['gender'] = int(0)
                        if data.text is None:
                            self.Context['gender'] = None
                    if data.tag == "arm":
                        if data.text == "right":
                            self.Context['arm'] = int(1)
                        if data.text == "left":
                            self.Context['arm'] = int(0)
                        if data.text is None:
                            self.Context['arm'] = None
                    if data.tag == "ftype":
                        if data.text == "Radio-Cephalic EndToEnd":
                            self.Context['ftype'] = int(0)
                        if data.text == "Radio-Cephalic EndToSide":
                            self.Context['ftype'] = int(1)
                        if data.text == "Radio-Cephalic SideToSide":
                            self.Context['ftype'] = int(2)
                        if data.text == "Brachio-Cephalic EndToSide":
                            self.Context['ftype'] = int(3)
                        if data.text == "Brachio-Cephalic SideToSide":
                            self.Context['ftype'] = int(4)
                        if data.text == "Brachio-Basilic EndToSide":
                            self.Context['ftype'] = int(5)
                        if data.text == "Brachio-Basilic SideToSide":
                            self.Context['ftype'] = int(6)
                        if data.text is None:
                            self.Context['ftype'] = None 
                    if data.tag == "height":
                        if data.text:
                            self.Context['height'] = float(data.text)
                        else:
                            self.Context['height'] = None
                    if data.tag == "weight":     
                        if data.text:
                            self.Context['weight'] = float(data.text)
                        else:
                            self.Context['weight'] = None
                    if data.tag == "bsa":
                        for value in data.findall(".//scalar"):
                            self.Context['bsa'] = float(value.text)
                        for exp in data.findall(".//expression"):
                            self.Context['bsa'] = exp.text
                    if data.tag == "sysp":     
                        if data.text:
                            self.Context['sysp'] = float(data.text)
                        else:
                            self.Context['sysp'] = None
                    if data.tag == "diap":     
                        if data.text:
                            self.Context['diap'] = float(data.text)
                        else:
                            self.Context['diap'] = None
                    if data.tag == "mean_pressure":
                        for value in data.findall(".//scalar"):
                            self.Context['mean_pressure'] = float(value.text)
                        for exp in data.findall(".//expression"):
                            self.Context['mean_pressure'] = exp.text
                    
                    if data.tag == "cardiac_output":
                        for value in data.findall(".//scalar"):
                            self.Context['cardiac_output'] = float(value.text)    
                        for exp in data.findall(".//expression"):
                            self.Context['cardiac_output'] = exp.text
                    if data.tag == "period":     
                        if data.text:
                            self.Context['period'] = float(data.text)
                        else:
                            self.Context['period'] = None
                    if data.tag == "brachial_flow":
                        if data.text:   
                            self.Context['brachial_flow'] = float(data.text)
                        else:
                            self.Context['brachial_flow'] = None
                    if data.tag == "radial_flow":     
                        if data.text:
                            self.Context['radial_flow'] = float(data.text)
                        else:
                            self.Context['radial_flow'] = None
                    if data.tag == 'ulnar_flow':     
                        if data.text:
                            self.Context['ulnar_flow'] = float(data.text)
                        else:
                            self.Context['ulnar_flow'] = None
                    if data.tag == 'ht':     
                        if data.text:
                            self.Context['ht'] = float(data.text)
                        else:
                            self.Context['ht'] = None
                    if data.tag == 'cp':     
                        if data.text:
                            self.Context['cp'] = float(data.text)
                        else:
                            self.Context['cp'] = None
                    if data.tag == "dynamic_viscosity":
                        for value in data.findall(".//scalar"):
                            self.Context['dynamic_viscosity'] = float(value.text)
                        for exp in data.findall(".//expression"):
                            self.Context['dynamic_viscosity'] = exp.text
                    if data.tag == 'hyp':     
                        if data.text:
                            if data.text == "no":
                                self.Context['hyp'] = int(0)
                            if data.text == "yes":
                                self.Context['hyp'] = int(1)
                        else:
                            self.Context['hyp'] = None 
                    if data.tag == 'diab':     
                        if data.text:
                            if data.text == "no":
                                self.Context['diab'] = int(0)
                            if data.text == "yes":
                                self.Context['diab'] = int(1)
                        else:
                            self.Context['diab'] = None 
                        
            for parameters in dataset.findall(".//simulation_parameters"):
                for param in parameters:
                    if param.tag == "blood_density":
                        self.Context['blood_density'] = float(param.text)
                    if param.tag == "aorta_ratio":
                        self.Context['aorta_ratio'] = float(param.text)
                    if param.tag == "carotid_ratio":
                        self.Context['carotid_ratio'] = float(param.text)
                    if param.tag == "upper_arm_ratio":
                        self.Context['upper_arm_ratio'] = float(param.text)
                    if param.tag == "vertebral_ratio":
                        self.Context['vertebral_ratio'] = float(param.text)
                    if param.tag == "lower_arm_ratio":
                        self.Context['lower_arm_ratio'] = float(param.text)
                    if param.tag == "upper_arm_2_ratio":
                        self.Context['upper_arm_2_ratio'] = float(param.text)
                    if param.tag == "vein_ratio":
                        self.Context['vein_ratio'] = float(param.text)
                    if param.tag == "poisson_ratio":
                        self.Context['poisson_ratio'] = float(param.text)
                    if param.tag == "timestep":
                        self.Context['timestep'] = float(param.text)
                    if param.tag == "cycles":
                        self.Context['cycles'] = int(param.text)
        
        #Default generic values
        self.Defaults['dos'] = '27/07/2010'
        self.Defaults['dob'] = '27/07/1960'
        self.Defaults['gender'] = int(1)
        self.Defaults['arm'] = int(1)
        self.Defaults['ftype'] = int(1)
        self.Defaults['height'] = float(175)
        self.Defaults['weight'] = float(70)
        self.Defaults['sysp'] = float(110)
        self.Defaults['diap'] = float(70)
        self.Defaults['brachial_flow'] = float(130) 
        self.Defaults['radial_flow'] = float(20)
        self.Defaults['ulnar_flow'] = float(30)
        self.Defaults['period'] = float(1)
        self.Defaults['ht'] = float(45)
        self.Defaults['cp'] = float(7)
        self.Defaults['hyp'] = int(0)
        self.Defaults['diab'] = int(0)
    
    def UpdateXML(self, genericXml, specificXml):
        '''
        This method updates Boundary Conditions XML File Parameters from ModelAdaptor.
        '''
        shutil.copy(genericXml, specificXml)
        self.xmlcontextpath = specificXml
        
        doccontextfile = open(self.xmlcontextpath)
        contexttree = etree.parse(doccontextfile)
        contextgraph = contexttree.getroot()
        
        contextgraph_dict = contextgraph.attrib
        self.Id = contextgraph_dict['id']
  
        for context in contextgraph.findall(".//specificdatasets"):
            for data in context.findall(".//patient_data"):
                for pdata in data:
                    if pdata.tag in self.Context:        
                        if pdata.tag == "age":
                            for exp in pdata.findall(".//expression"):
                                pdata.remove(exp)
                            parameter_value = etree.SubElement(pdata, "scalar")
                            parameter_value.text = str(self.Context[pdata.tag]) 
                        if pdata.tag == "bsa":
                            for exp in pdata.findall(".//expression"):
                                pdata.remove(exp)
                            parameter_value = etree.SubElement(pdata, "scalar")
                            parameter_value.text = str(self.Context[pdata.tag])
                        if pdata.tag == "mean_pressure":
                            for exp in pdata.findall(".//expression"):
                                pdata.remove(exp)
                            parameter_value = etree.SubElement(pdata, "scalar")
                            parameter_value.text = str(self.Context[pdata.tag]) 
                        if pdata.tag == "cardiac_output":
                            for exp in pdata.findall(".//expression"):
                                pdata.remove(exp)
                            parameter_value = etree.SubElement(pdata, "scalar")
                            parameter_value.text = str(self.Context[pdata.tag]) 
                        if pdata.tag == "dynamic_viscosity":
                            for exp in pdata.findall(".//expression"):
                                pdata.remove(exp)
                            parameter_value = etree.SubElement(pdata, "scalar")
                            parameter_value.text = str(self.Context[pdata.tag]) 
                        
                        if pdata.findall(".//scalar"):
                            parameter_value.text = str(self.Context[pdata.tag])
                        else:
                            pdata.text = str(self.Context[pdata.tag])
                            if pdata.tag == "gender":
                                if self.Context[pdata.tag] == 0.0:
                                    pdata.text = 'female'
                                if self.Context[pdata.tag] == 1.0:
                                    pdata.text = 'male'
                            if pdata.tag == "arm":
                                if self.Context[pdata.tag] == 0.0:
                                    pdata.text = 'left'
                                if self.Context[pdata.tag] == 1.0:
                                    pdata.text = 'right'
                            if pdata.tag == "ftype":
                                if self.Context[pdata.tag] == 0.0:
                                    pdata.text = 'Radio-Cephalic EndToEnd'
                                if self.Context[pdata.tag] == 1.0:
                                    pdata.text = 'Radio-Cephalic EndToSide'
                                if self.Context[pdata.tag] == 2.0:
                                    pdata.text = 'Radio-Cephalic SideToSide'
                                if self.Context[pdata.tag] == 3.0:
                                    pdata.text = 'Brachio-Cephalic EndToSide'
                                if self.Context[pdata.tag] == 4.0:
                                    pdata.text = 'Brachio-Cephalic SideToSide'
                                if self.Context[pdata.tag] == 5.0:
                                    pdata.text = 'Brachio-Basilic EndToSide'
                                if self.Context[pdata.tag] == 6.0:
                                    pdata.text = 'Brachio-Basilic SideToSide'
                                if self.Context[pdata.tag] == 7.0:
                                    pdata.text = 'Pre-Surgery'
                            if pdata.tag == "hyp":
                                if self.Context[pdata.tag] == 0.0:
                                    pdata.text = 'no'
                                if self.Context[pdata.tag] == 1.0:
                                    pdata.text = 'yes'
                            if pdata.tag == "diab":
                                if self.Context[pdata.tag] == 0.0:
                                    pdata.text = 'no'
                                if self.Context[pdata.tag] == 1.0:
                                    pdata.text = 'yes'    
            for parameters in context.findall(".//simulation_parameters"):
                for param in parameters:
                    if param.tag in self.Context:
                        param.text = str(self.Context[param.tag])            
     
        xmlcontext = etree.ElementTree(contextgraph)       
        xmlcontext.write(self.xmlcontextpath,encoding='iso-8859-1')
        self.ReadFromXML(self.xmlcontextpath)
        
        #Model Adaptor Parameters
        if self.Context['age'] < 40: 
            self.Context['K_ax'] = 1.34
            self.Context['K_sub'] = 1.675
            self.Context['K_ver'] = 3.3
        if self.Context['age'] >= 40 and self.Context['age'] <= 59: 
            self.Context['K_ax'] = 1.43
            self.Context['K_sub'] = 1.788
            self.Context['K_ver'] = 3.3
        if self.Context['age'] >= 60:
            self.Context['K_ax'] = 1.43   #era 1.64
            self.Context['K_sub'] = 2.05
            self.Context['K_ver'] = 3.6
        if self.Context['age'] <= 70:
            self.Context['K_C'] = 1
        if self.Context['age'] > 70:
            self.Context['K_C'] = 0
        if self.Context['hyp'] == 1 or self.Context['diab'] == 1:
            if self.Context['age'] <= 70:
                self.Context['K_C1'] = 1
            else:
                self.Context['K_C1'] = 0
        if self.Context['hyp'] == 0 and self.Context['diab'] == 0:
            self.Context['K_C1'] = 0	
            
        
class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass