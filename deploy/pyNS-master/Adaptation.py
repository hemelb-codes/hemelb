#!/usr/bin/env python

## Program:   PyNS
## Module:    Adaptation.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224

from numpy.lib.function_base import linspace

class Adaptation(object):
    """
    This class implements an algorithm for vascular remodelling
    due to a wall shear stress trigger (Manini S, CMBBE 2012)
    Each element involved will increase its diameter in response
    to an equation function of wall shear stress peak value.
    This class provides the following methods:
    SetSolutions: a method for setting solution dictionary.
    SetBoundaryConditions: a method for setting boundary conditions input
    SetRefValues: a method for setting wall shear stress peak reference values.
    Adapt: a method for applying adaptation algorithm.
    """

    def __init__(self):
        """Constructor"""
        self.solutions = {} #day:solution
        self.boundaryConditions = None
        self.simulationContext = None
        self.refValues = {} #meshName:refValue
        self.Coeff = 0.05 #numeric coefficient
    
    def SetSolutions(self, day, solution):
        """This method sets solutions dict, {day:solutions}"""
        self.solutions[day] = solution
        
    def SetBoundaryConditions(self, boundaryConditions):
        """This method sets boundaryConditions"""
        self.boundaryConditions = boundaryConditions
        
    def SetSimulationContext(self, simulationContext):
        """This method sets simulationContext"""
        self.simulationContext = simulationContext
    
    def SetRefValues(self, day, networkMesh):
        """
        This method sets wall shear stress peak values 
        from pre-operative simulation or to a user-defined value.
        This values are used as referral wss peak values for adaptation law.
        """
        if day == -1:
            for ent, elList in networkMesh.Entities.iteritems():   
                if ent.Id == "radial" or ent.Id == "ulnar" :
                    for el in elList:                        
                        #taoRef = max(self.solutions[-1].GetWssPeak(el)) #instead of using a value from literature taoRef is equal to pre-operative value.
                        taoRef = 4.
                        self.refValues[el.Name] = taoRef
                if ent.Id == "axillarian" or ent.Id == "brachial":
                    for el in elList:
                        taoRef = 3.
                        self.refValues[el.Name] = taoRef
        if day == 0:
            for ent, elList in networkMesh.Entities.iteritems():
                if ent.Id == "cephalic_vein" or ent.Id == "cubiti_vein" or ent.Id == "basilic_vein" or ent.Id == "subclavian_vein":
                    for el in elList:
                        if ent.Id == "cephalic_vein":
                            taoRef = 2.0
                        if ent.Id == "basilic_vein" or ent.Id == "cubiti_vein":
                            taoRef = 1.0
                        if ent.Id == "subclavian_vein":
                            taoRef = 0.5
                        self.refValues[el.Name] = taoRef    
        
    def Adapt(self, day):
        """This method will apply adaptation law for specified elements."""
        if day>0:
            if day == 1:    
                for el in self.solutions[day-1].NetworkMesh.Elements:
                    if el.Type == "WavePropagation":
                        kd1 = el.Radius[0]
                        kd2 = el.Radius[len(el.Radius)-1]
                        el.dayRadius[day-1]=[kd1,kd2]          
            for el in self.solutions[day-1].NetworkMesh.Elements:
                el.Initialized = False
            preRun = False
            
            for elem in self.solutions[day-1].NetworkMesh.Elements:
                if elem.Type == "Anastomosis":
                    proximalArtery = elem.Proximal
                    distalArtery = elem.Distal
                    proximalVein = elem.Vein
            
            #No adaptation in case of upper arm anastomosis and diabetic patient
            if self.simulationContext.Context["diab"] == 1 and self.simulationContext.Context["ftype"] > 2:
                print "Diabetic patient, no arterial adaptation"
                for ent, elList in self.solutions[day-1].NetworkMesh.Entities.iteritems():  
                    if ent.Id == "axillarian" or ent.Id == "brachial" or ent.Id == "radial" or ent.Id == "ulnar" or ent.Id == "cephalic_vein" or ent.Id == "cubiti_vein" or ent.Id == "basilic_vein" or ent.Id == "subclavian_vein":  
                        for el in elList:
                            kd1_n = el.Radius[0]
                            kd2_n = el.Radius[len(el.Radius)-1]     
                            el.dayRadius[day]=[kd1_n,kd2_n] 
                    if ent.Id == "cephalic_vein" or ent.Id == "cubiti_vein" or ent.Id == "basilic_vein" or ent.Id == "subclavian_vein": 
                        for el in elList:
                            taoRef = self.refValues[el.Name]
                            taoPeaks = self.solutions[day-1].GetWssPeak(el)
                            tao0 = taoPeaks[0]
                            tao1 = taoPeaks[1]
                            tao = linspace(tao0,tao1,len(el.Radius))
                            deltaTao = tao-taoRef
                            k = (1.0+(deltaTao*self.Coeff))
                            
                            if el == proximalVein:
                                #linear adaptation on anastomosis of 20%
                                k2 = 45./44.
                                if day > 10:
                                    k2 = 1.
                                x = linspace(0,len(el.Radius),len(el.Radius))
                                taoProxV = (tao0-tao1)*((x/len(el.Radius))**2-(2.0*(x/len(el.Radius))))+tao0   #y = (k1-k2)(x^2-2x)+k1 (quadratic decreasing wss)
                                tao = taoProxV    
                                deltaTaoProxV = taoProxV-taoRef
                                k = (1.0+(deltaTaoProxV*self.Coeff))
                                kProxV = (k-k2)*((x/len(el.Radius))**2)+k2     #y = (k2-k1)x^2+k1 (quadratic increasing radius)                                                                                
                                k=kProxV
                            if min(tao) > taoRef:
                                el.Radius*=k   
                            kd1_n = el.Radius[0]
                            kd2_n = el.Radius[len(el.Radius)-1]     
                            el.dayRadius[day]=[kd1_n,kd2_n] 
                            
            else:

                for ent, elList in self.solutions[day-1].NetworkMesh.Entities.iteritems():  
                    if ent.Id == "axillarian" or ent.Id == "brachial" or ent.Id == "radial" or ent.Id == "ulnar" or ent.Id == "cephalic_vein" or ent.Id == "cubiti_vein" or ent.Id == "basilic_vein" or ent.Id == "subclavian_vein":  
                        for el in elList:
                            taoRef = self.refValues[el.Name]
                            taoPeaks = self.solutions[day-1].GetWssPeak(el)
                            tao0 = taoPeaks[0]
                            tao1 = taoPeaks[1]
                            tao = linspace(tao0,tao1,len(el.Radius))
                            deltaTao = tao-taoRef
                            k = (1.0+(deltaTao*self.Coeff))
                            
                            if el == proximalArtery: 
                                x = linspace(0,len(el.Radius),len(el.Radius))
                                taoProxA = (tao1-tao0)*((x/len(el.Radius))**2)+tao0     #y = (k2-k1)x^2+k1 (quadratic increasing wss)
                                tao=taoProxA
                                deltaTaoProxA = taoProxA-taoRef
                                k = (1.0+(deltaTaoProxA*self.Coeff))                        
                                kProxA = (k-1.)*((x/len(el.Radius))**2-(2.0*(x/len(el.Radius))))+k   #y = (k1-k2)(x^2-2x)+k1 (quadratic decreasing radius)                                 
                                k=kProxA
                            if el == distalArtery: 
                                x = linspace(0,len(el.Radius),len(el.Radius))
                                taoDistA = (tao0-tao1)*((x/len(el.Radius))**2-(2.0*(x/len(el.Radius))))+tao0   #y = (k1-k2)(x^2-2x)+k1 (quadratic decreasing wss)
                                tao = taoDistA
                                deltaTaoDistA = taoDistA-taoRef
                                k = (1.0+(deltaTaoDistA*self.Coeff))
                                kDistA = (k-1.)*((x/len(el.Radius))**2)+1.     #y = (k2-k1)x^2+k1 (quadratic increasing radius)                                
                                k=kDistA
                            if el == proximalVein:
                                #linear adaptation on anastomosis of 20%
                                k2 = 1.02
                                #k2 = 45./44.
                                if day > 10:
                                    k2 = 1.
                                x = linspace(0,len(el.Radius),len(el.Radius))
                                taoProxV = (tao0-tao1)*((x/len(el.Radius))**2-(2.0*(x/len(el.Radius))))+tao0   #y = (k1-k2)(x^2-2x)+k1 (quadratic decreasing wss)
                                tao = taoProxV    
                                deltaTaoProxV = taoProxV-taoRef
                                k = (1.0+(deltaTaoProxV*self.Coeff))
                                kProxV = (k-k2)*((x/len(el.Radius))**2)+k2     #y = (k2-k1)x^2+k1 (quadratic increasing radius)                                                                                
                                print "PROXV K" , kProxV
                                k=kProxV
                            if min(tao) > taoRef:
                                el.Radius*=k                         
                                
                            kd1_n = el.Radius[0]
                            kd2_n = el.Radius[len(el.Radius)-1]     
                            el.dayRadius[day]=[kd1_n,kd2_n]
                            
                
        if day == 0:
            preRun = True
        if day == -1:
            preRun = False
        return preRun