#!/usr/bin/env python

## Program:   PyNS
## Module:    Evaluator.py
## Language:  Python
## Date:      $Date: 2012/09/04 10:21:12 $
## Version:   $Revision: 0.4.2 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

import re
from math import *
from datetime import *

class Evaluator(object):
    '''
    This class evaluates xml expression using regular expressions.
    This class provides the following methods:
    SetSimulationContext : a method for setting simulation context.
    SetNetworkGraph: a method for setting NetworkGraph input.
    SetNetworkMesh : a method for setting NetworkMesh.
    SetInfo: a method for setting info dictionary ({'DofMap':self.DofMap}...)
    SetAbscissa: a method for setting abscissa value.
    GetElement: a method for finding element from its edge and specified position (abscissa).
    GetEdge: a method for returning corresponding edge or superedge.
    GetVariableComponents: a method for splitting expression into variables.
    Evaluate: the main method of the class, it evaluates expression and returns result.
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.variableRe = re.compile('\$.*?[\}\]]')
        self.parameterRe = re.compile('\$.*?[\{\[]')
        self.elementRe = re.compile('\[.*?\]')
        self.edgeRe = re.compile('\{.*?\}')
        self.abscissa = None
        self.NetworkMesh = None
        self.NetworkGraph = None
        self.SimulationContext = None
        self.Info = None
        self.ExpressionCache = {}
        self.rhsCache = {}
        
    def SetSimulationContext(self, context):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = context
        
    def SetNetworkMesh(self, networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetNetworkGraph(self, networkGraph):
        '''
        Setting NetworkMesh
        '''
        self.NetworkGraph = networkGraph
        
    def SetInfo(self,info):
        '''
        Setting info dictionary
        '''
        self.Info = info
         
    def SetAbscissa(self,abscissa):
        '''
        Setting Abscissa Value
        '''
        self.abscissa = abscissa
        
    def GetElement(self, edge, abscissa = None):
        '''
        This method returns the element from its edge and
        specified position (abscissa)
        '''
        if abscissa is None:
            abscissa = 0.0
        if self.abscissa is None:
            pass
        else:
            abscissa = self.abscissa
        try:             
            ed = self.NetworkGraph.Edges[self.NetworkGraph.EdgeNamesToIds[edge]]
            for el in self.NetworkMesh.GraphEdgeToMesh[ed]:
                if abscissa <= el.values()[0][1] and abscissa >= el.values()[0][0]:
                    return self.NetworkMesh.ElementIdsToElements[str(el.keys()[0])]                            
        except KeyError:
            node = self.NetworkGraph.Nodes[self.NetworkGraph.NodeNamesToIds[edge]]
            return self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]]
                                             
    def GetEdge(self, edge, abscissa = None):
        '''
        This method returns corresponding edge or superedge
        '''
        if edge in self.NetworkGraph.EdgeNamesToIds:
            if abscissa is not None:
                self.NetworkGraph.Edges[self.NetworkGraph.EdgeNamesToIds[edge]].edgeAbscissa = abscissa
            else:
                abscissa = None
                self.NetworkGraph.Edges[self.NetworkGraph.EdgeNamesToIds[edge]].edgeAbscissa = abscissa
            
            return self.NetworkGraph.Edges[self.NetworkGraph.EdgeNamesToIds[edge]]
        elif edge in self.NetworkGraph.SuperEdgeNamesToIds:
            
            return self.NetworkGraph.SuperEdges[self.NetworkGraph.SuperEdgeNamesToIds[edge]]
        else:
            return None
              
    def GetVariableComponents(self,variable):
        '''
        This method split expression into variables.
        '''
        parameter = self.parameterRe.findall(variable)[0][1:-1]
        try:
            element = self.elementRe.findall(variable)[0][1:-1]
            abscissa = 0.0
            timeIndex = 0
            if len(element.split(',')) > 1:
                
                splitElement = element.split(',')[1:]
                try:
                    abscissa = float(splitElement[0])
                    timeIndex= int(splitElement[0])
                   
                except:
                    el = splitElement[0]
                    par = el.split('=')[0].strip()
                    if par == 's':
                        abscissa = float(el.split('=')[1])
                    elif par == 't':
                        timeIndex = int(el.split('=')[1])
                if len(splitElement) > 2:
                    try:
                        timeIndex = int(splitElement[1])
                    except:
                        el = splitElement[1]
                        par = el.split('=')[0].strip()
                        if par == 's':
                            abscissa = float(el.split('=')[1])
                        elif par == 't':
                            timeIndex = int(el.split('=')[1])
                element = element.split(',')[0]
            edge = None
        except:
            edge = self.edgeRe.findall(variable)[0][1:-1]
            element = None
            abscissa = None
            timeIndex = 0
            if len(edge.split(',')) > 1:
                abscissa = float(edge.split(',')[1])
                edge = edge.split(',')[0]
        
        return parameter, element, abscissa, timeIndex, edge
        
    def Evaluate(self,expression):
        '''
        Evaluate(expr,{'DofMap':self.DofMap}...)
        This method evaluates provided expression and returns result.
        '''
        info = self.Info
         
        if expression in self.ExpressionCache: 
            elEvals = self.ExpressionCache[expression]['elEvals']
            lhsEvalDict = self.ExpressionCache[expression]['lhsEvalDict']
            lhsExpr = self.ExpressionCache[expression]['lhsExpr']
            try:
                lhsElement = lhsEvalDict['lhsElement']
                lhsAbscissa = lhsEvalDict['lhsAbscissa']
            except KeyError:
                lhsEdge = lhsEvalDict['lhsEdge']
            for elEvalDict in elEvals:
                try:
                    rhsElement = elEvalDict['rhsElement']
                    rhsAbscissa = elEvalDict['rhsAbscissa']
                    rhsTimeIndex = elEvalDict['rhsTimeIndex']
                except KeyError:
                    rhsEdge = elEvalDict['rhsEdge']
                exec(elEvalDict['elEval'])
            exec(lhsEvalDict['lhsEval'])
            exec(lhsExpr)
            return
        
        splitExpression = expression.split('=')
        lhs = splitExpression[0]
        rhs = splitExpression[1]
        lhsVariable = self.variableRe.findall(lhs)[0]
        lhsParameter, lhsElement, lhsAbscissa, lhsTimeIndex, lhsEdge = self.GetVariableComponents(lhsVariable)
       
        rhsVariables = self.variableRe.findall(rhs)
        elCount = 0
        elEvals = []
        
        
        if lhsEdge is None:
            if self.rhsCache.has_key(lhsElement):
                self.rhsCache[lhsElement].append(rhs)
            else:
                self.rhsCache[lhsElement] = [rhs]
            if len(self.rhsCache[lhsElement]) == 2:
                if self.rhsCache.has_key(lhsElement) and rhs == self.rhsCache[lhsElement][0]:
                    rhs = self.rhsCache[lhsElement][1]
       
        for rhsVariable in rhsVariables:
            rhsParameter, rhsElement, rhsAbscissa, rhsTimeIndex, rhsEdge = self.GetVariableComponents(rhsVariable)
            if rhsEdge is not None:
                elEvals.append({'elEval': 'edge%d = self.GetEdge(rhsEdge,rhsAbscissa)' % elCount, 'rhsEdge':rhsEdge, 'rhsAbscissa':rhsAbscissa})
                rhs = self.variableRe.sub('edge%d.Get%s(rhsAbscissa)' % (elCount,rhsParameter),rhs,1)
            else:    
                if rhsElement == '':
                    rhsParameter = self.SimulationContext.Context[rhsParameter]
                    rhs = self.variableRe.sub('%s' % (rhsParameter),rhs,1)         
                else:     
                    elEvals.append({'elEval': 'el%d = self.GetElement(rhsElement,rhsAbscissa)' % elCount, 'rhsElement':rhsElement, 'rhsAbscissa':rhsAbscissa, 'rhsTimeIndex':rhsTimeIndex})
                    rhs = self.variableRe.sub('el%d.Get%s(info,%d)' % (elCount,rhsParameter,rhsTimeIndex),rhs,1)                
            elCount += 1  
            
        if lhsEdge is None:  
            self.rhsCache[lhsElement].append(rhs) 
            
        if lhsEdge is not None:
            lhsExpr = self.variableRe.sub('lhsEdge.Set%s(%s)' % (lhsParameter,rhs),lhs,1)
            lhsEvalDict = {'lhsEval':'lhsEdge = self.GetEdge(lhsEdge,lhsAbscissa)', 'lhsEdge':lhsEdge, 'lhsAbscissa':lhsAbscissa}     
        else:
            if lhsElement == '':            
                lhsExpr = self.variableRe.sub('self.SimulationContext.Context[%s]=%s' % ("'"+lhsParameter+"'",rhs),lhs,1)
                lhsEvalDict = {'lhsEval':'', 'lhsElement':None, 'lhsAbscissa':None}
            else:
                lhsExpr = self.variableRe.sub('lhsEl.Set%s(%s, info)' % (lhsParameter,rhs),lhs,1)
                lhsEvalDict = {'lhsEval':'lhsEl = self.GetElement(lhsElement,lhsAbscissa)', 'lhsElement':lhsElement, 'lhsAbscissa':lhsAbscissa}
        for elEvalDict in elEvals:         
            try:
                rhsEdge = elEvalDict['rhsEdge']
                rhsAbscissa = elEvalDict['rhsAbscissa']              
            except KeyError:
                rhsElement = elEvalDict['rhsElement']
                rhsAbscissa = elEvalDict['rhsAbscissa']           
            exec(elEvalDict['elEval'])          
        exec(lhsEvalDict['lhsEval'])
        exec(lhsExpr)  
        self.ExpressionCache[expression] = {'elEvals':elEvals, 'lhsEvalDict':lhsEvalDict, 'lhsExpr':lhsExpr}