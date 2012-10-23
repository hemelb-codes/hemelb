# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import numpy as N

from vtk import vtkPlaneWidget, vtkPolyDataMapper, vtkActor

from ..Util.Observer import Observable, ObservableList, NotifyOptions

import pdb
class PlacedIoletList(ObservableList):
    def __init__(self, *args, **kwargs):
        ObservableList.__init__(self, *args, **kwargs)

        # Add observers to add/remove observers for items. Note that
        # for removal we must trigger the action BEFORE the removal
        # happens, so that we can get the item to remove the observer.
        self.AddObserver("@INSERTION", self.HandleInsertion)
        self.AddObserver("@REMOVAL", self.HandlePreRemoval,
                         options=NotifyOptions(BEFORE_CHANGE=True,
                                               AFTER_CHANGE=False))
        return
    
    def SetItemEnabledChangeHandler(self, handler):
        self._ItemEnabledChangeHandler = handler
        return
    
    def SetInteractor(self, iact):
        self.Interactor = iact
        return
    
    def HandleInsertion(self, change):
        self[change.index].AddObserver('Enabled', self._ItemEnabledChangeHandler)
        self[change.index].widget.SetInteractor(self.Interactor)
        self[change.index].Enabled = True
        return
    
    def HandlePreRemoval(self, change):
        self[change.index].Enabled = False
        self[change.index].RemoveObserver('Enabled', self._ItemEnabledChangeHandler)
        return
    pass

class PlacedIolet(Observable):
    
    def __init__(self):
        # Working arrays, constructed once for speed.
        self._c = N.zeros(3)
        self._o = N.zeros(3)
        self._p1 = N.zeros(3)
        self._p2 = N.zeros(3)
        self._n = N.zeros(3)

        self.widget = vtkPlaneWidget()
        self.widget.KeyPressActivationOff()
        self.widget.SetRepresentationToOutline()
        # planeWidget.PlaceWidget(clickPos)
        
        # self.representation = vtkPolyData()
        # This is effectively a copy and is guaranteed to be up to
        # date when InteractionEvent or EndInteraction events are
        # invoked
        # self.widget.AddObserver("InteractionEvent", self._SyncRepresentation)
        self.representation = self.widget.GetPolyDataAlgorithm()
        
        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.representation.GetOutputPort())
        
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(self.colour)
        
        # Keep a cached copy of the radius etc to minimise the number of notifications we must send
        self._lastRadius = self.Radius
        self._lastCentre = self.Centre
        self._lastNormal = self.Normal
        
        # We can only enable after self.SetInteractor() has been called.
        self._Enabled = False
        
        
        self.widget.AddObserver("EndInteractionEvent", self.HandleInteraction)
        # self.AddObserver('Enabled', self.EnabledSet)
        return
    
    def HandleInteraction(self, obj, evt):
        """This handler is called from VTK when the widget has been
         interacted with by the user. We check each of the 
         properties against the cached values and if changed, notify
         any observers of this.
         """
        if self.Centre != self._lastCentre:
            self.DidChangeValueForKey('Centre')
            self._lastCentre = self.Centre
            pass
        if self.Normal != self._lastNormal:
            self.DidChangeValueForKey('Normal')
            self._lastNormal = self.Normal
            pass
        if self.Radius != self._lastRadius:
            self.DidChangeValueForKey('Radius')
            self._lastRadius = self.Radius
            pass
        return
    
    def GetEnabled(self):
        return self._Enabled
    def SetEnabled(self, enabled):
        if self.widget.GetInteractor() is None:
            return
        self._Enabled = enabled
        if enabled:
            self.widget.On()
        else:
            self.widget.Off()
            pass
        return
    Enabled = property(GetEnabled, SetEnabled)
    
    def _SyncRepresentation(self, obj, evt):
        self.widget.GetPolyData(self.representation)
        return
    
    def SetCentre(self, centre):
        self.widget.SetCenter(centre)
        # Force the plane to be updated
#        self.widget.InvokeEvent("InteractionEvent")
        return
    def GetCentre(self):
        return self.widget.GetCenter()
    Centre = property(GetCentre, SetCentre)
    
    def SetNormal(self, normal):
        self.widget.SetNormal(normal)
        # Force the plane to be updated
#        self.widget.InvokeEvent("InteractionEvent")
        return
    def GetNormal(self):
        return self.widget.GetNormal()
    Normal = property(GetNormal, SetNormal)
    
    def SetRadius(self, radius):
        # Get into numpy vectors
        self.widget.GetCenter(self._c)
        self.widget.GetOrigin(self._o)
        self.widget.GetPoint1(self._p1)
        self.widget.GetPoint2(self._p2)
        # Make corners relative to centre
        self._o -= self._c
        self._p1 -= self._c
        self._p2 -= self._c
        # Compute norms
        oNorm = N.dot(self._o, self._o)
        p1Norm = N.dot(self._p1, self._p1)
        p2Norm = N.dot(self._p2, self._p2)
        # Scale
        self._o *= radius * N.sqrt(2. / oNorm)
        self._p1 *= radius * N.sqrt(2. / p1Norm)
        self._p2 *= radius * N.sqrt(2. / p2Norm)
        # Add the centre back on
        self._o += self._c
        self._p1 += self._c
        self._p2 += self._c
        # Set
        self.widget.SetOrigin(self._o)
        self.widget.SetPoint1(self._p1)
        self.widget.SetPoint2(self._p2)
        # Force the plane to be updated
#        self.widget.InvokeEvent("InteractionEvent")
        return
    
    def GetRadius(self):
        # Get into numpy vectors
        self.widget.GetOrigin(self._o)
        self.widget.GetPoint1(self._p1)
        
        self._p1 -= self._o
        
        return 0.5 * N.sqrt(N.dot(self._p1, self._p1))
    
    Radius = property(GetRadius, SetRadius)
    
    def HandleWidgetSizeChange(self, change):
#        self.widget.InvokeEvent('LeftButtonPressEvent')
#        self.widget.InvokeEvent('LeftButtonReleaseEvent')
        return
    
    pass

class PlacedInlet(PlacedIolet):    
    def __init__(self):
        self.colour = (0., 1., 0.)
        PlacedIolet.__init__(self)
        return
    pass

class PlacedOutlet(PlacedIolet):    
    def __init__(self):
        self.colour = (1., 0., 0.)
        PlacedIolet.__init__(self)
        return
    pass

