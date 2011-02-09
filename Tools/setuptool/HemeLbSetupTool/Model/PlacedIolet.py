import numpy as N

from vtk import vtkPlaneWidget, vtkPolyData, vtkPolyDataMapper, \
     vtkActor

from ..Util.Observer import Observable, ObservableList, NotifyOptions

from .Iolets import Inlet, Outlet

class PlacedIoletList(ObservableList):
    def __init__(self, *args, **kwargs):
        ObservableList.__init__(self, *args, **kwargs)

        # Add observers to add/remove observers for items. Note that
        # for removal we must trigger the action BEFORE the removal
        # happens, so that we can get the item to remove the observer.
        self.AddObserver("@INSERTION", self.HandleInsertion)
        self.AddObserver("@REMOVAL", self.HandlePreRemoval,
                         options=NotifyOptions(BEFORE_CHANGE=True))
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
        return
    
    def HandlePreRemoval(self, change):
        self[change.index].RemoveObserver('Enabled', self._ItemEnabledChangeHandler)
        return
    pass

class PlacedIolet(Observable):
    
    def __init__(self, iolet):
        self.iolet = iolet
        if isinstance(iolet, Inlet):
            self.colour = (0.,1.,0.)
        elif isinstance(iolet, Outlet):
            self.colour = (1.,0.,0.)
        else:
            raise ValueError("IOlet must be either an inlet or outlet!")
        
        # Working arrays, constructed once for speed.
        self._c = N.zeros(3)
        self._p1 = N.zeros(3)
        self._p2 = N.zeros(3)

        self.widget = vtkPlaneWidget()
        self.widget.SetRepresentationToOutline()
        # planeWidget.PlaceWidget(clickPos)
        
        self.representation = vtkPolyData()
        # This is effectively a copy and is guaranteed to be up to
        # date when InteractionEvent or EndInteraction events are
        # invoked
        self.widget.AddObserver("InteractionEvent", self._SyncRepresentation)
        self.widget.GetPolyData(self.representation)

        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInput(self.representation)
        
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(self.colour)
        
        self.Enabled = False
        self.AddObserver('Enabled', self.EnabledSet)

        return
    
    def EnabledSet(self, change):
        if self.Enabled:
            self.widget.On()
        else:
            self.widget.Off()
            pass
        return
    
    def _SyncRepresentation(self, obj, evt):
        self.widget.GetPolyData(self.representation)
        return
    
    def SetCentre(self, centre):
        self.widget.SetCenter(centre)
        # Force the plane to be updated
        self.widget.InvokeEvent("InteractionEvent")
        return
    def GetCentre(self):
        return self.widget.GetCenter()
    Centre = property(GetCentre, SetCentre)
    
    def SetNormal(self, normal):
        # self.widget.SetNormal(normal)
        # # Force the plane to be updated
        # self.widget.InvokeEvent("InteractionEvent")
        return
    def GetNormal(self):
        return self.widget.GetNormal()
    Normal = property(GetNormal, SetNormal)
    
    def SetRadius(self, radius):
        # if not N.isfinite(radius):
        #     self.Enabled = False
        #     return
        
        # # Get into numpy vectors
        # self.widget.GetCenter(self._c)
        # self.widget.GetPoint1(self._p1)
        # self.widget.GetPoint2(self._p2)

        # # p1/2 now are relative to centre
        # self._p1 -= self._c
        # self._p2 -= self._c
        # # Rescale so that they have length == radius
        # self._p1 *= radius / N.sqrt(N.dot(self._p1, self._p1))
        # self._p2 *= radius / N.sqrt(N.dot(self._p2, self._p2))
        
        # # Set
        # self.widget.SetPoint1(self._p1)
        # self.widget.SetPoint1(self._p2)
        # # Force the plane to be updated
        # self.widget.InvokeEvent("InteractionEvent")
        return
    def GetRadius(self):
        # Get into numpy vectors
        self.widget.GetCenter(self._c)
        self.widget.GetPoint1(self._p1)
        self.widget.GetPoint2(self._p2)

        # p1/2 now are relative to centre
        self._p1 -= self._c
        self._p2 -= self._c
        # Square
        self._p1 **= 2
        self._p2 **= 2
        r1 = N.sqrt(self._p1.sum())
        r2 = N.sqrt(self._p1.sum())
        return 0.5*(r1+r2)
    Radius = property(GetRadius, SetRadius)
    
    pass
