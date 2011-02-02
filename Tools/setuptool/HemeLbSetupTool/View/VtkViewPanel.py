import wx

from vtk import vtkRenderer, vtkPolyDataMapper, vtkActor, vtkModifiedBSPTree

from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

from Delegator import Delegator, ParentDelegator
from Layout import H, V, StretchSpacer
from Placers import SeedPlacer

import pdb

class VtkViewPanel(wx.Panel, ParentDelegator):
    def __init__(self, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        ParentDelegator.__init__(self)
        
        self.exports.debugButton = wx.Button(self, label="DEBUG")
        
        self.exports.resetViewButton = wx.Button(self, label="Reset")
        
        self.exports.xButton = wx.Button(self, label="X")
        self.exports.yButton = wx.Button(self, label="Y")
        self.exports.zButton = wx.Button(self, label="Z")
        self.exports.rwi = RWI(self)

        layout = V((V((H(StretchSpacer(), self.exports.resetViewButton,
                         self.exports.xButton, self.exports.yButton,
                         self.exports.zButton, StretchSpacer()), 0, wx.EXPAND),
                      (self.rwi, 1, wx.EXPAND)), 1, wx.EXPAND))
        self.SetSizer(layout.create())

        return
    
class RWI(wxVTKRenderWindowInteractor, ParentDelegator):
    """Set up for the VTK window (on the RHS of the window).
    """
    def __init__(self, parent, *args, **kwargs):
        wxVTKRenderWindowInteractor.__init__(self, parent, wx.ID_ANY,
                          *args, **kwargs)
        ParentDelegator.__init__(self)
        self.AddObserver("ExitEvent", lambda o,e,f=parent: f.Close())
        
        self.exports.renderer = vtkRenderer()
        self.GetRenderWindow().AddRenderer(self.renderer)

        # Set the up direction and default to trackball mode for view control
        self.renderer.GetActiveCamera().SetViewUp(0.,0.,1.)
        self.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
        
        self.exports.pipeline = Pipeline(self)
        return

    def __getattr__(self, attr):
        """Since wxVTKRenderWindowInteractor has a __getattr__, we
        must be careful here.
        """
        try:
            return wxVTKRenderWindowInteractor.__getattr__(self, attr)
        except AttributeError:
            return ParentDelegator.__getattr__(self, attr)

    def ResetView(self):
        """Reset the view on the current scene.
        """
        self.renderer.ResetCamera()
        cam = self.renderer.GetActiveCamera()
        focus = cam.GetFocalPoint()
        dist = cam.GetDistance()
        cam.SetPosition(focus[0]+dist,focus[1], focus[2])
        cam.SetViewUp(0,0,1)
        return
    
    pass

class Pipeline(Delegator):
    """Represent the VTK pipeline.
    """
    def __init__(self, delegate):
        Delegator.__init__(self, delegate)
        self.stlMapper = vtkPolyDataMapper()
        self.stlActor = vtkActor()
        self.stlActor.SetMapper(self.stlMapper)
        
        self.seeder = SeedPlacer(self, self.stlActor)
         
        self.locator = vtkModifiedBSPTree()
        return
    
    def Start(self, surfaceOutputPort):
        self.stlMapper.SetInputConnection(surfaceOutputPort)
        self.locator.SetDataSet(surfaceOutputPort.GetProducer().GetOutput())
        self.locator.BuildLocator()
        return

    def IsActorAdded(self, actor):
        """Return whether the supplied argument is in the renderer's
        list of actors.
        """
        aList = self.renderer.GetActors()
        iterator = aList.NewIterator()
        
        while not iterator.IsDoneWithTraversal():
            a = iterator.GetCurrentObject()
            if a is actor:
                return True
            iterator.GoToNextItem()
            continue
        return False
    
    def Show(self):
        if not self.IsActorAdded(self.stlActor):
            self.renderer.AddActor(self.stlActor)
            self.rwi.Update()
            pass
        return
    
    def PlaceSeed(self):
        self.seeder.PlaceOn()
        return

    def SetSeedPoint(self, pos):
        return self.seeder.SetSeedPoint(pos)
    
    pass
