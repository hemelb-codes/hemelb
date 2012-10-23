# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import wx

from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

from HemeLbSetupTool.View.Layout import H, V, StretchSpacer
from HemeLbSetupTool.Bindings.Bindings import WxActionBinding

import pdb

class VtkViewPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        self.resetViewButton = wx.Button(self, label="Fit")
        controller.BindAction('ResetView',
                              WxActionBinding(self.resetViewButton, wx.EVT_BUTTON))
        
        self.xButton = wx.Button(self, label="X")
        controller.BindAction('SetViewX',
                              WxActionBinding(self.xButton, wx.EVT_BUTTON))
        self.yButton = wx.Button(self, label="Y")
        controller.BindAction('SetViewY',
                              WxActionBinding(self.yButton, wx.EVT_BUTTON))
        self.zButton = wx.Button(self, label="Z")
        controller.BindAction('SetViewZ',
                              WxActionBinding(self.zButton, wx.EVT_BUTTON))

        self.rwi = RWI(controller, self)
        
        layout = V((V((H(StretchSpacer(), self.resetViewButton,
                         self.xButton, self.yButton,
                         self.zButton, StretchSpacer()), 0, wx.EXPAND),
                      (self.rwi, 1, wx.EXPAND)), 1, wx.EXPAND))
        self.SetSizer(layout.create())

        return
    pass

class RWI(wxVTKRenderWindowInteractor):
    """Set up for the VTK window (on the RHS of the window).
    """
    def __init__(self, controller, parent, *args, **kwargs):
        wxVTKRenderWindowInteractor.__init__(self, parent, wx.ID_ANY,
                          *args, **kwargs)
        self.AddObserver("ExitEvent", lambda o,e,f=parent: f.Close())
        
        self.GetRenderWindow().AddRenderer(controller.GetValueForKey('Renderer'))
        
        # Set the up direction and default to trackball mode for view control
        # self.renderer.GetActiveCamera().SetViewUp(0.,0.,1.)
        # self.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
        
        controller.GetValueForKey('SetInteractor')(self)
        
        return
    
    pass

