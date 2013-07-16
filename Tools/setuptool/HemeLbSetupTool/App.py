# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

try:
    import wxversion
    wxversion.select('2.8')
except wxversion.VersionError:
    wxversion.select('2.9')
except ImportError:
    pass
import wx

from HemeLbSetupTool.Model.Profile import Profile
from HemeLbSetupTool.Model.Pipeline import Pipeline

from HemeLbSetupTool.Controller.ProfileController import ProfileController
from HemeLbSetupTool.Controller.PipelineController import PipelineController

from HemeLbSetupTool.View.MainWindow import MainWindow

class SetupTool(wx.App):
    def __init__(self, args={}, profile=None, **kwargs):
        self.cmdLineArgs = args
        self.cmdLineProfileFile = profile
        
        wx.App.__init__(self, **kwargs)
        return
    
    def OnInit(self):
        # Model
        self.profile = Profile()
        
        self.pipeline = Pipeline()
        
        # Controller
        self.controller = ProfileController(self.profile)
        self.controller.Pipeline = PipelineController(self.pipeline, self.controller)
        
        # View
        self.view = MainWindow(self.controller)

        if self.cmdLineProfileFile is not None:
            # Load the profile
            self.profile.LoadFromFile(self.cmdLineProfileFile)
            pass
        
        # override any keys that have been set on cmdline.
        self.UpdateAttributesBasedOnCmdLineArgs(self.cmdLineArgs)
        
        self.SetTopWindow(self.view)
        return True
    
    pass

