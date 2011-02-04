import wx

from HemeLbSetupTool.Model.Profile import Profile
from HemeLbSetupTool.Model.Pipeline import Pipeline

from HemeLbSetupTool.Controller.ProfileController import ProfileController
from HemeLbSetupTool.Controller.PipelineController import PipelineController

from HemeLbSetupTool.View.MainWindow import MainWindow

class SetupTool(wx.App):
    def OnInit(self):
        # Model
        self.profile = Profile()
        self.pipeline = Pipeline()
        
        # Controller
        self.controller = ProfileController(self.profile)
        self.controller.Pipeline = PipelineController(self.pipeline, self.controller)
        
        # View
        self.view = MainWindow(self.controller)

        self.SetTopWindow(self.view)
        return True
    
    pass

