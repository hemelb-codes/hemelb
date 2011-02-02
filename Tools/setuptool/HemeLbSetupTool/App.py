import wx

from HemeLbSetupTool.Model.Profile import Profile
from HemeLbSetupTool.Controller.ProfileController import ProfileController
from HemeLbSetupTool.View.MainWindow import MainWindow

class SetupTool(wx.App):
    def OnInit(self):
        # Model
        self.profile = Profile()
        # Controller
        self.controller = ProfileController(self.profile)
        
        # View
        self.view = MainWindow(self.controller)

        self.SetTopWindow(self.view)
        return True
    
    pass

