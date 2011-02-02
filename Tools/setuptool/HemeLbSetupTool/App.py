import wx

from Model import Profile
from View import WxView
from Controller import Controller

class SetupTool(wx.App):
    def OnInit(self):
        # Model
        profile = Profile()
        # View
        view = WxView()
        # Controller
        controller = Controller(profile, view)

        self.SetTopWindow(controller.InitWindow())
        return True
    
    pass

