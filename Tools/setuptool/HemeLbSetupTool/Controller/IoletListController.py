import wx

from HemeLbSetupTool.Bindings.ListController import ListController
from HemeLbSetupTool.Controller.VectorController import VectorController
from HemeLbSetupTool.Controller.IoletController import IoletController
from HemeLbSetupTool.Model.Iolets import Inlet, Outlet

class IoletListController(ListController):
    def __init__(self, delegate):
        ListController.__init__(self, delegate, SelectionControllerClass=IoletController)
        self.nInlets = 0
        self.nOutlets = 0
        return
    
    def AddInlet(self):
        self.nInlets += 1
        self.delegate.append(Inlet(Name='Inlet%d' % (self.nInlets)))
        self.SelectedIndex = len(self.delegate)-1
        return
    
    def AddOutlet(self):
        self.nOutlets += 1
        self.delegate.append(Outlet(Name='Outlet%d' % (self.nOutlets)))
        self.SelectedIndex = len(self.delegate)-1
        return
    
    pass
