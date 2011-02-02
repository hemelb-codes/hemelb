import wx

from HemeLbSetupTool.Bindings.ObjectController import ObjectController
from HemeLbSetupTool.Controller.VectorController import HasVectorKeys
import pdb

class IoletController(HasVectorKeys, ObjectController):
    def __init__(self, delegate):
        ObjectController.__init__(self, delegate)
        
        self.DefineVectorKey("Centre")
        self.DefineVectorKey("Normal")
        self.DefineVectorKey("Pressure")
        return
    pass
