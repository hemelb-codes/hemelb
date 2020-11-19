# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import weakref

from HemeLbSetupTool.Bindings.ObjectController import ObjectController
from HemeLbSetupTool.Controller.VectorController import HasVectorKeys
import pdb

class IoletController(HasVectorKeys, ObjectController):
    """Controller for inlets and outlets.
    
    The constructor
    """
    _instanceDict = weakref.WeakValueDictionary()
    @classmethod
    def New(cls, delegate):
        try:
            con = cls._instanceDict[delegate]
        except KeyError:
            con = cls._instanceDict[delegate] = cls(delegate)
            pass
        return con
    
    def __init__(self, delegate):
        ObjectController.__init__(self, delegate)
        
        self.DefineVectorKey("Centre")
        self.DefineVectorKey("Normal")
        self.DefineVectorKey("Pressure")
        return
    pass
