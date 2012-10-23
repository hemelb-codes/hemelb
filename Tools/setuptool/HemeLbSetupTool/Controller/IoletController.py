# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

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
