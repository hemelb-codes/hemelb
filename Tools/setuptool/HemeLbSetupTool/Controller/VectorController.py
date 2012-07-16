# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from HemeLbSetupTool.Bindings.ObjectController import ObjectController
# from HemeLbSetupTool.View.VectorCtrl import VectorCtrlMapper

class VectorController(ObjectController):
    """Controller for HemeLbSetupTool.Model.Vector objects.
    """
    def __init__(self, delegate):
        ObjectController.__init__(self, delegate)
        return
    
    
    pass

class HasVectorKeys(object):
    """Mixin for ObjectController subclasses with Vector keys.
    """
    BindFunctionDispatchTable = ((VectorController, 'BindVector'),)
    
    def BindVector(self, top, key, mapper):
        """Each component of the vector should be appropriately bound.
        """
        for coord in ('x', 'y', 'z'):
            # Get the key path to the component and bind the part of
            # the VectorCtrl to that path
            top.BindValue(key + '.' + coord,
                          mapper.CreateSubMapper(coord))
            continue
        return

    def DefineVectorKey(self, name):
        """Typically used in the subclass __init__ method to easily
        mark a key as being a Vector and hence needing a
        VectorController to manage it.
        """
        setattr(self, name,
                VectorController(getattr(self.delegate, name))
                )
        return
    
    pass

