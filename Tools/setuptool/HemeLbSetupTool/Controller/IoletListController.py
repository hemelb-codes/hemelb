# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from HemeLbSetupTool.Bindings.ListController import ListController, HasListKeys
from HemeLbSetupTool.Controller.IoletController import IoletController
from HemeLbSetupTool.Model.Iolets import PressureInlet, PressureOutlet, \
    VelocityInlet, VelocityOutlet
import pdb
class IoletListController(ListController):
    def __init__(self, delegate):
        ListController.__init__(self, delegate, SelectionControllerClass=IoletController.New)
        self.nInlets = 0
        self.nOutlets = 0
        self.InletTypeToAdd = 'Pressure'
        self.OutletTypeToAdd = 'Pressure'
        return
    
    def AddInlet(self):
        self.nInlets += 1
        if self.InletTypeToAdd == 'Pressure':
            InletClass = PressureInlet;
        elif self.InletTypeToAdd == 'Velocity':
            InletClass = VelocityInlet;
        else:
            raise Exception('Unknown inlet type!')
        self.append(InletClass(Name='Inlet%d' % (self.nInlets)))
        return
    
    def AddOutlet(self):
        self.nOutlets += 1
        if self.OutletTypeToAdd == 'Pressure':
            OutletClass = PressureOutlet;
        elif self.OutletTypeToAdd == 'Velocity':
            OutletClass = VelocityOutlet;
        else:
            raise Exception('Unknown outlet type!')
        self.append(OutletClass(Name='Outlet%d' % (self.nOutlets)))
        return

    def RemoveIolet(self):
        # TODO: Make sure PlacedIolet's widget + actor are removed from the scene. 
        if self.SelectedIndex is None:
            return
        del self[self.SelectedIndex]
        return
    
    pass

class HasIoletListKeys(HasListKeys):
    """Mixin for ObjectController subclasses with IoletList keys.
    """
    BindFunctionDispatchTable = ((IoletListController, 'BindList'),)
    
    def DefineIoletListKey(self, name):
        """Typically used in the subclass __init__ method to easily
        mark a key as being a List and hence needing a ListController
        to manage it.
        """
        setattr(self, name,
                IoletListController(getattr(self.delegate, name))
                )
        return
    
    pass
