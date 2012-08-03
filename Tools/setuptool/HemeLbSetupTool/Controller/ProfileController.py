# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import wx

from HemeLbSetupTool.Bindings.ObjectController import ObjectController
from HemeLbSetupTool.Bindings.Translators import QuickTranslator
from HemeLbSetupTool.Bindings.VtkObject import HasVtkObjectKeys
from HemeLbSetupTool.Controller.IoletListController import HasIoletListKeys
from HemeLbSetupTool.Controller.VectorController import HasVectorKeys

import pdb

class ProfileController(HasIoletListKeys, HasVectorKeys, HasVtkObjectKeys, ObjectController):
    """Proxy for HemeLbSetupTool.Mode.Profile objects.
    """
    def __init__(self, delegate):
        ObjectController.__init__(self, delegate)
        self.DefineVectorKey("SeedPoint")
        self.DefineIoletListKey("Iolets")
        self.DefineVtkObjectKey("StlReader")
        self.DefineVtkObjectKey('SideLengthCalculator')
        
        def trans(state):
            if state:
                return wx.WHITE
            return wx.Colour(255,128,128)
        
        self.validColourer = QuickTranslator(trans, lambda x: False) 
        return
    
    def Debug(self, ignored=None):
        """Drop into the debugger.
        """
        pdb.set_trace()
        return
    
    def ChooseStl(self, ignored=None):
        dialog = wx.FileDialog(None, style=wx.FD_OPEN, wildcard='*.stl')
        
        if dialog.ShowModal() == wx.ID_OK:
            self.delegate.StlFile = dialog.GetPath()
            pass
        
        dialog.Destroy()
        return
    
    def ChooseOutputGeometryFile(self, ignored=None):
        dialog = wx.FileDialog(None,
                               style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT,
                               wildcard='*.gmy')

        if dialog.ShowModal() == wx.ID_OK:
            self.delegate.OutputGeometryFile = dialog.GetPath()
            pass
        
        dialog.Destroy()
        return
    
    def ChooseOutputXmlFile(self, ignored=None):
        dialog = wx.FileDialog(None,
                               style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT,
                               wildcard='*.xml')
        
        if dialog.ShowModal() == wx.ID_OK:
            self.delegate.OutputXmlFile = dialog.GetPath()
            pass
        
        dialog.Destroy()
        return
        
    def ChooseSaveFile(self, ignored=None):
        dialog = wx.FileDialog(None,
                               style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT,
                               wildcard='*.pro')
        
        if dialog.ShowModal() == wx.ID_OK:
            try:
                self.delegate.Save(dialog.GetPath())
            except IOError as err:
                errDiag = wx.MessageDialog(None,
                                           'Cannot write profile file.\nMessage: ' + str(err),
                                           style=wx.OK | wx.ICON_ERROR)
                errDiag.ShowModal()
            pass
        
        dialog.Destroy()
        return
    
    def LoadFromFile(self, ignored=None):
        dialog = wx.FileDialog(None,
                               style=wx.FD_OPEN,
                               wildcard='*.pro')
        
        if dialog.ShowModal() == wx.ID_OK:
            self.delegate.LoadFromFile(dialog.GetPath())
            pass
        
        dialog.Destroy()
        return
    pass
