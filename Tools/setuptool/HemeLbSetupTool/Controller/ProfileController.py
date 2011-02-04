import wx

from HemeLbSetupTool.Bindings.ObjectController import ObjectController
from HemeLbSetupTool.Bindings.Translators import QuickTranslator
from HemeLbSetupTool.Bindings.VtkObject import HasVtkObjectKeys
from HemeLbSetupTool.Controller.IoletListController import IoletListController, HasIoletListKeys
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
    
            
    # def ResetView(self, event):
    #     """Tell the view to reset.
    #     """
    #     # Trigger the custom RenderWindowInteractor's Reset Method
    #     self.view.rwi.ResetView()
    #     return

    # # STL file selection events
    # @EventHandler
    def ChooseStl(self, ignored=None):
        dialog = wx.FileDialog(None, style=wx.FD_OPEN, wildcard='*.stl')
        
        if dialog.ShowModal() == wx.ID_OK:
            self.delegate.StlFile = dialog.GetPath()
            
            # self.view.pipeline.Show()
            # self.view.rwi.ResetView()
            pass
        
        dialog.Destroy()
        return

    # IOLets events


    # Seed placement events
    # @EventHandler
    # def OnViewPlaceSeed(self, event):
    #     if self._placeMode == 'off':
    #         self._placeMode = 'seed'
    #         self.view.pipeline.seeder.PlaceOn()
    #         self.view.seedPlaceButton.SetLabel('Cancel')
    #     elif self._placeMode == 'seed':
    #         self.EndSeedMode()
    #     return
    
    # def EndSeedMode(self):
    #     self._placeMode = 'off'
    #     self.view.pipeline.seeder.PlaceOff()
    #     self.view.seedPlaceButton.SetLabel('Place')
    #     return
    
    # @MessageHandler
    # def OnViewSeedPointSent(self, msg):
    #     self.EndSeedMode()
    #     self.model.SeedPoint = msg.data
    #     return
    
    # @EventHandler
    # def OnViewSeedPointEdited(self, event):
    #     seedVec = self.view.seedVector
    #     pos = [seedVec.x.GetValue(),
    #            seedVec.y.GetValue(),
    #            seedVec.z.GetValue()]
    #     self.model.SeedPoint = pos
    #     return
    
    # Output file events
    def ChooseOutputConfigFile(self, ignored=None):
        dialog = wx.FileDialog(None,
                               style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT,
                               wildcard='*.dat')

        if dialog.ShowModal() == wx.ID_OK:
            self.delegate.OutputConfigFile = dialog.GetPath()
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
        
    pass
