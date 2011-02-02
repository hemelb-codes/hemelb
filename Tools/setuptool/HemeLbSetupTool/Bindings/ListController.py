import wx

from HemeLbSetupTool.Util.Observer import ObservableList

from HemeLbSetupTool.Bindings.ObjectController import ObjectController
from HemeLbSetupTool.Bindings.EmptySelection import EmptySelection, isNone
from HemeLbSetupTool.Bindings.Translators import Translator
from HemeLbSetupTool.Bindings.Mappers import Mapper

class ListController(ObjectController):
    def __init__(self, delegate, SelectionControllerClass=ObjectController):
        assert isinstance(delegate, ObservableList)
        ObjectController.__init__(self, delegate)
        self.SelectionControllerClass = SelectionControllerClass
        self.SelectedIndex = None
        self.AddDependency('Selection', 'SelectedIndex')
        return
    
    def HandleInsertion(self, change):
        return

    
    @property
    def Selection(self):
        if self.SelectedIndex is None:
            return self.SelectionControllerClass(EmptySelection)
        return self.SelectionControllerClass(self.delegate[self.SelectedIndex])
    
    pass

class ListCtrlSelectionMapper(Mapper, Translator):
    def __init__(self, widget, inner=None):
        Mapper.__init__(self, translator=self)
        Translator.__init__(self, inner)
        
        self.widget = widget
        return

    def TranslateStage(self, val):
        if isNone(val):
            return -1
        if val < -1 or val >= self.widget.GetItemCount():
            raise IndexError('Index "%d" out of range for %s' % (val, str(widget)))
        
        return val
    
    def UntranslateStage(self, val):
        if val == -1:
            return None
        return val
    
    def Observe(self):
        self.widget.Bind(wx.EVT_LIST_ITEM_SELECTED, self.HandleUpdate)
        self.widget.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.HandleUpdate)
        return
    
    def Unobserve(self):
        self.widget.Unbind(wx.EVT_LIST_ITEM_SELECTED)
        self.widget.Unbind(wx.EVT_LIST_ITEM_DESELECTED)
        return
    
    def _Get(self):
        return self.widget.GetFirstSelected()
    
    def _Set(self, ind):
        self.Unobserve()
        try:
            prevSelected = self._Get()
            if ind != prevSelected:
                self.widget.SetItemState(newInd, wx.LIST_STATE_SELECTED, wx.LIST_STATE_SELECTED)
                self.widget.SetItemState(prevSelected, 0, wx.LIST_STATE_SELECTED)
                pass
            
            
        finally:
            self.Observe()
            pass

        return
    
