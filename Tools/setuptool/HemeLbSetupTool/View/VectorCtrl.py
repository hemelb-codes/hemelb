import wx
# from wx.lib.masked import NumCtrl, EVT_NUM

from HemeLbSetupTool.Bindings.Translators import FloatTranslator, NoneToValueTranslator
from HemeLbSetupTool.Bindings.Mappers import WxWidgetMapper, Mapper
from HemeLbSetupTool.View.Layout import H

class VectorCtrl(wx.Panel):
    """Simple container of three TextCtrl's for a vector quantity.
    """
    
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        
        # self.x = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, integerWidth=3, fractionWidth=3)
        # self.y = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, integerWidth=3, fractionWidth=3)
        # self.z = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, integerWidth=3, fractionWidth=3)
        self.x = wx.TextCtrl(self)
        self.y = wx.TextCtrl(self)
        self.z = wx.TextCtrl(self)
        sizer = H((self.x, 0, wx.EXPAND),
                  (self.y, 0, wx.EXPAND),
                  (self.z, 0, wx.EXPAND)).create()
        self.SetSizer(sizer)
        
        return
    
    pass

class VectorMapper(WxWidgetMapper):
    """Widget mapper for VectorCtrls.
    """
    def __init__(self, widget, key, event,
                 translator=NoneToValueTranslator(float('nan'),
                                                  inner=FloatTranslator())
                 ):
        # We want to skip the WxWidgetMapper's init for now as the
        # VectorCtrl typically won't have the required getters and
        # setters. On binding, this one mapper is turned into three
        # standard mappers anyway.
        Mapper.__init__(self, translator=translator)
        
        self.widget = widget
        self.key = key
        self.event = event
        
        return
    
    def CreateSubMapper(self, component):
        return WxWidgetMapper(getattr(self.widget, component),
                              self.key, self.event,
                              translator=self.translator)
    
    pass
