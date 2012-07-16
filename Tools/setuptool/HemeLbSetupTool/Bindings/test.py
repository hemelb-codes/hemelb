# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import wx

from HemeLbSetupTool.Util.Observer import Observable

from HemeLbSetupTool.Bindings.ObjectController import ObjectController
from HemeLbSetupTool.Bindings.WxMappers import WxWidgetMapper
from HemeLbSetupTool.Bindings.Translators import QuickTranslator
from HemeLbSetupTool.Bindings.Bindings import WxActionBinding

import pdb

class Model(Observable):
    def __init__(self, tKelvin):
        self.tKelvin = float(tKelvin)
        self.tKelvinDefault = self.tKelvin

        self.AddObserver('tKelvin', self.change)
        return

    def change(self, ignored):
        print 'T = %f K' % self.tKelvin
        return


    def Reset(self):
        self.tKelvin = self.tKelvinDefault
        return
    
    pass

class MyController(ObjectController):
    def Debug(self):
        pdb.set_trace()
        return
    pass

class View(object):
    def __init__(self, controller):
        self.controller = controller

        self.frame = wx.Frame(None, -1, "Hello from wxPython")
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.tCelciusCtrl = wx.TextCtrl(self.frame)
        cLabel = wx.StaticText(self.frame, label='Celcius')
        cSizer = wx.BoxSizer(wx.HORIZONTAL)
        cSizer.Add(cLabel, 0, wx.EXPAND)
        cSizer.Add(self.tCelciusCtrl, 0, wx.EXPAND)
        sizer.Add(cSizer, 0, wx.EXPAND)
        
        self.controller.BindValue(
            'tKelvin',
            WxWidgetMapper(
                self.tCelciusCtrl, 'Value', wx.EVT_TEXT,
                translator=QuickTranslator(lambda k: str(k-273.),
                                           lambda c: float(c)+273.)
                )
            )
        
        self.tFarenheitCtrl = wx.TextCtrl(self.frame)
        fLabel = wx.StaticText(self.frame, label='Farenheit')
        fSizer = wx.BoxSizer(wx.HORIZONTAL)
        fSizer.Add(fLabel, 0, wx.EXPAND)
        fSizer.Add(self.tFarenheitCtrl, 0, wx.EXPAND)
        sizer.Add(fSizer, 0, wx.EXPAND)
        self.controller.BindValue(
            'tKelvin',
            WxWidgetMapper(
                self.tFarenheitCtrl, 'Value', wx.EVT_TEXT,
                translator=QuickTranslator(lambda k: str(k*1.8-460.),
                                           lambda f: (float(f)+460.)/1.8)
                )
            )
        

        self.ResetButton = wx.Button(self.frame, label='Reset')
        sizer.Add(self.ResetButton, 0, wx.EXPAND)
        self.controller.BindAction('Reset',
                                   WxActionBinding(self.ResetButton, wx.EVT_BUTTON))

        
        self.DebugButton = wx.Button(self.frame, label='Debug')
        sizer.Add(self.DebugButton, 0, wx.EXPAND)
        self.controller.BindAction('Debug',
                                   WxActionBinding(self.DebugButton, wx.EVT_BUTTON))
        
        
        self.frame.SetSizer(sizer)

        self.frame.Show(True)
        return

class MyApp(wx.App):
    def OnInit(self):
        self.model = Model(273.)
        self.controller = MyController(self.model)
        self.view = View(self.controller)


        self.SetTopWindow(self.view.frame)
        return True
    pass

def test():
    app = MyApp(0)
    app.MainLoop()

if __name__ == "__main__":
    test()
    
