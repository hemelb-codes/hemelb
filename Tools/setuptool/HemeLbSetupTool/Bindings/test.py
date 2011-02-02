import wx

from Observer import Observable

from Controllers import ObjectController
from Mappers import WxWidgetMapper, SimpleObservingMapper
from Translators import QuickTranslator
from Bindings import WxActionBinding

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

    pass


class View(object):
    def __init__(self, controller):
        self.controller = controller

        self.frame = wx.Frame(None, -1, "Hello from wxPython")

        self.tCelciusCtrl = wx.TextCtrl(self.frame)
        cLabel = wx.StaticText(self.frame, label='Celcius')
        cSizer = wx.BoxSizer(wx.HORIZONTAL)
        cSizer.Add(cLabel, 0, wx.EXPAND)
        cSizer.Add(self.tCelciusCtrl, 0, wx.EXPAND)

        self.tFarenheitCtrl = wx.TextCtrl(self.frame)
        fLabel = wx.StaticText(self.frame, label='Farenheit')
        fSizer = wx.BoxSizer(wx.HORIZONTAL)
        fSizer.Add(fLabel, 0, wx.EXPAND)
        fSizer.Add(self.tFarenheitCtrl, 0, wx.EXPAND)

        self.ResetButton = wx.Button(self.frame, label='Reset')

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(cSizer, 0, wx.EXPAND)
        sizer.Add(fSizer, 0, wx.EXPAND)
        sizer.Add(self.ResetButton, 0, wx.EXPAND)

        self.frame.SetSizer(sizer)

        self.controller.BindValue(
            'tKelvin',
            WxWidgetMapper(self.tCelciusCtrl, 'Value', wx.EVT_TEXT),
            translator=QuickTranslator(lambda k: str(k-273.),
                                       lambda c: float(c)+273.)
            )

        self.controller.BindValue(
            'tKelvin',
            WxWidgetMapper(self.tFarenheitCtrl, 'Value', wx.EVT_TEXT),
            translator=QuickTranslator(lambda k: str(k*1.8-460.),
                                       lambda f: (float(f)+460.)/1.8)
            )

        self.controller.BindAction('Reset',
                                   WxActionBinding(self.ResetButton, wx.EVT_BUTTON))

        self.frame.Show(True)
        return

class MyApp(wx.App):
    def OnInit(self):
        self.model = Model(273.)
        self.controller = ObjectController(self.model)
        self.view = View(self.controller)


        self.SetTopWindow(self.view.frame)
        return True
    pass

def test():
    app = MyApp(0)
    app.MainLoop()
