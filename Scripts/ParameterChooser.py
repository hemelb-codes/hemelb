#! /usr/bin/env python2.7
import wx

from HemeLbSetupTool.Util.Observer import Observable
from HemeLbSetupTool.Bindings.ObjectController import ObjectController
from HemeLbSetupTool.Bindings.Translators import FloatTranslator
from HemeLbSetupTool.Bindings.WxMappers import WxWidgetMapper

class ParameterModel(Observable):
    def __init__(self):
        # Initialise EVERYTHING
        self.mu_phys  = 4e-3
        self.rho_phys = 1000
        self.rho_latt = 1

        self.deltaP_phys = 0
        self.UMax_phys   = 0
        
        self.delta_x = 1
        self.delta_t = 1

        self.AddDependency('mu_latt', 'nu_latt')
        self.AddDependency('mu_latt', 'rho_latt')

        self.AddDependency('nu_phys', 'mu_phys')
        self.AddDependency('nu_phys', 'rho_phys')

        self.AddDependency('nu_latt', 'nu_phys')
        self.AddDependency('nu_latt', 'delta_x')
        self.AddDependency('nu_latt', 'delta_t')

        self.AddDependency('deltaP_latt', 'deltaP_phys')
        self.AddDependency('deltaP_latt', 'delta_x')
        self.AddDependency('deltaP_latt', 'delta_t')
        self.AddDependency('deltaP_latt', 'rho_phys')
        self.AddDependency('deltaP_latt', 'rho_latt')

        self.AddDependency('UMax_latt', 'UMax_phys')
        self.AddDependency('UMax_latt', 'delta_x')
        self.AddDependency('UMax_latt', 'delta_t')

        self.AddDependency('tau', 'nu_latt')
        self.AddDependency('Mach_number', 'UMax_latt')
        self.AddDependency('rho_diff', 'deltaP_latt')

    @property
    def nu_phys(self):
        return self.mu_phys / self.rho_phys

    @property
    def nu_latt(self):
        return self.nu_phys / (self.delta_x**2 / self.delta_t)

    @property
    def mu_latt(self):
        return self.nu_latt * self.rho_latt

    @property
    def deltaP_latt(self):
        return self.deltaP_phys / ((self.delta_x/self.delta_t)**2 * (self.rho_phys/self.rho_latt))

    @property
    def UMax_latt(self):
        return self.UMax_phys / (self.delta_x / self.delta_t)

    @property
    def tau(self):
        return 0.5 + 3 * self.nu_latt

    @property
    def Mach_number(self):
        return 3**0.5 * self.UMax_latt

    @property
    def rho_diff(self):
        return self.deltaP_latt / (1./3.)

class ParameterPlay(object):
    def __init__(self, controller):
        self.controller = controller

        self.title  = "Parameter choice helper for lattice-Boltzmann"

        frame = wx.Frame(None, -1, self.title)
        self.frame = frame

        self.master_sizer = wx.BoxSizer(wx.VERTICAL)

        # Create components for the physical quantities
        self.mu_label = wx.StaticText(frame, label=u'\u03bc')
        self.mu_textbox = wx.TextCtrl(frame, value="4e-6")
        self.mu_latt = wx.StaticText(frame, label="Placeholder")

        self.rho_label = wx.StaticText(frame, label=u'\u03c1')
        self.rho_textbox = wx.TextCtrl(frame, value="4e-6")
        self.rho_latt = wx.StaticText(frame, label="Placeholder")        

        self.nu_label = wx.StaticText(frame, label=u'\u03bd')
        self.nu_phys = wx.StaticText(frame, label="4e-6")
        self.nu_latt = wx.StaticText(frame, label="Placeholder")

        # Add them to a grid
        self.phys_grid = wx.GridBagSizer(hgap=5, vgap=5)

        self.phys_grid.Add(self.nu_label, pos=(0,0))
        self.phys_grid.Add(self.nu_phys, pos=(0,1))
        self.phys_grid.Add(self.nu_latt, pos=(0,2))

        self.phys_grid.Add(self.mu_label, pos=(1,0))
        self.phys_grid.Add(self.mu_textbox, pos=(1,1))
        self.phys_grid.Add(self.mu_latt, pos=(1,2))

        self.phys_grid.Add(self.rho_label, pos=(2,0))
        self.phys_grid.Add(self.rho_textbox, pos=(2,1))
        self.phys_grid.Add(self.rho_latt, pos=(2,2))

        self.master_sizer.Add(self.phys_grid)

        # Create components for the flow, rather than the fluid
        self.deltaP_label = wx.StaticText(frame, label=u'\u0394 P')
        self.deltaP_textbox = wx.TextCtrl(frame, value="0")
        self.deltaP_latt = wx.StaticText(frame, label="Placeholder")

        self.UMax_label = wx.StaticText(frame, label="U max.")
        self.UMax_textbox = wx.TextCtrl(frame, value="0")
        self.UMax_latt = wx.StaticText(frame, label="Placeholder")

        self.flow_grid = wx.GridBagSizer(hgap=5, vgap=5)

        self.flow_grid.Add(self.deltaP_label, pos=(0,0))
        self.flow_grid.Add(self.deltaP_textbox, pos=(0,1))
        self.flow_grid.Add(self.deltaP_latt, pos=(0,2))

        self.flow_grid.Add(self.UMax_label, pos=(1,0))
        self.flow_grid.Add(self.UMax_textbox, pos=(1,1))
        self.flow_grid.Add(self.UMax_latt, pos=(1,2))

        self.master_sizer.Add(self.flow_grid)

        # Create components for the phys / latt translation
        self.deltaX_label = wx.StaticText(frame, label=u'\u0394 x')
        self.deltaX_textbox = wx.TextCtrl(frame, value="0")

        self.deltaT_label = wx.StaticText(frame, label=u'\u0394 t')
        self.deltaT_textbox = wx.TextCtrl(frame, value="0")

        self.trans_grid = wx.GridBagSizer(hgap=5, vgap=5)

        self.trans_grid.Add(self.deltaX_label, pos=(0,0))
        self.trans_grid.Add(self.deltaX_textbox, pos=(0,1))
        self.trans_grid.Add(self.deltaT_label, pos=(1,0))
        self.trans_grid.Add(self.deltaT_textbox, pos=(1,1))

        self.master_sizer.Add(self.trans_grid)

        # Create components for the critical values
        self.tau_label = wx.StaticText(frame, label=u'\u03c4')
        self.tau_text = wx.StaticText(frame, label="Placeholder")

        self.Mach_number_label = wx.StaticText(frame, label='Ma')
        self.Mach_number_text = wx.StaticText(frame, label="Placeholder")

        self.rho_diff_label = wx.StaticText(frame, label=u'\u03c1 diff')
        self.rho_diff_text = wx.StaticText(frame, label="Placeholder")

        self.computed_grid = wx.GridBagSizer(hgap=5, vgap=5)

        self.computed_grid.Add(self.tau_label, pos=(0,0))
        self.computed_grid.Add(self.tau_text, pos=(0,1))

        self.computed_grid.Add(self.Mach_number_label, pos=(1,0))
        self.computed_grid.Add(self.Mach_number_text, pos=(1,1))

        self.computed_grid.Add(self.rho_diff_label, pos=(2,0))
        self.computed_grid.Add(self.rho_diff_text, pos=(2,1)) 

        self.master_sizer.Add(self.computed_grid)

        # Now we go through and bind all the values
        trans = FloatTranslator()

        self.controller.BindValue('mu_phys', WxWidgetMapper(self.mu_textbox, 'Value', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('mu_latt', WxWidgetMapper(self.mu_latt, 'Label', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('nu_phys', WxWidgetMapper(self.nu_phys, 'Label', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('nu_latt', WxWidgetMapper(self.nu_latt, 'Label', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('rho_phys', WxWidgetMapper(self.rho_textbox, 'Value', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('rho_latt', WxWidgetMapper(self.rho_latt, 'Label', wx.EVT_TEXT,translator=trans))

        self.controller.BindValue('deltaP_phys', WxWidgetMapper(self.deltaP_textbox, 'Value', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('deltaP_latt', WxWidgetMapper(self.deltaP_latt, 'Label', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('UMax_phys', WxWidgetMapper(self.UMax_textbox, 'Value', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('UMax_latt', WxWidgetMapper(self.UMax_latt, 'Label', wx.EVT_TEXT,translator=trans))

        self.controller.BindValue('delta_x', WxWidgetMapper(self.deltaX_textbox, 'Value', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('delta_t', WxWidgetMapper(self.deltaT_textbox, 'Value', wx.EVT_TEXT,translator=trans))

        self.controller.BindValue('tau', WxWidgetMapper(self.tau_text, 'Label', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('Mach_number', WxWidgetMapper(self.Mach_number_text, 'Label', wx.EVT_TEXT,translator=trans))
        self.controller.BindValue('rho_diff', WxWidgetMapper(self.rho_diff_text, 'Label', wx.EVT_TEXT,translator=trans))

        frame.SetSizer(self.master_sizer)
        self.frame.Show()
        self.master_sizer.Fit(self.frame)

class ParameterChooser(wx.App):
    def OnInit(self):
        model = ParameterModel()
        controller = ObjectController(model)
        view =  ParameterPlay(controller)

        self.SetTopWindow(view.frame)

        return 1

if __name__ == "__main__":
    programme = ParameterChooser(0)
    programme.MainLoop()
