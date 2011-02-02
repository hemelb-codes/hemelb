import wx
from wx.lib.masked import NumCtrl, EVT_NUM

from HemeLbSetupTool.View.Layout import H, V, StretchSpacer
from HemeLbSetupTool.View.VectorCtrl import VectorCtrl, VectorMapper

from HemeLbSetupTool.Bindings.Mappers import WxWidgetMapper, WxWidgetEnabledMapper
from HemeLbSetupTool.Bindings.Translators import NoneToValueTranslator, FloatTranslator
from HemeLbSetupTool.Bindings.Bindings import WxActionBinding

import pdb
class ToolPanel(wx.Panel):
    """Tools Panel for the LHS of the window.
    """
    
    def __init__(self, controller, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        self.controlPanel = ControlPanel(controller, self)
        self.inputPanel = InputPanel(controller, self)
        self.ioletsPanel = IoletsPanel(controller, self)
        controller.BindValue('HaveValidStlFile',
                             WxWidgetEnabledMapper(self.ioletsPanel))
        self.voxelPanel = VoxelPanel(controller, self)
        controller.BindValue('HaveValidStlFile',
                             WxWidgetEnabledMapper(self.voxelPanel))
        self.seedPanel = SeedPanel(controller, self)
        controller.BindValue('HaveValidStlFile',
                             WxWidgetEnabledMapper(self.seedPanel))
        self.outputPanel = OutputPanel(controller, self)
        controller.BindValue('HaveValidStlFile',
                             WxWidgetEnabledMapper(self.outputPanel))
        
        layout = V(
            (self.controlPanel, 1, wx.EXPAND),
            (self.inputPanel, 1, wx.EXPAND),
            (self.ioletsPanel, 3, wx.EXPAND),
            (self.voxelPanel, 1, wx.EXPAND),
            (self.seedPanel, 1, wx.EXPAND),
            (self.outputPanel, 1, wx.EXPAND),
            )
        sizer = layout.create()
        self.SetSizer(sizer)
        return
    pass

class ControlPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup the control panel.
        """
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        self.openProfileButton = wx.Button(self, label='Open Profile')
        self.saveProfileButton = wx.Button(self, label='Save Profile')
        self.generateButton = wx.Button(self, label='Generate')
        self.progressGauge = wx.Gauge(self)
        
        layout = V(
            (H(
                # (self._openProfileButton, 0, wx.ALIGN_LEFT),
                # (self._saveProfileButton, 0, wx.ALIGN_CENTRE),
                # (self._generateButton, 0, wx.ALIGN_RIGHT)
                StretchSpacer(),
                self.openProfileButton,
                self.saveProfileButton,
                self.generateButton,
                StretchSpacer(),
                ), 0, wx.EXPAND),
            self.progressGauge)
        self.SetSizer(layout.create())
        
        return
    pass

class InputPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Set up the input file part of the tool panel.
        """
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        iLabel = wx.StaticText(self, label='Input STL File')
        self.inputStlField = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.controller.BindValue('StlFile',
                                  WxWidgetMapper(
                                      self.inputStlField, 'Value', wx.EVT_TEXT,
                                      translator=NoneToValueTranslator('')
                                      )
                                  )
        self.inputStlButton = wx.Button(self, label='Choose')
        self.controller.BindAction('ChooseStl', WxActionBinding(self.inputStlButton, wx.EVT_BUTTON))
        
        layout = V(
            iLabel,
            (H(
                (V((self.inputStlField, 0, wx.EXPAND)),1, wx.EXPAND),
                self.inputStlButton
                ), 0, wx.EXPAND
             )
            )
        self.SetSizer(layout.create())
        return
    pass

class IoletsDetailPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        centreLabel = wx.StaticText(self, label='Centre (mm)')
        self.centreVector = VectorCtrl(self)
        controller.BindValue(
            'Selection.Centre',
            VectorMapper(self.centreVector, 'Value', wx.EVT_TEXT)
            )
        
        normalLabel = wx.StaticText(self, label='Normal')
        self.normalVector = VectorCtrl(self)
        controller.BindValue(
            'Selection.Normal',
            VectorMapper(self.normalVector, 'Value', wx.EVT_TEXT)
            )
        
        pressureLabel = wx.StaticText(self, label='Pressure (mmHg)')
        self.pressureVector = VectorCtrl(self)
        controller.BindValue(
            'Selection.Pressure',
            VectorMapper(self.pressureVector, 'Value', wx.EVT_TEXT)
            )
        

        self.pressureExpressionLabel = wx.StaticText(self, label='p =')
        controller.BindValue('Selection.PressureEquation',
                                  WxWidgetMapper(self.pressureExpressionLabel,
                                                 'Label', wx.EVT_TEXT,
                                                 translator=NoneToValueTranslator('')))
        layout = V(
            V(centreLabel,
              self.centreVector),
            
            V(normalLabel,
              self.normalVector),
            
            V(pressureLabel,
              self.pressureVector,
              self.pressureExpressionLabel
              )
            )
        self.SetSizer(layout.create())
        
        controller.BindValue('SelectedIndex', WxWidgetEnabledMapper(self))
        
        return
    
    pass

class IoletsPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Set up the inlets and outlets panel.
        """
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        ioletsLabel = wx.StaticText(self, label='Inlets and Outlets')
        
        self.ioletsListCtrl = wx.ListCtrl(self, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        self.ioletsListCtrl.InsertColumn(0, 'Name')
        
        # self.controller.BindValue('Iolets.SelectedIndex',
        #                           ListCtrlSelectionBinding(self.ioletsListCtrl))
        self.addInletButton = wx.Button(self, label='Add Inlet')
        self.addOutletButton = wx.Button(self, label='Add Outlet')
        self.removeIoletButton = wx.Button(self, label='Remove')

        self.detail = IoletsDetailPanel(controller.Iolets, self)
        
        layout = V(ioletsLabel,
                   H( V( (self.ioletsListCtrl, 1, wx.EXPAND) ),
                      self.detail
                     ),
                   H(self.addInletButton,
                     self.addOutletButton,
                     self.removeIoletButton
                     )
                   )
        self.SetSizer(layout.create())
        return
    pass

class VoxelPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup up the voxel size panel.
        """
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        voxelLabel = wx.StaticText(self, label='Voxel size (mm)')
        self.voxelSizeField = NumCtrl(self, style=wx.TE_PROCESS_ENTER, integerWidth=2, fractionWidth=4, allowNegative=False)
        self.controller.BindValue('VoxelSize',
                                  WxWidgetMapper(self.voxelSizeField, 'Value', EVT_NUM,
                                                 translator=NoneToValueTranslator(0.)))
        
        self.voxelResetButton = wx.Button(self, label='Reset')
        self.controller.BindAction('ResetVoxelSize', WxActionBinding(self.voxelResetButton, wx.EVT_BUTTON))
        
        layout = V(
            voxelLabel,
            (H((V((self.voxelSizeField, 0, wx.EXPAND)), 1, wx.EXPAND),
               self.voxelResetButton),0, wx.EXPAND)
            )
        self.SetSizer(layout.create())
        return
    pass

class SeedPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup the seed point panel.
        """
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        seedLabel = wx.StaticText(self, label='Seed Position')
        self.seedVector = VectorCtrl(self)

        self.controller.BindValue('SeedPoint',
                                  VectorMapper(self.seedVector, 'Value', wx.EVT_TEXT))
                                  
        self.seedPlaceButton = wx.Button(self, label='Place')
        
        layout = V(seedLabel,
                   H(self.seedVector,
                   # H(self.seedVector.x,
                   #   self.seedVector.y,
                   #   self.seedVector.z,
                     self.seedPlaceButton)
                   )
        self.SetSizer(layout.create())
        return
    pass

class OutputPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Setup the output file panel.
        """
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        configLabel = wx.StaticText(self, label='Output config')
        self.configField = wx.TextCtrl(self)
        self.controller.BindValue('OutputConfigFile',
                                  WxWidgetMapper(
                                      self.configField, 'Value', wx.EVT_TEXT,
                                      translator=NoneToValueTranslator('')
                                      )
                                  )
        
        self.configChooseButton = wx.Button(self, label='Choose')
        self.controller.BindAction('ChooseOutputConfigFile',
                                   WxActionBinding(self.configChooseButton, wx.EVT_BUTTON))

        xmlLabel = wx.StaticText(self, label='Output xml')
        self.xmlField = wx.TextCtrl(self)
        self.controller.BindValue('OutputXmlFile',
                                  WxWidgetMapper(
                                      self.xmlField, 'Value', wx.EVT_TEXT,
                                      translator=NoneToValueTranslator('')
                                      )
                                  )
        
        self.xmlChooseButton = wx.Button(self, label='Choose')
        self.controller.BindAction('ChooseOutputXmlFile',
                                   WxActionBinding(self.xmlChooseButton, wx.EVT_BUTTON))
        
        layout = V(
            xmlLabel,
            ( H( (V((self.xmlField, 0, wx.EXPAND)), 1, wx.EXPAND),
                 self.xmlChooseButton
                 ), 0, wx.EXPAND),
            configLabel,
            ( H( (V((self.configField, 0, wx.EXPAND)), 1, wx.EXPAND),
                 self.configChooseButton
                 ), 0, wx.EXPAND),
            )

        self.SetSizer(layout.create())
        return 
    
    
    pass
