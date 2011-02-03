import wx
from wx.lib.masked import NumCtrl, EVT_NUM

from HemeLbSetupTool.View.Layout import H, V, StretchSpacer, RectSpacer
from HemeLbSetupTool.View.VectorCtrl import VectorCtrl, VectorCtrlMapper
from HemeLbSetupTool.View.IoletListCtrl import IoletListCtrl

from HemeLbSetupTool.Bindings.WxMappers import WxWidgetMapper, WxWidgetEnabledMapper, NonObservingWxWidgetMapper, WxListCtrlMapper
from HemeLbSetupTool.Bindings.Translators import NoneToValueTranslator, FloatTranslator, QuickTranslator
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
            (self.controlPanel, 0, wx.EXPAND),
            (self.inputPanel, 0, wx.EXPAND),
            (self.ioletsPanel, 1, wx.EXPAND),
            (self.voxelPanel, 0, wx.EXPAND),
            (self.seedPanel, 0, wx.EXPAND),
            (self.outputPanel, 0, wx.EXPAND),
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
        controller.BindValue('IsReadyToGenerate', WxWidgetEnabledMapper(self.generateButton))
        
        self.progressGauge = wx.Gauge(self)
        
        layout = V(
            (H(
                StretchSpacer(),
                self.openProfileButton,
                self.saveProfileButton,
                self.generateButton,
                StretchSpacer(),
                ), 0, wx.EXPAND),
            (self.progressGauge, 0, wx.EXPAND)
            )
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
        controller.BindValue('StlFile',
                                  WxWidgetMapper(
                                      self.inputStlField, 'Value', wx.EVT_TEXT,
                                      translator=NoneToValueTranslator('')
                                      )
                                  )
        controller.BindValue('HaveValidStlFile',
                             NonObservingWxWidgetMapper(self.inputStlField, 'BackgroundColour',
                                                        translator=controller.validColourer)
                             )
        
        self.inputStlButton = wx.Button(self, label='Choose')
        controller.BindAction('ChooseStl', WxActionBinding(self.inputStlButton, wx.EVT_BUTTON))
        
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
    @staticmethod
    def translatorHelper(val):
        if isNone(val):
            return False
        return True
    
    def __init__(self, controller, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        centreLabel = wx.StaticText(self, label='Centre (mm)')
        self.centreVector = VectorCtrl(self)
        controller.BindValue(
            'Selection.Centre',
            VectorCtrlMapper(self.centreVector, 'Value', wx.EVT_TEXT)
            )
        
        normalLabel = wx.StaticText(self, label='Normal')
        self.normalVector = VectorCtrl(self)
        controller.BindValue(
            'Selection.Normal',
            VectorCtrlMapper(self.normalVector, 'Value', wx.EVT_TEXT)
            )
        
        pressureLabel = wx.StaticText(self, label='Pressure (mmHg)')
        self.pressureVector = VectorCtrl(self)
        controller.BindValue(
            'Selection.Pressure',
            VectorCtrlMapper(self.pressureVector, 'Value', wx.EVT_TEXT)
            )
        

        self.pressureExpressionLabel = wx.StaticText(self, label='p =')
        controller.BindValue('Selection.PressureEquation',
                                  WxWidgetMapper(self.pressureExpressionLabel,
                                                 'Label', wx.EVT_TEXT,
                                                 translator=NoneToValueTranslator('')))
        layout = V(
            (V(centreLabel,
              (self.centreVector, 1, wx.EXPAND)),
             0, wx.EXPAND),
            
            (V(normalLabel,
               (self.normalVector, 1, wx.EXPAND)),
            0, wx.EXPAND),
            
            (V(pressureLabel,
               (self.pressureVector, 1, wx.EXPAND),
               self.pressureExpressionLabel
               ),
             0, wx.EXPAND),
            RectSpacer(0,2)
            )
        self.SetSizer(layout.create())
        
        controller.BindValue(
            'SelectedIndex',
            WxWidgetEnabledMapper(
                self,
                translator=QuickTranslator(self.translatorHelper,
                                           lambda x: None)
                )
            )
        
        return
    
    pass

class IoletsPanel(wx.Panel):
    def __init__(self, controller, *args, **kwargs):
        """Set up the inlets and outlets panel.
        """
        wx.Panel.__init__(self, *args, **kwargs)
        self.controller = controller
        
        ioletsLabel = wx.StaticText(self, label='Inlets and Outlets')
        
        self.ioletsListCtrl = IoletListCtrl(self, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        
        controller.BindValue('Iolets',
                             WxListCtrlMapper(self.ioletsListCtrl))
        # controller.BindValue('Iolets.SelectedIndex',
        #                           ListCtrlSelectionBinding(self.ioletsListCtrl))
        self.addInletButton = wx.Button(self, label='Add Inlet')
        controller.BindAction('Iolets.AddInlet',
                              WxActionBinding(self.addInletButton, wx.EVT_BUTTON))

        self.addOutletButton = wx.Button(self, label='Add Outlet')
        controller.BindAction('Iolets.AddOutlet',
                              WxActionBinding(self.addOutletButton, wx.EVT_BUTTON))
        
        self.removeIoletButton = wx.Button(self, label='Remove')
        controller.BindValue('Iolets.SelectedIndex',
                             WxWidgetEnabledMapper(self.removeIoletButton))
        controller.BindAction('Iolets.RemoveIolet',
                              WxActionBinding(self.removeIoletButton, wx.EVT_BUTTON))
        
        self.detail = IoletsDetailPanel(controller.Iolets, self)
        
        layout = V(ioletsLabel,
                   (H( (V( (self.ioletsListCtrl, 1, wx.EXPAND) ), 1, wx.EXPAND),
                      (V(StretchSpacer(), (self.detail, 0, wx.EXPAND)), 2, wx.EXPAND)
                     ), 1, wx.EXPAND),
                   (H(self.addInletButton,
                     self.addOutletButton,
                     self.removeIoletButton
                     ), 0, wx.EXPAND),
                   RectSpacer(0,1)
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
            (H((V(self.voxelSizeField), 1, wx.EXPAND),
               self.voxelResetButton), 1, wx.EXPAND)
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

        controller.BindValue('SeedPoint',
                             VectorCtrlMapper(self.seedVector, 'Value', wx.EVT_TEXT))
        controller.BindValue('HaveValidSeedPoint',
                             NonObservingWxWidgetMapper(self.seedVector, 'BackgroundColour',
                                                        translator=controller.validColourer)
                             )
        
        self.seedPlaceButton = wx.Button(self, label='Place')
        
        layout = V(seedLabel,
                   (H((self.seedVector, 1, wx.EXPAND),
                     self.seedPlaceButton), 1, wx.EXPAND)
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
        controller.BindValue('OutputConfigFile',
                             WxWidgetMapper(
                                 self.configField, 'Value', wx.EVT_TEXT,
                                 translator=NoneToValueTranslator('')
                                 )
                             )
        controller.BindValue('HaveValidOutputConfigFile',
                             NonObservingWxWidgetMapper(self.configField, 'BackgroundColour',
                                                        translator=controller.validColourer)
                             )

        self.configChooseButton = wx.Button(self, label='Choose')
        controller.BindAction('ChooseOutputConfigFile',
                              WxActionBinding(self.configChooseButton, wx.EVT_BUTTON))

        xmlLabel = wx.StaticText(self, label='Output xml')
        self.xmlField = wx.TextCtrl(self)
        controller.BindValue('OutputXmlFile',
                             WxWidgetMapper(
                                 self.xmlField, 'Value', wx.EVT_TEXT,
                                 translator=NoneToValueTranslator('')
                                 )
                             )
        controller.BindValue('HaveValidOutputXmlFile',
                             NonObservingWxWidgetMapper(self.xmlField, 'BackgroundColour',
                                                        translator=controller.validColourer)
                             )

        self.xmlChooseButton = wx.Button(self, label='Choose')
        controller.BindAction('ChooseOutputXmlFile',
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
