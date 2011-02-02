import wx
from wx.lib.masked import NumCtrl, EVT_NUM

from Layout import H, V, StretchSpacer, VectorCtrl
from Delegator import ParentDelegator

class ToolPanel(wx.Panel, ParentDelegator):
    """Tools Panel for the LHS of the window.
    """
    
    def __init__(self, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        ParentDelegator.__init__(self)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        layout = V(
            (self._InitControlPanel(), 1, wx.EXPAND),
            (self._InitInputPanel(), 1, wx.EXPAND),
            (self._InitIoletsPanel(), 3, wx.EXPAND),
            (self._InitVoxelPanel(), 1, wx.EXPAND),
            (self._InitSeedPanel(), 1, wx.EXPAND),
            (self._InitOutputPanel(), 1, wx.EXPAND),
            )
        sizer = layout.create()
        self.SetSizer(sizer)
        return
    
    def _InitInputPanel(self):
        """Set up the input file part of the tool panel.
        """
        iLabel = wx.StaticText(self, label='Input STL File')
        self.exports.inputStlField = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.exports.inputStlButton = wx.Button(self, label='Choose')
        
        layout = V(
            iLabel,
            (H(
                (V((self.exports.inputStlField, 0, wx.EXPAND)),1, wx.EXPAND),
                self.exports.inputStlButton
                ), 0, wx.EXPAND
             )
            )
        return layout

    def _InitIoletsPanel(self):
        """Set up the inlets and outlets panel.
        """
        ioletsLabel = wx.StaticText(self, label='Inlets and Outlets')
        self.ioletsList = wx.ListCtrl(self)
        
        centreLabel = wx.StaticText(self, label='Centre')
        self.centreVector = VectorCtrl(self)
        
        normalLabel = wx.StaticText(self, label='Normal')
        self.normalVector = VectorCtrl(self)
        
        pressureLabel = wx.StaticText(self, label='Pressure')
        self.pressureVector = VectorCtrl(self)

        self.exports.addInletButton = wx.Button(self, label='Add Inlet')
        self.exports.addOutletButton = wx.Button(self, label='Add Outlet')
        self.exports.removeIoletButton = wx.Button(self, label='Remove')
        
        self.pressureExpressionLabel = wx.StaticText(self, label='p =')
        layout = V(ioletsLabel,
                   H(self.ioletsList, V(
                       V(centreLabel,
                         H(self.centreVector.x,
                           self.centreVector.y,
                           self.centreVector.z)
                         ),
                       
                       V(normalLabel,
                         H(self.normalVector.x,
                           self.normalVector.y,
                           self.normalVector.z)
                         ),
                       
                       V(pressureLabel,
                         H(self.pressureVector.x,
                           self.pressureVector.y,
                           self.pressureVector.z),
                         self.pressureExpressionLabel
                         ),
                       
                       H(self.addInletButton,
                         self.addOutletButton,
                         self.removeIoletButton
                         )
                       )
                     )
                   )
        return layout

    def _InitVoxelPanel(self):
        """Setup up the voxel size panel.
        """
        voxelLabel = wx.StaticText(self, label='Voxel size (mm)')
        self.exports.voxelSizeField = NumCtrl(self, style=wx.TE_PROCESS_ENTER, fractionWidth=10, allowNegative=False)
        self.exports.voxelResetButton = wx.Button(self, label='Reset')

        layout = V(
            voxelLabel,
            (H((V((self.voxelSizeField, 0, wx.EXPAND)), 1, wx.EXPAND),
               self.voxelResetButton),0, wx.EXPAND)
            )
            
        return layout

    def _InitSeedPanel(self):
        """Setup the seed point panel.
        """
        seedLabel = wx.StaticText(self, label='Seed Position')
        self.exports.seedVector = VectorCtrl(self)
        self.exports.seedPlaceButton = wx.Button(self, label='Place')
        
        layout = V(seedLabel,
                   H(self.seedVector.x,
                     self.seedVector.y,
                     self.seedVector.z,
                     self.seedPlaceButton)
                   )
        return layout

    def _InitOutputPanel(self):
        """Setup the output file panel.
        """
        configLabel = wx.StaticText(self, label='Output config')
        self.exports.configField = wx.TextCtrl(self)
        self.exports.configChooseButton = wx.Button(self, label='Choose')
        xmlLabel = wx.StaticText(self, label='Output xml')
        self.exports.xmlField = wx.TextCtrl(self)
        self.exports.xmlChooseButton = wx.Button(self, label='Choose')
        
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
        
        return layout
    
    def _InitControlPanel(self):
        """Setup the control panel.
        """
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
        return layout
    
    pass
