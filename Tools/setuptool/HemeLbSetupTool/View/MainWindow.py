import wx
from ToolPanel import ToolPanel
from VtkViewPanel import VtkViewPanel
from Delegator import Delegator, ParentDelegator

class DSplitter(wx.SplitterWindow, ParentDelegator):
    def __init__(self, *args, **kwargs):
        wx.SplitterWindow.__init__(self, *args, **kwargs)
        ParentDelegator.__init__(self)
        return
    pass

class MainWindow(Delegator, wx.Frame):
    def __init__(self, delegate):
        Delegator.__init__(self, delegate)
        wx.Frame.__init__(self, None, title="HemeLB Setup Tool")
                
        self._InitMenu()
        #self._InitStatusBar() # A Statusbar in the bottom of the window
        
        # Going to have a vertically split window; tools on the left,
        # render window on the right.
        self.splitter = DSplitter(self)
        self.exports.toolPanel = ToolPanel(self.splitter)
        self.exports.vtkPanel = VtkViewPanel(self.splitter)
        self.splitter.SplitVertically(self.toolPanel, self.vtkPanel)
        
        self.Show(True)
        
        return
    
    def _InitMenu(self):
        """Do stuff for the menu.
        """
        filemenu = wx.Menu()
        
        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        filemenu.Append(wx.ID_ABOUT, "&About"," Information about this program")
        filemenu.AppendSeparator()
        filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")
        
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        return
     

    pass
