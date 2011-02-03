import wx

# from HemeLbSetupTool.Bindings.Mappers import ListCtrlSelectionBinding

class IoletListCtrl(wx.ListCtrl):
    
    def __init__(self, *args, **kwargs):
        wx.ListCtrl.__init__(self, *args, **kwargs)
        
        self.InsertColumn(0, "Name")
        # controller.BindValue("SelectedIndex", ListCtrlSelectionBinding(self))
        
        return

    
    def InsertItemAtIndex(self, index, item):
        self.InsertStringItem(index, item.Name)
        return
    
