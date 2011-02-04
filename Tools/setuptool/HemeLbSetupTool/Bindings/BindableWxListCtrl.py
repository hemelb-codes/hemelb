import wx

from .WxMappers import WxListCtrlItemMapper

class BindableWxListCtrl(wx.ListCtrl):
    """A subclass of wx.ListCtrl that an more easily be bound to a
    ListController. This basic one requires that the controller passed
    to the constructor have keys matching every column name inserted,
    and further that the associated values are strings. Subclasses can
    be used to specialise this behaviour.
    """
    
    def __init__(self, controller, *args, **kwargs):
        wx.ListCtrl.__init__(self, *args, **kwargs)
        self.contentController = controller
        self.itemControllers = []
        self.colData = []
        return
    
    def InsertColumn(self, col, heading, format=wx.LIST_FORMAT_LEFT, width=-1):
        wx.ListCtrl.InsertColumn(self, col, heading, format, width)
        self.colData.insert(col, heading)
        return
    
    def InsertItemAtIndex(self, row, item):
        itemCon = self.contentController[row]
        self.itemControllers.insert(row, itemCon)
        
        self.InsertStringItem(row, '')
        for col, colName in enumerate(self.colData):
            itemCon.BindValue(colName,
                              WxListCtrlItemMapper(self, row, col))
            continue
        
        return
    
    def DeleteItem(self, row):
        del self.itemController[row]
        return wx.ListCtrl.DeleteItem(self, row)
    pass
