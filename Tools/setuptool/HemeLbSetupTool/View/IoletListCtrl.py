from HemeLbSetupTool.Bindings.BindableWxListCtrl import BindableWxListCtrl

class IoletListCtrl(BindableWxListCtrl):
    def __init__(self, *args, **kwargs):
        BindableWxListCtrl.__init__(self, *args, **kwargs)
        self.InsertColumn(0, "Name")
        return
    
    pass
    
