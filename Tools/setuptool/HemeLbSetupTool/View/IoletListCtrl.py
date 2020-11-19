# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from HemeLbSetupTool.Bindings.BindableWxListCtrl import BindableWxListCtrl

class IoletListCtrl(BindableWxListCtrl):
    def __init__(self, *args, **kwargs):
        BindableWxListCtrl.__init__(self, *args, **kwargs)
        self.InsertColumn(0, "Name")
        return
    
    pass
    
