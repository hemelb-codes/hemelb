import wx
from wx.lib.masked import NumCtrl

class VectorCtrl(object):
    """Simple container of three TextCtrl's for a vector quantity.
    """
    
    def __init__(self, parent):
        self.x = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, fractionWidth=6)
        self.y = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, fractionWidth=6)
        self.z = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, fractionWidth=6)
        return
    
    pass

class LayoutContainer(object):
    """Represents a BoxSizer. It's children can be other
    LayoutContainers, wx.Window subclasses or Spacer subclasses. The
    first two cases can also be given as tuples of (object,
    proportion, flags).
    
    """
    def __init__(self, *args):
        self.children = args
        return

    def create(self):
        """Create the sizer hierarchy."""
        sizer = wx.BoxSizer(self.orient)
        
        for c in self.children:
            if isinstance(c, tuple):
                child = c[0]
                args = c[1:]
            else:
                child = c
                args = tuple()
                pass
            
            if isinstance(child, LayoutContainer):
                subElem = child.create()
            else:
                subElem = child
                pass

            if isinstance(subElem, StretchSpacer):
                sizer.AddStretchSpacer(prop=subElem.prop)
            else:
                sizer.Add(subElem, *args)
                pass
            
            continue
        return sizer
    
    pass

class StretchSpacer(object):
    def __init__(self, prop=1):
        self.prop = prop
        return
    pass

class V(LayoutContainer):
    orient = wx.VERTICAL
    pass

class H(LayoutContainer):
    orient = wx.HORIZONTAL
    pass

