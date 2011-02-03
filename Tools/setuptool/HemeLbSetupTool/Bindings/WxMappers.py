import types

from .Mappers import Mapper, WriteOnlyMapper
from .Translators import Translator, UnitTranslator, QuickTranslator

class WxWidgetMapper(Mapper):
    def __init__(self, widget, key, event, translator=UnitTranslator()):
        Mapper.__init__(self, translator=translator)
        
        self.widget = widget
        self.key = key
        self.event = event
        
        try:
            # First try the Get/Set pair
            self._Get = getattr(widget, 'Get'+self.key)
            self._Set = getattr(widget, 'Set'+self.key)
        except AttributeError:
            try:
                # Now try using properties
                prop = getattr(type(widget), key)
                assert isinstance(prop, property)
                self._Get = types.MethodType(prop.fget, widget, type(widget))
                self._Set = types.MethodType(prop.fset, widget, type(widget))
            except (AttributeError, AssertionError):
                raise AttributeError("'%s' object has no property '%s' or Get%s/Set%s pair" \
                                     % (self.__class__.__name__, key,key,key))
            pass
        return
    
    def Observe(self):
        self.widget.Bind(self.event, self.HandleUpdate)
        return
    
    def Unobserve(self):
        self.widget.Unbind(self.event)
        return

    pass

class NonObservingWxWidgetMapper(WriteOnlyMapper, WxWidgetMapper):
    def __init__(self, widget, key, translator=UnitTranslator()):
        WxWidgetMapper.__init__(self, widget, key, None, translator=translator)
        return
    
    def Observe(self):
        return
    
    def Unobserve(self):
        return
    pass

class WxWidgetEnabledMapper(NonObservingWxWidgetMapper):
    """Simple binding for enabledness. This is a one way binding, so
    don't want to watch for events.
    """
    def setterHelper(value):
        if value:
            return True
        return False
    
    def __init__(self, widget, translator=QuickTranslator(setterHelper, lambda x: None)):
        NonObservingWxWidgetMapper.__init__(self, widget, 'Enabled', translator=translator)
        return

    setterHelper = staticmethod(setterHelper)
    pass


class WxListCtrlMapper(WriteOnlyMapper):
    """A mapper for the contents of a wx.ListCtrl

    It requires the model end of the binding to be a
    HemeLbSetupTool.Bindings.ListController.ListContentsMapper (or
    subclass thereof) object that will produce appropriate data for
    it. The widget you bind must be a custom subclass of wx.ListCtrl
    that defines the InsertItemAtIndex(int index, item) method to
    actually add the text etc.

    The translator supplied operates on individual items of the
    managed list.
    """
    def __init__(self, widget, translator=UnitTranslator()):
        Mapper.__init__(self, translator=translator)
        
        self.widget = widget
        return

    def Set(self, val):
        """This Set method skips translation here; that is dealt with
        by the Insert handler.
        """
        self.Unobserve()
        self._Set(val)
        self.Observe()

    def _Set(self, val):
        model, change = val
        if change is None:
            self.Clean(model)
            return
        
        if change.key == '@DELETION':
            self.Delete(model, change)
        elif change.key == '@INSERTION':
            self.Insert(model, change)
        elif change.key == '@REPLACEMENT':
            self.Replace(model, change)
        return

    def Clean(self, model):
        self.widget.DeleteAllItems()
        for i, item in enumerate(model):
            self.widget.InsertItemAtIndex(i, self.translator.Translate(item))
            continue
        return
    
    def Insert(self, model, change):
        index = change.index
        self.widget.InsertItemAtIndex(index, self.translator.Translate(model[index]))
        return
    
    def Delete(self, model, change):
        index = change.index
        self.widget.DeleteItem(index)
        return
    
    def Replace (self, model, change):
        self.Delete(model, change)
        self.Insert(model, change)
        return
    
    pass

class WxListCtrlSelectionMapper(Mapper, Translator):
    def __init__(self, widget, inner=None):
        Mapper.__init__(self, translator=self)
        Translator.__init__(self, inner)
        
        self.widget = widget
        return

    def TranslateStage(self, val):
        if isNone(val):
            return -1
        if val < -1 or val >= self.widget.GetItemCount():
            raise IndexError('Index "%d" out of range for %s' % (val, str(widget)))
        
        return val
    
    def UntranslateStage(self, val):
        if val == -1:
            return None
        return val
    
    def Observe(self):
        self.widget.Bind(wx.EVT_LIST_ITEM_SELECTED, self.HandleUpdate)
        self.widget.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.HandleUpdate)
        return
    
    def Unobserve(self):
        self.widget.Unbind(wx.EVT_LIST_ITEM_SELECTED)
        self.widget.Unbind(wx.EVT_LIST_ITEM_DESELECTED)
        return
    
    def _Get(self):
        return self.widget.GetFirstSelected()
    
    def _Set(self, ind):
        self.Unobserve()
        try:
            prevSelected = self._Get()
            if ind != prevSelected:
                self.widget.SetItemState(newInd, wx.LIST_STATE_SELECTED, wx.LIST_STATE_SELECTED)
                self.widget.SetItemState(prevSelected, 0, wx.LIST_STATE_SELECTED)
                pass
            
            
        finally:
            self.Observe()
            pass

        return
    pass
