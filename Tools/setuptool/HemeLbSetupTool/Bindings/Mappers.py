from HemeLbSetupTool.Util.Observer import Observable
from HemeLbSetupTool.Bindings.EmptySelection import isNone
from HemeLbSetupTool.Bindings.Translators import UnitTranslator

class Mapper(object):
    def __init__(self, translator=UnitTranslator()):
        self.translator = translator
        return
    
    def SetTranslator(self, trans):
        self.translator = trans
        return
    
    def SetBinding(self, b):
        self.binding = b
        return
    
    def HandleUpdate(self, ignored):
        self.binding.MapperWasUpdated(self)
        return

    def Get(self):
        ans = self._Get()
        ans = self.translator.Untranslate(ans)
        return ans
    
    def Set(self, val):
        self.Unobserve()
        val = self.translator.Translate(val)
        self._Set(val)
        self.Observe()
    pass

class SimpleObservingMapper(Mapper):
    def __init__(self, model, key, translator=UnitTranslator()):
        assert isinstance(model, Observable)
        Mapper.__init__(self, translator=translator)
        self.model = model
        self.key = key
        return
    
    def Observe(self):
        self.model.AddObserver(self.key, self.HandleUpdate)
        return
    
    def Unobserve(self):
        self.model.RemoveObserver(self.key, self.HandleUpdate)
        return
    
    def _Get(self):
        return getattr(self.model, self.key)
    
    def _Set(self, val):
        setattr(self.model, self.key, val)
        return
    pass

class WxWidgetMapper(Mapper):
    def __init__(self, widget, key, event, translator=UnitTranslator()):
        Mapper.__init__(self, translator=translator)
        
        self.widget = widget
        self.key = key
        self.event = event

        self._Get = getattr(widget, 'Get'+self.key)
        self._Set = getattr(widget, 'Set'+self.key)
        return
    
    def Observe(self):
        self.widget.Bind(self.event, self.HandleUpdate)
        return
    
    def Unobserve(self):
        self.widget.Unbind(self.event)
        return

    pass

class WxWidgetEnabledMapper(Mapper):
    """Simple binding for enabledness. This is a one way binding, so
    don't want to watch for events.
    """
    def __init__(self, widget):
        Mapper.__init__(self)
        self.widget = widget
        return
    
    def _Get(self):
        return self.widget.Enabled
    def _Set(self, value):
        if value:
            self.widget.Enable()
        else:
            self.widget.Disable()
            pass
        return
    
    def Observe(self):
        return
    
    def Unobserve(self):
        return
    
