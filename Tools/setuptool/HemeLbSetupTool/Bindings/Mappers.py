from Observer import Observable

class Mapper(object):
    
    def SetBinding(self, b):
        self.binding = b
        return
    
    def HandleUpdate(self, ignored):
        self.binding.MapperWasUpdated(self)
        return
    
    pass

class SimpleObservingMapper(Mapper):
    def __init__(self, model, key):
        assert isinstance(model, Observable)
        
        self.model = model
        self.key = key
        return
    
    def Observe(self):
        self.model.AddObserver(self.key, self.HandleUpdate)
        return
    
    def Unobserve(self):
        self.model.RemoveObserver(self.key, self.HandleUpdate)
        return
    
    def Get(self):
        return getattr(self.model, self.key)
    
    def Set(self, val):
        self.Unobserve()
        setattr(self.model, self.key, val)
        self.Observe()
        return
    pass

class WxWidgetMapper(Mapper):
    def __init__(self, widget, key, event):
        self.widget = widget
        self.key = key
        self.event = event
        
        self.Get = getattr(widget, 'Get'+self.key)
        self.setter = getattr(widget, 'Set'+self.key)
        return
    
    def Observe(self):
        self.widget.Bind(self.event, self.HandleUpdate)
        return
    
    def Unobserve(self):
        self.widget.Bind(self.event, self.HandleUpdate)
        return

    def Set(self, val):
        self.Unobserve()
        self.setter(val)
        self.Observe()
        return
    
    pass
