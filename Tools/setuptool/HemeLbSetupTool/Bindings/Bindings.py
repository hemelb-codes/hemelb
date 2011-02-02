from Translators import NullTranslator

class ValueBinding(object):
    def __init__(self, modelMapper, widgetMapper, translator=NullTranslator()):
        modelMapper.SetBinding(self)
        widgetMapper.SetBinding(self)
        self.modelMapper = modelMapper
        self.widgetMapper = widgetMapper
        self.translator = translator
        
        modelMapper.Observe()
        widgetMapper.Observe()
        
        self.MapperWasUpdated(self.modelMapper)
        return
    
    def MapperWasUpdated(self, which):
        if which is self.modelMapper:
            try:
                self.widgetMapper.Set(
                    self.translator.ModelToWidget(
                        self.modelMapper.Get()
                        )
                    )
            except FormattingError:
                pass
            
        elif which is self.widgetMapper:
            try:
                self.modelMapper.Set(
                    self.translator.WidgetToModel(
                        self.widgetMapper.Get()
                        )
                    )
            except ValidationError:
                self.widgetMapper.Set(
                    self.translator.ModelToWidget(
                        self.modelMapper.Get()
                        )
                    )
        else:
            raise ValueError('Mapper is not one of mine!')
        return
    
    pass

class ActionBinding(object):
    pass

class WxActionBinding(ActionBinding):
    def __init__(self, widget, event):
        self.widget = widget
        self.event = event
        return

    def Bind(self, callback):
        self.widget.Bind(self.event, callback)
        return
    pass
