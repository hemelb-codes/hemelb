class ValidationError(Exception):
    pass
class FormattingError(Exception):
    pass

class Translator(object):
    def ModelToWidget(self, val):
        raise NotImplementedError
    def WidgetToModel(self, val):
        raise NotImplementedError
    pass

class QuickTranslator(Translator):
    def __init__(self, m2w, w2m):
        self.m2w = m2w
        self.w2m = w2m
        return
    
    def ModelToWidget(self, val):
        try:
            return self.m2w(val)
        except:
            raise FormattingError()
        
    def WidgetToModel(self, val):
        try:
            return self.w2m(val)
        except:
            raise ValidationError()
        
    pass

class NullTranslator(Translator):
    def ModelToWidget(self, val):
        return val
    def WidgetToModel(self, val):
        return val
    pass

