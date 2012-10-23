# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from .EmptySelection import isNone

class ValidationError(Exception):
    pass
class FormattingError(Exception):
    pass

class Translator(object):
    """Base class for all translators.
    """
    def __init__(self, inner=None):
        self.inner = inner
        return
    
    def Translate(self, val):
        if self.inner is None:
            return self.TranslateStage(val)
        else:
            return self.inner.Translate(self.TranslateStage(val))
        return
    def Untranslate(self, val):
        if self.inner is None:
            return self.UntranslateStage(val)
        else:
            return self.UntranslateStage(self.inner.Untranslate(val))
        return
    
    def TranslateStage(self, val):
        raise NotImplementedError
    def UntranslateStage(self, val):
        raise NotImplementedError
    pass

class QuickTranslator(Translator):
    """Simple translator using functions supplied to the constructor
    to go forwards and backwards.
    """
    def __init__(self, forwards, backwards, inner=None):
        Translator.__init__(self, inner)
        self.forwards = forwards
        self.backwards = backwards
        return
    
    def TranslateStage(self, val):
        try:
            return self.forwards(val)
        except:
            raise FormattingError()
        
    def UntranslateStage(self, val):
        try:
            return self.backwards(val)
        except:
            raise ValidationError()
        
    pass

class UnitTranslator(Translator):
    """This translator just passes its values through.
    """
    def TranslateStage(self, val):
        return val
    def UntranslateStage(self, val):
        return val
    pass

class NoneToValueTranslator(UnitTranslator):
    """This translator converts None (or the EmptySelection object) to
    the object specified in constructor.
    """
    def __init__(self, valueForNone, inner=None):
        UnitTranslator.__init__(self, inner)
        self.default = valueForNone
        return
    
    def TranslateStage(self, val):
        if isNone(val):
            return self.default
        return val
    
    pass

class FloatTranslator(Translator):
    """Translates from a float to strings.
    """
    def __init__(self, format='%g', inner=None):
        Translator.__init__(self, inner)
        self.format = format
        return

    def TranslateStage(self, value):
        try:
            return self.format % value
        except TypeError:
            return 'nan'
        return
    
    def UntranslateStage(self, value):
        try:
            return float(value)
        except ValueError:
            raise ValidationError('Cannot convert value "%s" to a float' % value)
        
    pass


        
