from Observer import Observable

from Translators import NullTranslator
from Mappers import SimpleObservingMapper
from Bindings import ValueBinding

class ObjectController(Observable):
    """Acts as an intermediary between the model object and view objects.
    """
    def __init__(self, delegate):
        """The delegate will typically be the model or another controller.
        """
        
        assert isinstance(delegate, Observable)
        self.__delegate = delegate
        self.__values = set()
        self.__actions = set()
        return
    
    def BindValue(self, modelKey, widgetMapper, translator=NullTranslator()):
        modelMapper = SimpleObservingMapper(self.__delegate, modelKey)
        self.__values.add(ValueBinding(modelMapper, widgetMapper, translator))
        return
    
    def BindAction(self, key, action):
        self.__actions.add(action)
        action.Bind(self.__getCallbackWrapper(key))
        return
    
    def __getCallbackWrapper(self, key):
        obj = self
        cb = None
        while cb is None:
            if hasattr(obj, key):
                cb = getattr(obj, key)
                break
                pass
            
            try:
                obj = obj.__delegate
            except AttributeError:
                raise AttributeError('No matching key found')
            
            continue

        def cbwrapper(self, ignored=None):
            cb()
            return
        return cbwrapper
