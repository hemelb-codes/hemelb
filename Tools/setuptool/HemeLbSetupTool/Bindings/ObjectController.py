from HemeLbSetupTool.Util.Observer import Observable

from HemeLbSetupTool.Bindings.Mappers import Mapper, SimpleObservingMapper
from HemeLbSetupTool.Bindings.Bindings import ValueBinding

class ObjectController(Observable):
    """Acts as an intermediary between the model object and view objects.
    """
    WidgetMapperClassDispatchTable = ((Mapper, 'BindSimpleValue'), )
    
    def ChooseBindMethod(self, widgetMapper):
        """Here we walk up the class hierarchy looking in the
        WidgetMapperClassDispatchTable for each class to see if the
        given widget mapper is an instance of any entry. We return the
        corresponding method for the first match.

        ObjectController.WidgetMapperClassDispatchTable contains a
        single entry that forwards any Mapper subclass to
        'BindSimpleValue'.
        """
        # Walk up the hierarchy
        for controllerClass in type(self).__mro__:
            try:
                table = controllerClass.WidgetMapperClassDispatchTable
            except AttributeError:
                # If the class doesn't have a WidgetMapperClassDispatchTable
                continue
        
            for mapperClass, methodName in table:
                # Loop through the entries
                if isinstance(widgetMapper, mapperClass):
                    # We found a match so return the method
                    return getattr(self, methodName)
                continue
            continue
        
        # Search terminated without matching
        raise ValueError('No matching Bind method for Mapper of type "%s"' % str(type(widgetMapper)))
    
    def __init__(self, delegate):
        """The delegate will typically be the model or another controller.
        """
        
        assert isinstance(delegate, Observable)
        self.delegate = delegate
        self.__values = dict()
        self.__actions = set()
        return
    
    def BindValue(self, modelKey, widgetMapper):
        """Bind a mapper to a value. First find the responsible
        controller by walking the key path, then dispatched on the
        type of the mapper.
        """
        parts = modelKey.split('.', 1)
        if len(parts) == 1:
            self.ChooseBindMethod(widgetMapper)(modelKey, widgetMapper)
        else:
            local = parts[0]
            rest = parts[1]
            subController = getattr(self, local)
            subController.BindValue(rest, widgetMapper)
            pass
        return
    
    def BindSimpleValue(self, modelKey, widgetMapper):
        """Do a simple value binding.
        """
        try:
            # If we've already got a BindingManager for this attribute
            # of the model, use that.
            bindingMgr = self.__values[modelKey]
        except KeyError:
            # If not, create it, on self if we have that attribute,
            # otherwise on the delegate
            if hasattr(self, modelKey):                
                modelMapper = SimpleObservingMapper(self, modelKey)
            else:
                modelMapper = SimpleObservingMapper(self.delegate, modelKey)
                pass
            bindingMgr = self.__values[modelKey] = ValueBinding(modelMapper)
            pass
        # Bind our widget to the manager
        bindingMgr.BindWidget(widgetMapper)
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
                obj = obj.delegate
            except AttributeError:
                raise AttributeError('No matching key found')
            
            continue

        def cbwrapper(self, ignored=None):
            cb()
            return
        return cbwrapper
    pass
