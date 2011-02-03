from HemeLbSetupTool.Util.Observer import Observable

from HemeLbSetupTool.Bindings.Mappers import Mapper, SimpleObservingMapper
from HemeLbSetupTool.Bindings.Bindings import ValueBinding

class ObjectController(Observable):
    """Acts as an intermediary between the model object and view objects.
    """
    BindMethodDispatchTable = ((object, 'BindSimpleValue'), )
    
    def ChooseBindMethod(self, key):
        """Here we walk up the class hierarchy looking in the
        BindMethodDispatchTable for each class to see if the value
        corresponding to the key is an instance of any entry. We
        return the corresponding method for the first match.

        ObjectController.BindMethodDispatchTable contains a single
        entry that forwards any value to 'BindSimpleValue'.
        """
        # First, get the object, possibly from the delegate
        try:
            value = getattr(self, key)
        except AttributeError:
            value = getattr(self.delegate, key)
            pass
        
        # Walk up the hierarchy
        for controllerClass in type(self).__mro__:
            try:
                table = controllerClass.BindMethodDispatchTable
            except AttributeError:
                # If the class doesn't have a BindMethodDispatchTable
                continue
            
            for cls, methodName in table:
                # Loop through the entries
                if isinstance(value, cls):
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
            self.ChooseBindMethod(modelKey)(modelKey, widgetMapper)
        else:
            local = parts[0]
            rest = parts[1]
            subController = getattr(self, local)
            subController.BindValue(rest, widgetMapper)
            pass
        return
    
    def BindComplexValue(self, modelKey, modelMapperFactory, modelFactoryArgs,
                         bindMgrFactory, widgetMapper):
        try:
            # If we've already got a BindingManager for this attribute
            # of the model, use that.
            bindingMgr = self.__values[modelKey]
        except KeyError:
            # If not, create it, on self if we have that attribute,
            # otherwise on the delegate
            if hasattr(self, modelKey):                
                modelMapper = modelMapperFactory(self, modelKey, *modelFactoryArgs)
            else:
                modelMapper = modelMapperFactory(self.delegate, modelKey, *modelFactoryArgs)
                pass
            bindingMgr = self.__values[modelKey] = bindMgrFactory(modelMapper)
            pass
        # Bind our widget to the manager
        bindingMgr.BindWidget(widgetMapper)
        return
    
    def BindSimpleValue(self, modelKey, widgetMapper):
        """Do a simple value binding.
        """
        self.BindComplexValue(modelKey, SimpleObservingMapper, (), ValueBinding, widgetMapper)
        return
    
    def BindAction(self, modelKey, action):
        parts = modelKey.split('.', 1)
        if len(parts) == 1:
            self.__actions.add(action)
            action.Bind(self.__getCallbackWrapper(modelKey))
        else:
            local = parts[0]
            rest = parts[1]
            subController = getattr(self, local)
            subController.BindAction(rest, action)
            pass

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
