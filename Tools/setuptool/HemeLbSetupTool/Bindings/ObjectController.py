# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from ..Util.Observer import Observable

from .Mappers import SimpleObservingMapper
from .Bindings import ValueBinding

class ObjectController(Observable):
    """Acts as an intermediary between the model object and view objects.
    """
    BindFunctionDispatchTable = ((object, 'BindSimpleValue'), )
    
    def ChooseBindFunction(self, key):
        """Here we walk up the class hierarchy looking in the
        BindFunctionDispatchTable for each class to see if the value
        corresponding to the key is an instance of any entry. We
        return the corresponding method for the first match.

        ObjectController.BindFunctionDispatchTable contains a single
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
                table = controllerClass.BindFunctionDispatchTable
            except AttributeError:
                # If the class doesn't have a BindFunctionDispatchTable
                continue
            
            for cls, methodName in table:
                # Loop through the entries
                if isinstance(value, cls):
                    # We found a match so return the method
                    return getattr(self, methodName)
                continue
            continue
        
        # Search terminated without matching
        raise ValueError('No matching Bind method for Mapper of type "%s"' % str(type(self.widgetMapper)))
    
    def __init__(self, delegate):
        """The delegate will typically be the model or another controller.
        """
        
        assert isinstance(delegate, Observable)
        self.delegate = delegate
        self._values = dict()
        self._actions = set()
        return
    
    def _GetLocalValueForKey(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            return self.delegate._GetLocalValueForKey(key)
        return

    def _SetLocalValueForKey(self, key, value):
        if hasattr(self, key):
            setattr(self, key, value)
        else:
            self.delegate._SetLocalValueForKey(key, value)
            pass
        return
    
    def _AddObserverToLocalKey(self, keyPath, callback, options):
        if hasattr(self, keyPath) or '.' in keyPath:
            # If its our attribute or a dotted path, add to self
            Observable._AddObserverToLocalKey(self, keyPath, callback, options)
        else:
            # Delegate it
            self.delegate._AddObserverToLocalKey(keyPath, callback, options)
            pass
        return
    
    def BindValue(self, modelKey, widgetMapper):
        """Bind a mapper to a value. First find the responsible
        controller by walking the key path, then dispatched on the
        type of the mapper.
        """
        con, conKey = self.FindResponsibleControllerAndKey(modelKey)
        
        con.ChooseBindFunction(conKey)(self, modelKey, widgetMapper)
        return
    
    def FindResponsibleControllerAndKey(self, modelKey):
        """Find the responsible controller by walking along the key
        path, return that and the key path from that point to the full
        key.
        """
        parts = modelKey.split('.', 1)
        if len(parts) == 1:
            return self, parts[0]
        
        local = parts[0]
        rest = parts[1]
        subController = getattr(self, local)
        return subController.FindResponsibleControllerAndKey(rest)
    
    def BindComplexValue(self, topController, modelKey,
                         modelMapperFactory, modelFactoryArgs,
                         bindMgrFactory, widgetMapper):
        try:
            # If the topController's already got a BindingManager for
            # this attribute of the model, use that.
            bindingMgr = topController._values[modelKey]
        except KeyError:
            # If not, create it, based on self if we have that attribute,
            # otherwise on the delegate
            # if hasattr(self, modelKey):
            modelMapper = modelMapperFactory(topController, modelKey, *modelFactoryArgs)
            # else:
            #     modelMapper = modelMapperFactory(self.delegate, modelKey, *modelFactoryArgs)
            #     pass
            bindingMgr = topController._values[modelKey] = bindMgrFactory(modelMapper)
            pass
        # Bind our widget to the manager
        bindingMgr.BindWidget(widgetMapper)
        return
    
    def BindSimpleValue(self, topController, modelKey, widgetMapper):
        """Do a simple value binding.
        """
        self.BindComplexValue(topController, modelKey,
                              SimpleObservingMapper, (), ValueBinding,
                              widgetMapper)
        
        return
    
    def BindAction(self, modelKey, action):
        parts = modelKey.split('.', 1)
        if len(parts) == 1:
            self._actions.add(action)
            action.Bind(self._getCallbackWrapper(modelKey))
        else:
            local = parts[0]
            rest = parts[1]
            subController = getattr(self, local)
            subController.BindAction(rest, action)
            pass

        return
    
    def _getCallbackWrapper(self, key):
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
