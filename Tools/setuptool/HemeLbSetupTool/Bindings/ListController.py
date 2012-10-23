# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import collections

from ..Util.Observer import ObservableList

from .ObjectController import ObjectController
from .EmptySelection import EmptySelection
from .Translators import UnitTranslator
from .Bindings import ValueBinding
from .Mappers import ReadOnlyMapper, WriteOnlyMapper

import pdb

class ListController(ObjectController, collections.MutableSequence):
    def __init__(self, delegate, SelectionControllerClass=ObjectController):
        assert isinstance(delegate, ObservableList)
        ObjectController.__init__(self, delegate)
        self.SelectionControllerClass = SelectionControllerClass
        self.SelectedIndex = None
        self.AddDependency('Selection', 'SelectedIndex')
        return
    
    @property
    def Selection(self):
        if self.SelectedIndex is None:
            return self.SelectionControllerClass(EmptySelection)
        return self.SelectionControllerClass(self.delegate[self.SelectedIndex])
    
    def insert(self, index, object):
        self.delegate.insert(index, object)
        self.SelectedIndex = index
        return
    
    def __getitem__(self, index):
        return self.SelectionControllerClass(self.delegate[index])
    def __setitem__(self, index, obj):
        if isinstance(obj, self.SelectionControllerClass):
            self.delegate[index] = obj.delegate
        else:
            self.delegate[index] = obj
            pass
        return
    
    def __delitem__(self, index):
        curLen = len(self.delegate)
        newLen = curLen - 1
        # Now take care to update the selected index such that we take
        # the next item, unless we just popped the item at the end of
        # the list. In that case take the new end. If the list is
        # empty, set it to None.
        newInd = index
        if index >= newLen:
            newInd = newLen - 1
            pass
        if newLen == 0:
            newInd = None
            pass
        self.WillChangeValueForKey('SelectedIndex')
        
        del self.delegate[index]
        # Now, we've pre notified of the SelectedIndex change, so skip the Observable setattr        
        object.__setattr__(self, 'SelectedIndex', newInd)
        # And post the change
        self.DidChangeValueForKey('SelectedIndex')
        return
    
    def __len__(self):
        return self.delegate.__len__()
    
    pass

class HasListKeys(object):
    """Mixin for ObjectController subclasses with ObservableList keys.
    """
    BindFunctionDispatchTable = ((ListController, 'BindList'),)
    
    def BindList(self, top, modelKey, widgetMapper):
        """We need to bind the selection and deal with add/remove/update.
        """
        self.BindComplexValue(top, modelKey, ListContentsSourceMapper, (),
                              ValueBinding, widgetMapper)
        
        return
    
    def DefineListKey(self, name):
        """Typically used in the subclass __init__ method to easily
        mark a key as being a List and hence needing a ListController
        to manage it.
        """
        setattr(self, name,
                ListController(getattr(self.delegate, name))
                )
        return
    
    pass


class ListContentsSourceMapper(ReadOnlyMapper):
    """This is a mapper for list add/remove/replace events.
    """
    def __init__(self, controller, key):
        ReadOnlyMapper.__init__(self)
        listController = getattr(controller, key)
        self.model = listController.delegate
        self.key = key
        self.currentChange = None
        return
    
    def HandleListChange(self, change):
        self.currentChange = change
        try:
            self.HandleUpdate()
        finally:
            self.currentChange = None
            pass
        return

    def _Get(self):
        return self.model, self.currentChange
    
    def Observe(self):
        self.model.AddObserver('@INSERTION', self.HandleListChange)
        self.model.AddObserver('@REMOVAL', self.HandleListChange)
        self.model.AddObserver('@REPLACEMENT', self.HandleListChange)
        return
    
    def Unobserve(self):
        self.model.RemoveObserver('@INSERTION', self.HandleListChange)
        self.model.RemoveObserver('@REMOVAL', self.HandleListChange)
        self.model.RemoveObserver('@REPLACEMENT', self.HandleListChange)
        return
    pass


class ListContentsDestMapper(WriteOnlyMapper):
    """A mapper for the contents of list that you want to shadow that
    of an ObservableList.

    It requires the model end of the binding to be a
    HemeLbSetupTool.Bindings.ListController.ListContentsMapper (or
    subclass thereof) object that will produce appropriate data for
    it. The list you bind must implement the MutableSequence ABC.
    
    The translator supplied operates on individual items of the
    managed list.
    """
    def __init__(self, target, translator=UnitTranslator()):
        WriteOnlyMapper.__init__(self, translator=translator)
        
        self.target = target
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
        
        if change.key == '@REMOVAL':
            self.Delete(model, change)
        elif change.key == '@INSERTION':
            self.Insert(model, change)
        elif change.key == '@REPLACEMENT':
            self.Replace(model, change)
        return

    def Clean(self, model):
        while len(self.target) > 0:
            self.target.pop()
            continue
        
        for item in model:
            self.target.append(self.translator.Translate(item))
            continue
        return
    
    def Insert(self, model, change):
        index = change.index
        # This effectively creates a new control- we will probably
        # have to bind it to a model key
        self.target.insert(index, self.translator.Translate(model[index]))
        return
    
    def Delete(self, model, change):
        index = change.index
        # This effectively deletes a control. We should throw away any
        # bindings
        del self.target[index]
        return
    
    def Replace (self, model, change):
        index = change.index
        self.target[index] = model[index]
        return
    
    pass
