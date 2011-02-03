import wx
import pdb
from ..Util.Observer import ObservableList

from .ObjectController import ObjectController
from .EmptySelection import EmptySelection, isNone
from .Translators import Translator
from .Bindings import ValueBinding
from .Mappers import Mapper, ReadOnlyMapper

class ListController(ObjectController):
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
    
    pass

class ListMapper(Mapper):
    def __init__(self, model):
        Mapper.__init__(self)
        self.model = model
        return
    
    pass

class HasListKeys(object):
    """Mixin for ObjectController subclasses with ObservableList keys.
    """
    BindMethodDispatchTable = ((ListController, 'BindList'),)
    
    def BindList(self, modelKey, widgetMapper):
        """We need to bind the selection and deal with add/remove/update.
        """
        self.BindComplexValue(modelKey, ListContentsMapper, (),
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


class ListContentsMapper(ReadOnlyMapper):
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


