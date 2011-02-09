from ..Bindings.Translators import QuickTranslator
from ..Bindings.ListController import ObjectController, ListController, HasListKeys
from ..Bindings.Mappers import SimpleObservingMapper
from ..Bindings.VtkObject import HasVtkObjectKeys

from ..Model.PlacedIolet import PlacedIolet

from .IoletController import IoletController

import pdb


class PlacedIoletController(HasVtkObjectKeys, ObjectController):
    def __init__(self, delegate):
        ObjectController.__init__(self, delegate)
        self.DefineVtkObjectKey("widget")
        self.DefineVtkObjectKey("representation")
        return
    
    pass

class PlacedIoletListController(ListController):
    """Controller for a list of PlacedIolet objects.
    """
    def __init__(self, delegate):
        ListController.__init__(self, delegate, SelectionControllerClass=PlacedIoletController)
        self.translator = QuickTranslator(self.IoletToPlacedIolet, lambda x: x)
        return

    def IoletToPlacedIolet(self, iolet):
        pi = PlacedIolet(iolet)
        controller = IoletController(iolet)
        # controller.BindValue('Radius', SimpleObservingMapper(pi, 'Radius'))
        return pi
    
    pass

class HasPlacedIoletListKeys(HasListKeys):
    """Mixin for ObjectController subclasses with PlacedIoletList
    keys.
    """
    BindFunctionDispatchTable = ((PlacedIoletListController, 'BindList'),)

    def BindList(self, top, modelKey, widgetMapper):
        HasListKeys.BindList(self, top, modelKey, widgetMapper)
        # widgetMapper.controller = top
        # widgetMapper.key = modelKey
        return
    
    def DefinePlacedIoletListKey(self, name):
        """Typically used in the subclass __init__ method to easily
        mark a key as being a List and hence needing a ListController
        to manage it.
        """
        setattr(self, name,
                PlacedIoletListController(getattr(self.delegate, name))
                )
        return
    
    
    pass

# class VtkPropertyMapper(SimpleObservingMapper):
#     """Custom bindings mapper for the SeedPosition attribute to the
#     'Center' of the VTK pipeline's representation of the seed
#     location.
#     """
#     def __init__(self, model, key, translator=UnitTranslator()):
#         SimpleObservingMapper.__init__(self, model, key, translator)
#         self._vtkObsId = None
#         return
    
#     def Observe(self):
#         SimpleObservingMapper.Observer(self)
#         if self.obsId is not None:
#             self.PlacedSeed.representation.RemoveObserver(self.obsId)
#             pass
        
#         self.obsId = self.PlacedSeed.representation.AddObserver('ModifiedEvent', self.HandleUpdate)
#         return
#     def Unobserve(self):
#         SimpleObservingMapper.Unobserve(self)
        
#         if self.obsId is None:
#             return
        
#         self.PlacedSeed.representation.RemoveObserver(self.obsId)
#         self.obsId = None
#         return
    
#     pass
