from ..Bindings.Translators import QuickTranslator, UnitTranslator
from ..Bindings.ListController import ObjectController, ListController, HasListKeys
from ..Bindings.Mappers import SimpleObservingMapper
from ..Bindings.VtkObject import HasVtkObjectKeys

from ..Model.PlacedIolet import PlacedIolet, PlacedInlet, PlacedOutlet
from ..Model.Iolets import Inlet, Outlet

from .IoletController import IoletController

import pdb


class PlacedIoletController(HasVtkObjectKeys, ObjectController):
    def __init__(self, delegate):
        ObjectController.__init__(self, delegate)
        self.DefineVtkObjectKey("widget")
        self.DefineVtkObjectKey("representation")
        return
    
    pass

# class ComponentTranslator(Translator):
#     def __init__(self, model, component, inner=None):
#         Translator.__init__(self, inner)
#         self.component = component
#         self.model = model
#         self.index = {'x': 0,
#                       'y': 1,
#                       'z': 2}[component.lower()]
#         return
#     def TranslateStage(self, val):
#         orig = self.model.GetValueForKey(self.key)
#         orig[self.index] = val
#         self.model.SetValueForKey(self.key, orig)
#         return
    
class VectorMapper(SimpleObservingMapper):
    """Mapper for Vectors to components of a list.
    """
    
    # def __init__(self, model, key, translator=UnitTranslator()):
    #     SimpleObservingMapper.__init__(self, translator=translator)
        
    #     self.model = model
    #     self.key = key
        
    #     return
    
    def CreateSubMapper(self, component):
        return VectorComponentMapper(self.model,
                                     self.key+'.'+component,
                                     translator=self.translator)
    
    pass

class VectorComponentMapper(SimpleObservingMapper):
    def __init__(self, model, key, translator=UnitTranslator()):
        self.model = model
        self.fullkey = key
        baseKey, self.component = key.rsplit('.', 1)
        self.index = {'x': 0,
                      'y': 1,
                      'z': 2}[self.component.lower()]
        
        SimpleObservingMapper.__init__(self, model, baseKey, translator)
        return
    
    def _Get(self):
        return self.model.GetValueForKey(self.key)[self.index]
    
    def _Set(self, val):
        # Can't set a single component. So get the coords (it's a
        # tuple so must convert to a list to allow item assignment)
        # update and set them all.
        new = list(self.model.GetValueForKey(self.key))
        new[self.index] = val
        self.model.SetValueForKey(self.key, new)
        return

    # We add the index to the change options here so our sibling
    # VCMappers can ignore us.
    def WillChangeValueForKey(self, key):
        return SimpleObservingMapper.WillChangeValueForKey(self, key, index=self.index)
    def DidChangeValueForKey(self, key):
        return SimpleObservingMapper.DidChangeValueForKey(self, key, index=self.index)
    # Here we ignore the update if we have an index in the
    # ChangeOptions object and it's not ours.
    def HandleUpdate(self, change=None, alsoIgnored=None):
        try:
            index = change.index
            if index != self.index:
                return
        except:
            pass
        return SimpleObservingMapper.HandleUpdate(self, change, alsoIgnored)
    pass

class PlacedIoletListController(ListController):
    """Controller for a list of PlacedIolet objects.
    """
    def __init__(self, delegate):
        ListController.__init__(self, delegate, SelectionControllerClass=PlacedIoletController)
        self.translator = QuickTranslator(self.IoletToPlacedIolet, lambda x: x)
        return

    def IoletToPlacedIolet(self, iolet):
        if isinstance(iolet, Inlet):
            pi = PlacedInlet()
        elif isinstance(iolet, Outlet):
            pi = PlacedOutlet()
        else:
            raise ValueError("Must be an Inlet or Outlet")
        
        controller = IoletController(iolet)
        controller.BindValue('Radius', SimpleObservingMapper(pi, 'Radius'))
        controller.BindValue('Centre', VectorMapper(pi, 'Centre'))
        controller.BindValue('Normal', VectorMapper(pi, 'Normal'))
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
