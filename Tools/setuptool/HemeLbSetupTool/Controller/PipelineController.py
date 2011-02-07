import operator
from vtk import vtkSphereSource, vtkPolyDataMapper, vtkActor, vtkInteractorStyleTrackballCamera, vtkInteractorStyleUser

from ..Util.Observer import Observable
from ..Bindings.ObjectController import ObjectController
from ..Bindings.VtkObject import HasVtkObjectKeys
from ..Bindings.Mappers import Mapper
from ..Bindings.Translators import UnitTranslator

from .VectorController import VectorController
import pdb

# class VtkPointController(VectorController):
#     def __init__(self, delegate):
#         VectorController.__init__(self, delegate)
        
#         self.representation = vtkSphereSource()
#         self.mapper = vtkPolyDataMapper()
#         self.mapper.SetInputConnection(self.representation.GetOutputPort())
        
#         self.actor = vtkActor()
#         self.actor.SetMapper(self.mapper)
#         # Make it blue
#         self.actor.GetProperty().SetColor(0,0,1)
        
#         return

#     pass

# class HasVtkPointKeys(object):
#     """Mixin for ObjectController subclasses with VtkPoint (Vectors
#     that will shown as points in the view) keys.
#     """
#     BindFunctionDispatchTable = ((VtkPointController, 'BindVector'),)
#     pass
class SeedCoordMapper(Mapper):
    def __init__(self, i, placer, translator=UnitTranslator()):
        Mapper.__init__(self, translator)
        
        self.placer = placer
        self.i = i
        self.obsId = None
        return
    
    def _Get(self):
        return self.placer.representation.GetCenter()[self.i]
    
    def _Set(self, val):
        if val != val:
            # It's a nan
            self.placer.On = False
            pass
        
        old = list(self.placer.representation.GetCenter())
        old[self.i] = val
        self.placer.representation.SetCenter(old)

        allFinite = reduce(operator.and_, map(lambda x: x==x, old))
        if allFinite:
            self.placer.On = True
            pass
        return

    def Observe(self):
        if self.obsId is not None:
            self.Unobserve()
            pass
        
        self.obsId = self.placer.representation.AddObserver('ModifiedEvent', self.HandleUpdate)
        return
    def Unobserve(self):
        if self.obsId is None:
            return
        
        self.placer.representation.RemoveObserver(self.obsId)
        self.obsId = None
        return
    
    pass


class SeedPlacer(Observable):
    def __init__(self, controller, colour=(0,0,1)):
        self.controller = controller
        
        self.representation = vtkSphereSource()
        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.representation.GetOutputPort())
        
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        # Make it blue
        self.actor.GetProperty().SetColor(colour)

        self.On = False
        self.AddObserver('On', self.OnSet)
        
        return
    
    def SetCentre(self, centre):
        self.representation.SetCenter(centre)
        return
    
    def OnSet(self, change):
        if self.On:
            if not self.controller.IsActorAdded(self.actor):
                self.GetValueForKey('controller.Renderer').AddActor(self.actor)
                pass
        else:
            if self.controller.IsActorAdded(self.actor):
                self.GetValueForKey('controller.Renderer').RemoveActor(self.actor)
                pass
            pass
        return
    pass

    
class PipelineController(HasVtkObjectKeys, ObjectController):
    """Represent the VTK pipeline.
    """
    def __init__(self, delegate, profileController):
        ObjectController.__init__(self, delegate)
        
        self.profileController = profileController
        
        self.mode = 'view'
        
        self.DefineVtkObjectKey('StlMapper')
        self.DefineVtkObjectKey('StlActor')
        self.DefineVtkObjectKey('Locator')
        
        self.GetValueForKey('StlMapper.SetInputConnection')(
            profileController.GetValueForKey('StlReader.GetOutputPort')()
            )
        
        self.SeedPlacer = SeedPlacer(self)
        
        profileController.BindValue('SeedPoint.x',
                                    SeedCoordMapper(0, self.SeedPlacer))
        profileController.BindValue('SeedPoint.y',
                                    SeedCoordMapper(1, self.SeedPlacer))
        profileController.BindValue('SeedPoint.z',
                                    SeedCoordMapper(2, self.SeedPlacer))
        
        profileController.AddObserver('StlReader.Modified', self.HandleStlReaderModified)
        self.AddDependency('SeedPlaceButtonEnabled', 'mode')
        self.AddDependency('SeedPlaceButtonLabel', 'mode')
        return
    
    def IsActorAdded(self, actor):
        """Return whether the supplied argument is in the renderer's
        list of actors.
        """
        aList = self.delegate.Renderer.GetActors()
        iterator = aList.NewIterator()
        
        while not iterator.IsDoneWithTraversal():
            a = iterator.GetCurrentObject()
            if a is actor:
                return True
            iterator.GoToNextItem()
            continue
        return False
    
    @property
    def SeedPlaceButtonLabel(self):
        if self.mode == 'seed':
            return 'Finish'
        return 'Place'
    
    @property
    def SeedPlaceButtonEnabled(self):
        if self.mode == 'view' or self.mode == 'seed':
            return True
        return False
    
    def SeedPlaceClicked(self):
        if self.mode == 'view':
            self.SeedPlacer.On = True
            self.mode = 'seed'
            
        elif self.mode == 'seed':
            self.mode = 'view'
            
            pass
        
        return
    
    def SetInteractor(self, iact):
        self.delegate.SetInteractor(iact)
        
        self.viewInteractorStyle = vtkInteractorStyleTrackballCamera()
        iact.SetInteractorStyle(self.viewInteractorStyle)
        
        self.leftClickObserverId = iact.AddObserver(
            "LeftButtonPressEvent", self.HandleLeftClick
            )
        return
    
    def HandleStlReaderModified(self, change):
        self.GetValueForKey('Locator.SetDataSet')(
            change.obj.GetValueForKey('StlReader.GetOutput')()
            )
        self.GetValueForKey('Locator.BuildLocator')()
        self.delegate.ResetView()
        return

    def HandleLeftClick(self, obj, evt):
        """Callback to handle VTK click events in the RWI.
        """
        if self.mode == 'seed':
            self.HandleLeftButtonPressForSeedPlacement(obj, evt)
            pass
        return
    
    
    def HandleLeftButtonPressForSeedPlacement(self, obj, evt):
        mousePos = obj.GetEventPosition()
        didClickSurface, worldPos = self.MouseToWorld(mousePos)
        
        if didClickSurface:
            self.SeedPlacer.SetCentre(worldPos)
            
            # Want to abort further handling of this event, but
            # currently can't do this from Python. Grr.
            pass

        return
    
    def MouseToWorld(self, mousePos, worldPos=None):
        worldRefPos = [0.,0.,0.]
        if worldPos is None:
            worldPos = [0., 0., 0.]
            pass
        worldOrient = [0.,0.,0.,
                       0.,0.,0.,
                       0.,0.,0.]
        didClickSurface = self.GetValueForKey('SurfacePlacer').ComputeWorldPosition(
            self.GetValueForKey('Renderer'), mousePos, worldRefPos,
            worldPos, worldOrient)
        return didClickSurface, worldPos


    
    # def ChangeStlValidity(self, change):
    #     if self.GetValueForKey('HaveValidStlFile'):
    #         # We have a valid STL
    #         newStlFile = self.GetValueForKey('StlFile')
    #         if newStlFile != self._lastStlFile:
    #             # The new path is different from the last one
    #             self.locator.SetDataSet(self.GetValueForKey('StlReader').GetOutput())
    #             self.locator.BuildLocator()
    #             self._lastStlFile = newStlFile
    #         else:
    #             # Its the same file
    #             pass
    #     else:
    #         # STL is invalid, do nothin
    #         pass
    #     return
    
    # def IsActorAdded(self, actor):
    #     """Return whether the supplied argument is in the renderer's
    #     list of actors.
    #     """
    #     aList = self.renderer.GetActors()
    #     iterator = aList.NewIterator()
        
    #     while not iterator.IsDoneWithTraversal():
    #         a = iterator.GetCurrentObject()
    #         if a is actor:
    #             return True
    #         iterator.GoToNextItem()
    #         continue
    #     return False
    
    # def Show(self):
    #     if not self.IsActorAdded(self.stlActor):
    #         self.renderer.AddActor(self.stlActor)
    #         self.rwi.Update()
    #         pass
    #     return
    
    # def PlaceSeed(self):
    #     self.seeder.PlaceOn()
    #     return

    # def SetSeedPoint(self, pos):
    #     return self.seeder.SetSeedPoint(pos)
    
    pass

    # def GetActorReady(self, placer):
    #     # Update the centre
    #     if self.SeedPoint is not None:
    #         self.representation.SetCenter(self.SeedPoint)
    #         bb = placer.surfaceActor.GetBounds()
    #         # Set the radius to 1% of the average bounding box side
    #         self.representation.SetRadius(sum([bb[2*i+1]-bb[2*i] for i in range(3)])/300.)
    #         return self.actor
    #     return None

