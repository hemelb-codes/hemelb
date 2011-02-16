import operator
from vtk import vtkInteractorStyleTrackballCamera

from ..Bindings.ObjectController import ObjectController
from ..Bindings.ListController import ListContentsDestMapper
from ..Bindings.VtkObject import HasVtkObjectKeys
from ..Bindings.Mappers import Mapper, SimpleObservingMapper
from ..Bindings.Translators import UnitTranslator

from .PlacedIoletController import HasPlacedIoletListKeys
import pdb


class SeedCoordMapper(Mapper):
    """Custom bindings mapper for the SeedPosition attribute to the
    'Center' of the VTK pipeline's representation of the seed
    location.
    """
    def __init__(self, i, placedSeed, translator=UnitTranslator()):
        Mapper.__init__(self, translator)
        
        self.PlacedSeed = placedSeed
        self.i = i
        self.obsId = None
        return
    
    def _Get(self):
        return self.PlacedSeed.GetCentre()[self.i]
    
    def _Set(self, val):
        if val != val:
            # It's a nan
            self.PlacedSeed.Enabled = False
            pass
        
        old = list(self.PlacedSeed.GetCentre())
        old[self.i] = val
        self.PlacedSeed.SetCentre(old)

        allFinite = reduce(operator.and_, map(lambda x: x==x, old))
        if allFinite:
            self.PlacedSeed.Enabled = True
            pass
        return

    def Observe(self):
        if self.obsId is not None:
            self.Unobserve()
            pass
        
        self.obsId = self.PlacedSeed.representation.AddObserver('ModifiedEvent', self.HandleUpdate)
        return
    def Unobserve(self):
        if self.obsId is None:
            return
        
        self.PlacedSeed.representation.RemoveObserver(self.obsId)
        self.obsId = None
        return
    
    pass


class PipelineController(HasVtkObjectKeys, HasPlacedIoletListKeys, ObjectController):
    """Represent the VTK pipeline.
    """
    def __init__(self, delegate, profileController):
        ObjectController.__init__(self, delegate)
        
        self.profileController = profileController
        
        self.mode = 'view'
        
        self.DefineVtkObjectKey('StlMapper')
        self.DefineVtkObjectKey('StlActor')
        self.DefineVtkObjectKey('Locator')
        self.DefinePlacedIoletListKey('PlacedIolets')
        
        profileController.BindValue('Iolets',
                                    ListContentsDestMapper(self.delegate.PlacedIolets,
                                                           translator=self.PlacedIolets.translator))
        profileController.BindValue('Iolets.SelectedIndex',
                                    SimpleObservingMapper(self.PlacedIolets, 'SelectedIndex'))
        
        self.GetValueForKey('StlMapper.SetInputConnection')(
            profileController.GetValueForKey('StlReader.GetOutputPort')()
            )
        
        # self.PlacedSeed = PlacedSeed(self)
        
        profileController.BindValue('SeedPoint.x',
                                    SeedCoordMapper(0, self.GetValueForKey('PlacedSeed')))
        profileController.BindValue('SeedPoint.y',
                                    SeedCoordMapper(1, self.GetValueForKey('PlacedSeed')))
        profileController.BindValue('SeedPoint.z',
                                    SeedCoordMapper(2, self.GetValueForKey('PlacedSeed')))

        
        profileController.AddObserver('StlReader.Modified', self.HandleStlReaderModified)
        self.AddDependency('SeedPlaceButtonEnabled', 'mode')
        self.AddDependency('SeedPlaceButtonLabel', 'mode')
        self.AddDependency('IoletPlaceButtonEnabled', 'mode')
        self.AddDependency('IoletPlaceButtonLabel', 'mode')
        return
        
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
            self.SetValueForKey('PlacedSeed.Enabled', True)
            self.mode = 'seed'
            
        elif self.mode == 'seed':
            self.mode = 'view'
            
            pass
        
        return
    
    @property
    def IoletPlaceButtonLabel(self):
        if self.mode == 'iolet':
            return 'Finish'
        return 'Place'
    
    @property
    def IoletPlaceButtonEnabled(self):
        if self.mode =='view' or self.mode == 'iolet':
            return True
        return False
    
    def IoletPlaceClicked(self):
        if self.mode == 'view':
            self.SetValueForKey('PlacedSeed.Enabled', True)
            self.mode = 'iolet'
            
        elif self.mode == 'iolet':
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
        elif self.mode == 'iolet':
            self.HandleLeftButtonPressForIoletPlacement(obj, evt)
            pass
        
        return
    
    
    def HandleLeftButtonPressForSeedPlacement(self, obj, evt):
        mousePos = obj.GetEventPosition()
        didClickSurface, worldPos = self.MouseToWorld(mousePos)
        
        if didClickSurface:
            self.SetValueForKey('PlacedSeed.Centre', worldPos)
            
            # Want to abort further handling of this event, but
            # currently can't do this from Python. Grr.
            pass

        return
    
    def HandleLeftButtonPressForIoletPlacement(self, obj, evt):
        mousePos = obj.GetEventPosition()
        didClickSurface, worldPos = self.MouseToWorld(mousePos)
        
        if didClickSurface:
            self.SetValueForKey('PlacedIolets.Selection.Centre', worldPos)
            self.SetValueForKey('PlacedIolets.Selection.Enabled', True)
#            self.SetValueForKey('PlacedIolets.Selection.Normal', (0.,0.,1.))
#            self.SetValueForKey('PlacedIolets.Selection.Radius', 1.)
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

