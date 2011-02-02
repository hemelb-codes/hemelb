from vtk import vtkPolyDataMapper, vtkActor, \
     vtkSphereSource, vtkPolygonalSurfacePointPlacer
from wx.lib.pubsub import Publisher


import pdb

class SurfacePlacer(object):
    """Base class for dealing with the placing of entities in the scene.
    States of this object:

    enabled:

        The actor(s) representing the entity(ies) should be shown if and
        only if we have a valid position.

        placing:

            We should respond to clicks within the RWI and then call the
            HandlePlacementAtPosition method with the world coordinates if
            they were valid.

    notenabled:

        Shouldn't show the actor(s) and don't do anything at all.

    """
    def __init__(self, pipeline, surfaceActor):
        """We'll delegate unknown attributes to the pipeline.

        surfaceActor is the vtkActor representing the surface on which
        to place things.
        """
        self.pipeline = pipeline
        self.surfaceActor = surfaceActor

        self.entities = []

        self.enabled = False
        self.placing = False

        self.placer = vtkPolygonalSurfacePointPlacer()
        self.placer.AddProp(self.surfaceActor)

        return

    def On(self):
        if not self.enabled:
            self.enabled = True
            self.Show()

        return

    def Off(self):
        if self.enabled:
            self.PlaceOff()
            self.Hide()
            self.enabled = False
        return

    def Show(self):
        """Show the actors in the scene if we're enabled and and have
        a valid location.
        """
        if not self.enabled:
            return

        for entity in self.entities:
            actor = entity.GetActorReady(self)
            if actor is None: continue
            
            if not self.pipeline.IsActorAdded(actor):
                # Our actor isn't on stage, so add it
                self.pipeline.renderer.AddActor(actor)
                pass
            continue
        
        return

    def Hide(self):
        """Hide any of our on stage actors.
        """
        for entity in self.entities:
            if self.pipeline.IsActorAdded(entity.actor):
                # If the actor's on stage, remove it
                self.pipeline.renderer.RemoveActor(entity.actor)
                pass
            continue
        
        self.pipeline.rwi.Update()
        return

    def PlaceOn(self):
        """Enter placing mode, where we watch for clicks in the RWI.
        """
        if not self.enabled:
            self.On()

        self.Hide()

        if not self.placing:
            self.placing = True
            self.leftClickObserverId = self.pipeline.rwi.AddObserver(
                "LeftButtonPressEvent",self.HandleLeftClick
                )
        return

    def PlaceOff(self):
        """Finish placing mode, so remove our observer.
        """
        if self.placing:
            self.placing = False
            self.pipeline.rwi.RemoveObserver(self.leftClickObserverId)
            self.Show()
        return

    def __del__(self):
        """Turn ourself off before deletion.
        """
        self.Off()
        return

    def HandleLeftClick(self, obj, event):
        """Callback to handle VTK click events in the RWI.
        """
        mousePos = self.pipeline.rwi.GetEventPosition()
        worldRefPos = [0.,0.,0.]
        worldPos = [0., 0., 0.]
        worldOrient = [0.,0.,0.,
                       0.,0.,0.,
                       0.,0.,0.]

        didClickSurface = self.placer.ComputeWorldPosition(
            self.pipeline.renderer, mousePos, worldRefPos,
            worldPos, worldOrient)
        if didClickSurface:
            self.HandlePlacementAtPosition(worldPos)
            # Want to abort further handling of this event, but
            # currently can't do this from Python. Grr.
            pass

        return

    def HandlePlacementAtPosition(self, worldPos):        
        raise NotImplementedError(
            "Subclasses of SurfacePlacer must override the 'HandlePlacementAtPosition' method")

    pass

class SeedPlacer(SurfacePlacer):
    """Deal with the placing of the seed point.
    States of this object:

    enabled:
        The seed actor should be shown if and only if we have a valid position.

        placing:
            We should respond to clicks within the RWI and emit a
            'view.SeedPoint' message if we can get a valid seed point
            from it

    notenabled:
        Shouldn't show the seed actor and don't do anything at all.

    """

    def __init__(self, pipeline, surfaceActor):
        SurfacePlacer.__init__(self, pipeline, surfaceActor)
        self.entities.append(SeedEntity())
        self.seed = self.entities[0]
        return
    
    def SetSeedPoint(self, pos):
        # Keep our copy of the point
        self.seed.SetSeedPoint(pos)
        
        if self.enabled:
            # If we're enabled make sure the seed is shown (Show()
            # checks that the pos is valid so this is always safe)
            self.Show()

    def HandlePlacementAtPosition(self, worldPos):
        Publisher.sendMessage('view.SeedPoint',
                              data=worldPos)
        self.PlaceOff()

    pass

class IoletPlacer(SurfacePlacer):
    def __init__(self, pipeline, surfaceActor):
        self.ioletWidgets = {}
        return
    
    def ComputeWidgetBoundsAtPosition(self, worldPos):
        return (worldPos[0]-5., worldPos[0]+5.,
                worldPos[1]-5., worldPos[1]+5.,
                worldPos[2]-5., worldPos[2]+5.)

    def ComputePlaneCenter(self, worldPos):
        return worldPos

    def ComputePlaneNormal(self, worldPos):
        return (0., 0., 1.)
    
    def HandlePlacementAtPosition(self, worldPos):
        newWidget = self.ioletWidgets[-1] = vtkPlaneWidget()
        newWidget.PlaceWidget(self.ComputeWidgetBounds(worldPos))
        newWidget.SetCenter(self.ComputePlaneCenter(worldPos))
        newWidget.SetNormal(self.ComputePlaneNormal(worldPos))
        
    pass

class Entity(object):
    pass

class SeedEntity(Entity):

    def __init__(self):
        self.SeedPoint = None
        
        self.representation = vtkSphereSource()
        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.representation.GetOutputPort())
        
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        # Make it red
        self.actor.GetProperty().SetColor(1,0,0)
        return

    def SetSeedPoint(self, pos):
        self.SeedPoint = pos
        return
    
    def GetActorReady(self, placer):
        # Update the centre
        if self.SeedPoint is not None:
            self.representation.SetCenter(self.SeedPoint)
            bb = placer.surfaceActor.GetBounds()
            # Set the radius to 1% of the average bounding box side
            self.representation.SetRadius(sum([bb[2*i+1]-bb[2*i] for i in range(3)])/300.)
            return self.actor
        return None
    pass
