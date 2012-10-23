# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from vtk import vtkPolyDataMapper, vtkActor, vtkModifiedBSPTree, \
     vtkPolygonalSurfacePointPlacer, vtkRenderer, vtkAnnotatedCubeActor, \
     vtkAxesActor, vtkTextProperty, vtkPropAssembly, \
     vtkOrientationMarkerWidget, vtkSphereSource

from ..Util.Observer import Observable

from .PlacedIolet import PlacedIoletList

import pdb

class Pipeline(Observable):
    def __init__(self):
        # 1 mm is probably about right.
        # TODO: force this to be recalculated when the StlReader updates
        self.WidgetSize = 1
        
        self.SurfaceMapper = vtkPolyDataMapper()
        self.SurfaceActor = vtkActor()
        self.SurfaceActor.SetMapper(self.SurfaceMapper)
        # self.seeder = SeedPlacer(self, self.stlActor)
        self.Locator = vtkModifiedBSPTree()
        
        self.SurfacePlacer = vtkPolygonalSurfacePointPlacer()
        self.SurfacePlacer.AddProp(self.SurfaceActor)

        self.Renderer = vtkRenderer()
        self.Renderer.AddActor(self.SurfaceActor)
        
        self.CreateMarker()

        self.PlacedSeed = PlacedSeed(self)
        self.PlacedSeed.AddObserver('Enabled', self.HandlePlacedItemEnabledChange)
        
        self.PlacedIolets = PlacedIoletList()
        self.PlacedIolets.SetItemEnabledChangeHandler(self.HandlePlacedItemEnabledChange)
        
        return
            
    def ResetView(self):
        """Reset the view on the current scene.
        """
        self.Renderer.ResetCamera()
        return
    
    def SetViewX(self):
        cam = self.Renderer.GetActiveCamera()
        focus = cam.GetFocalPoint()
        dist = cam.GetDistance()
        cam.SetPosition(focus[0]+dist,focus[1], focus[2])
        cam.SetViewUp(0,0,1)
        return
    def SetViewY(self):
        cam = self.Renderer.GetActiveCamera()
        focus = cam.GetFocalPoint()
        dist = cam.GetDistance()
        cam.SetPosition(focus[0],focus[1]+dist, focus[2])
        cam.SetViewUp(0,0,1)
        return
    def SetViewZ(self):
        cam = self.Renderer.GetActiveCamera()
        focus = cam.GetFocalPoint()
        dist = cam.GetDistance()
        cam.SetPosition(focus[0],focus[1], focus[2]+dist)
        cam.SetViewUp(1,0,0)
        return

    def SetInteractor(self, iact):
        self.Interactor = iact
        self.SetMarkerInteractor(iact)
        self.PlacedIolets.SetInteractor(iact)
        
        return
    
    def SetInteractorForObject(self, obj):
        obj.SetInteractor(self.Interactor)
        return
    
    def SetMarkerInteractor(self, iact):
        self.OrientationMarker.SetInteractor(iact)
        self.OrientationMarker.SetEnabled(1)
        self.OrientationMarker.InteractiveOff()
        return
    
    def CreateMarker(self):
        # Create a composite orientation marker using
        # vtkAnnotatedCubeActor and vtkAxesActor.
        #
        cube = vtkAnnotatedCubeActor()
        cube.SetXPlusFaceText("R")
        cube.SetXMinusFaceText("L")
        cube.SetYPlusFaceText("A")
        cube.SetYMinusFaceText("P")
        cube.SetZPlusFaceText("I")
        cube.SetZMinusFaceText("S")
        cube.SetXFaceTextRotation(0)
        cube.SetYFaceTextRotation(0)
        cube.SetZFaceTextRotation(-90)
        cube.SetFaceTextScale(0.65)

        prop = cube.GetCubeProperty()
        prop.SetColor(0.5, 1, 1)

        prop = cube.GetTextEdgesProperty()
        prop.SetLineWidth(1)
        prop.SetDiffuse(0)
        prop.SetAmbient(1)
        prop.SetColor(0.18, 0.28, 0.23)

        for axis, colour in (('X', (0,0,1)),
                             ('Y', (0,1,0)),
                             ('Z', (1,0,0))):
            for orient in ('Plus', 'Minus'):
                prop = getattr(cube, 'Get'+axis+orient+'FaceProperty')()
                prop.SetColor(*colour)
                prop.SetInterpolationToFlat()
                continue
            continue

        axes = vtkAxesActor()
        axes.SetShaftTypeToCylinder()
        axes.SetXAxisLabelText("x")
        axes.SetYAxisLabelText("y")
        axes.SetZAxisLabelText("z")
        axes.SetTotalLength(1.5, 1.5, 1.5)

        tpropx = vtkTextProperty()
        tpropx.ItalicOn()
        tpropx.ShadowOn()
        tpropx.SetFontFamilyToTimes()
        axes.GetXAxisCaptionActor2D().SetCaptionTextProperty(tpropx)
        tpropy = vtkTextProperty()
        tpropy.ShallowCopy(tpropx)
        axes.GetYAxisCaptionActor2D().SetCaptionTextProperty(tpropy)
        tpropz = vtkTextProperty()
        tpropz.ShallowCopy(tpropx)
        axes.GetZAxisCaptionActor2D().SetCaptionTextProperty(tpropz)

        # Combine the two actors into one with vtkPropAssembly ...
        #
        assembly = vtkPropAssembly()
        assembly.AddPart(axes)
        assembly.AddPart(cube)

        # Add the composite marker to the widget.  The widget
        # should be kept in non-interactive mode and the aspect
        # ratio of the render window taken into account explicitly, 
        # since the widget currently does not take this into 
        # account in a multi-renderer environment.
        # 
        marker = vtkOrientationMarkerWidget()
        marker.SetOutlineColor(0.93, 0.57, 0.13)
        marker.SetOrientationMarker(assembly)
        marker.SetViewport(0.0, 0.0, 0.15, 0.3)
        self.OrientationMarker = marker
        return

    def IsActorAdded(self, actor):
        """Return whether the supplied argument is in the renderer's
        list of actors.
        """
        aList = self.Renderer.GetActors()
        iterator = aList.NewIterator()
        
        while not iterator.IsDoneWithTraversal():
            a = iterator.GetCurrentObject()
            if a is actor:
                return True
            iterator.GoToNextItem()
            continue
        return False

    def HandlePlacedItemEnabledChange(self, change):
        item = change.obj
        if item.Enabled:
            if not self.IsActorAdded(item.actor):
                self.AddObserver('WidgetSize', item.HandleWidgetSizeChange)
                self.Renderer.AddActor(item.actor)
                pass
        else:
            if self.IsActorAdded(item.actor):
                self.Renderer.RemoveActor(item.actor)
                self.RemoveObserver('WidgetSize', item.HandleWidgetSizeChange)
                pass
            pass
        return

    pass

class PlacedSeed(Observable):
    def __init__(self, pipeline, colour=(0,0,1)):
        self.representation = vtkSphereSource()
        self.representation.SetRadius(pipeline.WidgetSize)
#        pipeline.AddObserver('WidgetSize', self.HandleWidgetSizeChange)
        
        self.mapper = vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.representation.GetOutputPort())
        
        self.actor = vtkActor()
        self.actor.SetMapper(self.mapper)
        # Make it blue
        self.actor.GetProperty().SetColor(colour)

        self.Enabled = False
        
        return
    
    def HandleWidgetSizeChange(self, change):
        pro = change.obj
        self.representation.SetRadius(pro.WidgetSize)
        return
    
    def SetCentre(self, centre):
        self.representation.SetCenter(centre)
        return
    def GetCentre(self):
        return self.representation.GetCenter()
    Centre = property(GetCentre, SetCentre)
    
    pass

