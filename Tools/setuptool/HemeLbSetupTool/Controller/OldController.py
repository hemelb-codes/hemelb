from contextlib import contextmanager
from functools import wraps

import wx
from wx.lib.pubsub import Publisher
from wx.lib.masked import EVT_NUM

from Notifier import Notifier

import pdb

def MessageHandler(method):
    """Decorator that wraps instance methods. It checks whether the
    topic for the supplied message is ignored and only calls the
    method if it is not.
    """
    @wraps(method)
    def methodwrapper(self, message):
        if not self.shouldIgnoreMessage(message):
            return method(self, message)
        return
    # return methodwrapper
    return method

def EventHandler(method):
    """Decorator that wraps instance methods. It checks whether the
    event triggers is ignored and only calls the method if it is not.
    """
    def methodwrapper(self, event):
        if not self.shouldIgnoreEvent(event):
            return method(self, event)
        event.Skip()
        return
    # return methodwrapper

    # For the time being don't bother wrapping as I haven't figured
    # out the best way to ignore WX events.
    return method

class Controller(Notifier):
    """The App's controller. Stitches together the View and Model.
    """
    
    def __init__(self, model, view):
        self.model = model
        self.view = view
        
        self._ignoredTopics = set()
        self._ignoredEvents = {}

        self._placeMode = 'off'
        pass
    
    def BindToUiEvents(self):
        """Bind our EventHandler methods to UI events.
        """
        # Profile buttons:
        view = self.view.main
        # STL file selection events
        view.main.Bind(wx.EVT_BUTTON, self.OnChooseStl,
                       source=self.view.inputStlButton)
        # IOLet events
        view.Bind(EVT_BUTTON, self.OnViewAddInlet,
                  source=self.view.addInletButton)
        view.Bind(EVT_BUTTON, self.OnViewAddOutlet,
                  source=self.view.addOutletButton)
        view.Bind(EVT_BUTTON, self.OnViewRemoveIolet,
                  source=self.view.removeIoletButton)
        
        # Voxel size events
        view.Bind(EVT_NUM, self.OnViewVoxelSizeEdited,
                  source=self.view.voxelSizeField)
        view.Bind(wx.EVT_BUTTON, self.model.ResetVoxelSize,
                  source=self.view.voxelResetButton)
        # Seed events
        view.Bind(wx.EVT_BUTTON, self.OnViewPlaceSeed,
                  source=self.view.seedPlaceButton)
        
        seedVec = self.view.seedVector
        view.Bind(EVT_NUM, self.OnViewSeedPointEdited,
                  source=seedVec.x)
        view.Bind(EVT_NUM, self.OnViewSeedPointEdited,
                  source=seedVec.y)
        view.Bind(EVT_NUM, self.OnViewSeedPointEdited,
                  source=seedVec.z)
        
        # Get the SeedPoint from the VTK code that figures it out from
        # the UI event
        Publisher.subscribe(self.OnViewSeedPointSent,
                            'view.SeedPoint')
        
        # Output events
        view.Bind(wx.EVT_BUTTON, self.OnChooseOutputXmlFile,
                  source=self.view.xmlChooseButton)
        view.Bind(wx.EVT_TEXT, self.OnViewOutputXmlFileEdited,
                  source=self.view.xmlField)
        view.Bind(wx.EVT_BUTTON, self.OnChooseOutputConfigFile,
                  source=self.view.configChooseButton)
        view.Bind(wx.EVT_TEXT, self.OnViewOutputConfigFileEdited,
                  source=self.view.configField)
        
        # View controls
        view.Bind(wx.EVT_BUTTON, self.OnDebug,
                  source=self.view.debugButton)
        view.Bind(wx.EVT_BUTTON, self.OnResetView,
                  source=self.view.resetViewButton)
        return
    
    def SubscribeToModelTopics(self):
        """Subscribe our MessageHandler methods to changes in the
        model data.
        """
        # Publisher.subscribe(self.OnModelStlFileChanged,
        #                     self.model.SelfTopic('StlFileChanged'))
        # Publisher.subscribe(self.OnModelVoxelSizeChanged,
        #                     self.model.SelfTopic('VoxelSizeChanged'))
        # Publisher.subscribe(self.OnModelSeedPointChanged,
        #                     self.model.SelfTopic('SeedPointChanged'))
        # Publisher.subscribe(self.OnModelOutputXmlFileChanged,
        #                     self.model.SelfTopic('OutputXmlFileChanged'))
        # Publisher.subscribe(self.OnModelOutputConfigFileChanged,
        #                     self.model.SelfTopic('OutputConfigFileChanged'))
        self.model.addObserver('StlFile', self.OnModelStlFileChanged)
        self.model.addObserver('VoxelSize', self.OnModelVoxelSizeChanged)
        self.model.addObserver('SeedPoint', self.OnModelSeedPointChanged)
        self.model.addObserver('OutputXmlFile', self.OnModelOutputXmlFileChanged)
        self.model.addObserver('OutputConfigFile', self.OnModelOutputConfigFileChanged)
        return
    
    def shouldIgnoreEvent(self, event):
        """Work out if we should ignore this event.
        """
        return False
    
    def shouldIgnoreMessage(self, msg):
        """Indicate if we should ignore the delivered message.
        """
        if msg.topic in self._ignoredTopics:
            return True
        return False
    
    def IgnoreTopic(self, topic):
        """Mark a topic as to be ignored.
        """
        self._ignoredTopics.add(topic)
        return
    
    def UnignoreTopic(self, topic):
        """If the topic is ignored, remove it from the list of ignored
        ones.
        """
        self._ignoredTopics.discard(topic)
        return

    @contextmanager
    def IgnoredTopicContext(self, topic):
        """Context manager to temporarily ignore messages sent to a
        topic. Use as
        
        > with controller.IgnoredTopicContext('yadda'):
        >     doSomething()
        
        """
        self.IgnoreTopic(topic)
        yield
        self.UnignoreTopic(topic)
        return
    
    def InitWindow(self):
        """Kick off the app.
        """
        # Create the windows
        frame = self.view.Create(self)
        
        self.view.pipeline.Start(self.model.stlReader.GetOutputPort())
        self.BindToUiEvents()
        self.SubscribeToModelTopics()
        
        # Show it
        frame.Show()
        frame.Maximize()
        
        return frame

    ##############################
    # Event and Message handlers #
    ##############################
    @EventHandler
    def OnDebug(self, event):
        """Drop into the debugger.
        """
        pdb.set_trace()
        return

    # View Controls
    @EventHandler
    def OnResetView(self, event):
        """Tell the view to reset.
        """
        # Trigger the custom RenderWindowInteractor's Reset Method
        self.view.rwi.ResetView()
        return

    # STL file selection events
    @EventHandler
    def OnChooseStl(self, event):
        dialog = wx.FileDialog(None, style=wx.FD_OPEN, wildcard='*.stl')

        if dialog.ShowModal() == wx.ID_OK:
            self.model.StlFile = dialog.GetPath()
            
            self.view.pipeline.Show()
            self.view.rwi.ResetView()
            pass
        
        dialog.Destroy()
        return

    @MessageHandler
    def OnModelStlFileChanged(self, msg):
        self.view.inputStlField.SetValue(msg.data)
        return

    # IOLets events

    # Voxel size events
    @MessageHandler
    def OnModelVoxelSizeChanged(self, msg):
        newSize = msg.data
        self.view.voxelSizeField.ChangeValue(newSize)
        return
    
    @EventHandler
    def OnViewVoxelSizeEdited(self, event):
        # with self.IgnoredTopicContext(self.model.SelfTopic('VoxelSizeChanged')):
        self.model.VoxelSize = self.view.voxelSizeField.GetValue()
        # pass
        
        return

    # Seed placement events
    @EventHandler
    def OnViewPlaceSeed(self, event):
        if self._placeMode == 'off':
            self._placeMode = 'seed'
            self.view.pipeline.seeder.PlaceOn()
            self.view.seedPlaceButton.SetLabel('Cancel')
        elif self._placeMode == 'seed':
            self.EndSeedMode()
        return
    
    def EndSeedMode(self):
        self._placeMode = 'off'
        self.view.pipeline.seeder.PlaceOff()
        self.view.seedPlaceButton.SetLabel('Place')
        return
    
    @MessageHandler
    def OnViewSeedPointSent(self, msg):
        self.EndSeedMode()
        self.model.SeedPoint = msg.data
        return

    @MessageHandler
    def OnModelSeedPointChanged(self, msg):
        if msg.data is None:
            pos = (0.,0.,0.)
        else:
            pos = msg.data
            
        seedVec = self.view.seedVector
        seedVec.x.ChangeValue(pos[0])
        seedVec.y.ChangeValue(pos[1])
        seedVec.z.ChangeValue(pos[2])
        
        self.view.pipeline.seeder.SetSeedPoint(msg.data)
        return
    
    @EventHandler
    def OnViewSeedPointEdited(self, event):
        seedVec = self.view.seedVector
        pos = [seedVec.x.GetValue(),
               seedVec.y.GetValue(),
               seedVec.z.GetValue()]
        self.model.SeedPoint = pos
        return
    
    # Output file events
    @EventHandler
    def OnChooseOutputConfigFile(self, event):
        dialog = wx.FileDialog(None,
                               style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT,
                               wildcard='*.dat')

        if dialog.ShowModal() == wx.ID_OK:
            self.model.OutputConfigFile = dialog.GetPath()
            pass
        
        dialog.Destroy()
        return
    
    @EventHandler
    def OnViewOutputConfigFileEdited(self, event):
        # with self.IgnoredTopicContext(self.model.SelfTopic('OutputConfigFileChanged')):
        self.model.OutputConfigFile = self.view.configField.GetValue()
        # pass
        
        return
    
    @MessageHandler
    def OnModelOutputConfigFileChanged(self, msg):
        self.view.configField.SetValue(msg.data)
        return
    
    @EventHandler
    def OnChooseOutputXmlFile(self, event):
        dialog = wx.FileDialog(None,
                               style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT,
                               wildcard='*.xml')

        if dialog.ShowModal() == wx.ID_OK:
            self.model.OutputXmlFile = dialog.GetPath()
            pass
        
        dialog.Destroy()
        return
    
    @EventHandler
    def OnViewOutputXmlFileEdited(self, event):
        # with self.IgnoredTopicContext(self.model.SelfTopic('OutputXmlFileChanged')):
        self.model.OutputXmlFile = self.view.xmlField.GetValue()
        # pass
        
        return

    @MessageHandler
    def OnModelOutputXmlFileChanged(self, msg):
        self.view.xmlField.SetValue(msg.data)
        return

    pass
