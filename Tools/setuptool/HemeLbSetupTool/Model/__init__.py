from math import sqrt

from vtk import vtkSTLReader, vtkPolyDataMapper, vtkActor, vtkProgrammableFilter, vtkPoints, vtkFloatArray, vtkPolyData
# from wx.lib.pubsub import Publisher

# from Notifier import Notifier
from Observer import Observable, ObservableList
import pdb

# class NotifyingList(Notifier, list):
#     def __new__(cls, *args, **kwargs):
#         """Do required initialisation to ensure that we insert the
#         private attribute that stores the topics that will be messaged
#         on append.
#         """
#         new = Notifier.__new__(cls, *args, **kwargs)
#         object.__setattr__(new, '_topicsToNotifyOnSet', [])
#         return new
    
#     def addNotifiedTopicOnSet(self, topic):
#         """Append a topic to the list to be notified when an element
#         changes. Returns a boolean indicating whether the topic was
#         added or not.
#         """
#         assert hasattr(self, attr)
        
#         topic = self._canonicaliseTopic(topic)
        
#         if topic in self._topicsToNotifyOnSet:
#             return False
#         else:
#             attrTopics.append(topic)
#             return True
#         return

#     def removeNotifiedTopic(self, attr, topic):
#         """Remove a topic from the list to be notified when an element
#         changes. Returns a boolean indicating whether the topic was
#         removed or if it wasn't present.
#         """
#         assert hasattr(self, attr)
#         topic = self._canonicaliseTopic(topic)
        
#         try:
#             self._topicsToNotifyOnSet.remove(topic)
#             return True
#         except ValueError:
#             return False
#         return
 
#     def __setitem__(self, i, value):
#         """Set the attribute and notify topics if the value's changed.
#         """
#         try:
#             # If the value hasn't changed, don't need to notify 
#             if value == self[i]:
#                 return
#         except KeyError:
#             # But if this is the first assignment to the index, the
#             # above will be an error; one that we can swallow.
#             pass
        
#         for topic in self._topicsToNotifyOnSet:
#             # Send the notifications to the topics
#             Publisher.sendMessage(topic, data=value)
#             continue
#         # Set attribute and return
#         return object.__setattr__(self, name, value)
    
#     pass

class Profile(Observable):
    _Args = {'StlFile': None,
             'Iolets': ObservableList(),
             'VoxelSize': None,
             'SeedPoint': None,
             'OutputConfigFile': None,
             'OutputXmlFile': None}
    
    def __init__(self, **kwargs):
        it = Profile._Args.iteritems()
        for a, default in it:
            setattr(self, a,
                    kwargs.pop(a, default))
            continue
        
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)
        
        self.stlReader = vtkSTLReader()
        
        self.sider = AverageSideLengthCalculator()
        self.sider.SetInputConnection(self.stlReader.GetOutputPort())
        
        # # Set up notifications for all key properties
        # for attr in self._Args:
        #     self.addNotifiedTopic(attr,
        #                           self.SelfTopic(attr+'Changed'))
        #     continue
        
        # #  Subscribe our method to StlFile changes
        # Publisher.subscribe(self.OnStlFileChanged, 
        #                     topic=self.SelfTopic('StlFileChanged'))

        self.addObserver('StlFile', self.OnStlFileChanged)
        return
    
    # def OnStlFileChanged(self, msg):
    #     self.stlReader.SetFileName(msg.data)
    #     self.VoxelSize = self.sider.GetOutputValue()
    #     return
    def OnStlFileChanged(self, change):
        self.stlReader.SetFileName(self.StlFile)
        self.VoxelSize = self.sider.GetOutputValue()
        return
    
    @classmethod
    def NewFromFile(cls, filename):
        new = cls()
        
    def Save(self, filename):
        return

    def Generate(self):
        return
    
    # def SetStlFile(self):
    #     Publisher.sendMessage('setuptool.Model.Profile.SetStlFile', data=self)
    #     return
    
    # def GetStlFile(self, stl):
    #     return
    
    # def AddIolet(self, iolet):
    #     return
    
    # def RemoveIolet(self, id):
    #     return

    def ResetVoxelSize(self, ignored=None):
        self.VoxelSize = self.sider.GetOutputValue()
        return
    
    pass

class AverageSideLengthCalculator(vtkProgrammableFilter):
    def __init__(self):
        self.SetExecuteMethod(self.Execute)
        
    def Execute(self, *args):
        polydata = self.GetPolyDataInput()
        nTris = polydata.GetNumberOfCells()
        
        totalPerim = 0.
        for i in xrange(nTris):
            tri = polydata.GetCell(i)
            perim = 0.
            for j in range(3):
                perim += sqrt(tri.GetEdge(j).GetLength2())
                continue
            totalPerim += perim
            continue
        
        aveSide = totalPerim / (3*nTris)

        p = vtkPoints()
        p.InsertPoint(0, 0.,0.,0.)
        
        val = vtkFloatArray()
        val.InsertValue(0, aveSide)
        
        out = vtkPolyData()
        out.SetPoints(p)
        out.GetPointData().SetScalars(val)
        
        self.GetPolyDataOutput().ShallowCopy(out)
        return

    def GetOutputValue(self):
        self.Update()
        vals = self.GetOutput().GetPointData().GetScalars()
        if vals is None:
            return None
        return vals.GetValue(0)
    
    pass

class Iolet(Observable):
    _Args = {'Centre': None,
             'Normal': [],
             'Radius': None}
    
    def __init__(self, **kwargs):
        it = Iolet._Args.iteritems()
        for a, default in it:
            setattr(self, a,
                    kwargs.pop(a, default))
            continue
        
        for k in kwargs:
            raise TypeError("__init__() got an unexpected keyword argument '%'" % k)
        
        # Set up notifications for all key properties
        for attr in self._Args:
            self.addNotifiedTopic(attr,
                                  self.SelfTopic(attr+'Changed'))
            continue
        return
    
class Inlet(Iolet):
    pass

class Outlet(Iolet):
    pass
