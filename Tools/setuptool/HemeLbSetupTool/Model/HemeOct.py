import Oct
import numpy as np
import os.path
import vtk

class NodeBase(Oct.NodeBase):
    """Base class for Heme-specific stuff.
    """
    def ToVtk(self):
        if self.levels > 1:
            if self.children is None:
                return None
            
            kids = []
            for child in self.children.flat:
                if child is not None:
                    # child must be created
                    kid = child.ToVtk()
                    if kid is not None:
                        kids.append(kid)
                        pass
                    pass
                continue
            
            if len(kids) == 0:
                return None
            
            myDS = vtk.vtkMultiBlockDataSet()
            myDS.SetNumberOfBlocks(len(kids))
            for i, kid in enumerate(kids):
                myDS.SetBlock(i, kid)
            return myDS
        elif self.levels == 1:
            # children are leaves, I therefore contain their data
            if self.children is None:
                return None
            
            myDS = vtk.vtkUniformGrid()
            extent = np.concatenate((self.offset,
                                     self.offset + 1)).reshape((2,3)).transpose()
            myDS.SetExtent(extent.flatten())
            myDS.SetScalarTypeToUnsignedChar()
            pd = myDS.GetPointData()
            data = vtk.vtkUnsignedCharArray()
            data.SetNumberOfTuples(8)
            
            strides = np.array((1,2,4))
            for ijk, child in np.ndenumerate(self.children):
                i = np.sum(ijk * strides)
                data.SetTuple1(i, child is None)
            return myDS
        else:
            # am leaf  - should never be here
            raise ValueError("Should not convert a leaf node to VTK")
            pass
        
    def Write(self, streamMap):
        if self.fluid_count == 0:
            return
        
        # Each node in the tree has 8 flags for the presence of children, child stream position and a fluid/solid flag 
        child_state = np.uint8(0)
        child_stream_pos = 0
        if self.levels:
            child_stream_pos = streamMap[self.levels - 1].tell()
            if self.children is not None:
                for i, child in enumerate(self.children.flat):
                    if child is not None and child.fluid_count > 0:
                        child_state ^= (1 << i)
                        child.Write(streamMap)
        
        s = streamMap[self.levels]
        line = '{} {} {}\n'.format(child_state, child_stream_pos, self.fluid_count)
        s.write(line)

class Node(NodeBase, Oct.Node):
    """Standard node in the tree.
    """
    pass

class Tree(Oct.Tree):
    NodeClass = Node
    def ToVtk(self, name):
        ds = self.root.ToVtk()
        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer.SetFileName(name)
        writer.SetInput(ds)
        writer.Write()
        
    def Write(self, name):
        os.mkdir(name)
        streams = []
        names = []
        for i in xrange(self.levels+1):
            fn = os.path.join(name, '{:02d}'.format(i))
            names.append(fn)
            streams.append(NodeStream(fn))
            continue
        self.root.Write(streams)
        for f in streams:
            f.close()
        return
    pass

class UniformNode(NodeBase):
    """Node of constant fluid/solid type.
    """
    def GetNode(self, level, index, relative=False):
        if level == 0:
            assert np.all(index == (0, 0, 0)), "Indexing error"
            return self
        raise ValueError("Should never try to get a subnode of a uniform block!")
    
    def IterDepthFirst(self, *levels):
        if len(levels) == 0:
            bot_level = 0
            top_level = self.levels
        elif len(levels) == 1:
            bot_level = levels[0]
            top_level = self.levels
        elif len(levels) == 2:
            bot_level, top_level = levels
        else:
            raise ValueError("Max two arguments!")
        if bot_level <= self.levels and top_level >= self.levels:
            yield self
            
    pass

class SolidNode(UniformNode):
    fluid_count = 0
    def __eq__(self, other):
        if isinstance(other, SolidNode):
            return super(SolidNode, self).__eq__(other)
    pass
class FluidNode(UniformNode):
    @property
    def fluid_count(self):
        return 8**self.levels
    def __eq__(self, other):
        if isinstance(other, FluidNode):
            return super(FluidNode, self).__eq__(other)
    pass

class NodeStream(file):
    """Wrapper around __builtins__.file to track the number of lines written for
    reference by the Write() methods above.
    """
    def __init__(self, fn):
        file.__init__(self, fn, 'w')
        self.n = 0
    def write(self, data):
        file.write(self, data)
        self.n += 1
    def tell(self):
        return self.n
    pass
