import numpy as np

class NodeBase(object):
    """Base class for Octree nodes.
    """
    shape = (2,2,2)
    def __init__(self, maxLevels, offset):
        self.levels = maxLevels
        self.offset = np.array(offset)
        self.children = None
        self.rel_offset = self.offset >> self.levels
        return
    
    def __eq__(self, other):
        if not isinstance(other, NodeBase):
            return False
        
        if self.levels != other.levels:
            return False
        if not np.all(self.offset == other.offset):
            return False
        return True
    
    def __ne__(self, other):
        return not (self == other)
    pass

class NodeError(LookupError):
    def __init__(self, parent, child_index):
        self.parent = parent
        self.child_index = child_index
    @property
    def message(self):
        return "Node at level {} index {} has no child {}".format(self.parent.levels, self.parent.offset, self.child_index)
    pass

class Node(NodeBase):
    
    def GetNode(self, level, index, relative=False, create=False):
        """Get the node at given level with a given index.
        
        The relative flag specifies whether we are counting from the top or the bottom.
        False (the default) counts from the bottom, so level = 0 => get leaf nodes
        True counts from the top, so level = 0 => return self
        
        Index must be a numpy array of shape == (3,) dtype == int
        """
        if relative:
            rel_level = level
            abs_level = self.levels - level
            rel_index = index
            abs_index = (index << abs_level) + self.offset
        else:
            abs_level = level
            rel_level = self.levels - level
            abs_index = index
            rel_index = (index - self.offset) >> abs_level
            pass
        
        if rel_level == 0:
            if not np.all(rel_index == (0, 0, 0)):
                raise IndexError("Invalid index {} at level {} on OctNode".format(rel_index, rel_level))
            return self
        
        shift = rel_level - 1 
        localIndex = (rel_index & (1 << shift)) >> shift
        lInd = tuple(localIndex)
        chAr = self.children
        if chAr is None:
            if create:
                chAr = self.children = np.empty(self.shape, dtype=object)
            else:
                raise NodeError(self, localIndex)
            pass
        
        child = chAr[lInd]
        if child is None:
            if create:
                halfsize = 2**(self.levels - 1)
                offset = self.offset + halfsize * localIndex
                child = chAr[lInd] = type(self)(self.levels - 1, offset)
            else:
                raise NodeError(self, localIndex)
            pass
        
        childIndex = rel_index ^ (localIndex << shift)
        ans = child.GetNode(rel_level - 1, childIndex, relative=True, create=create)
        assert np.all(ans.offset == abs_index)
        return ans
    
    def DelNode(self, level, index, relative=False):
        if relative:
            rel_level = level
            abs_level = self.levels - level
            rel_index = index
            # abs_index = (index << abs_level) + self.offset
        else:
            abs_level = level
            rel_level = self.levels - level
            # abs_index = index
            rel_index = index >> abs_level
            pass
        
        if rel_level == 0:
            raise ValueError("Node cannot delete itself")
        
        shift = rel_level - 1 
        localIndex = (rel_index & (1 << shift)) >> shift
        lInd = tuple(localIndex)
        chAr = self.children
        if  chAr is None:
            raise NodeError(self, localIndex)
        child = chAr[lInd]
        if child is None:
            raise NodeError(self, localIndex)

        if rel_level == 1:
            # delete my child
            self.children[lInd] = None
            if np.all(np.equal(self.children, None)):
                self.children = None
                pass
            
        else:
            # delegate to appropriate child
            shift = rel_level - 1 
            localIndex = (rel_index & (1 << shift)) >> shift
            lInd = tuple(localIndex)
            chAr = self.children
            if  chAr is None:
                raise NodeError(self, localIndex)
            child = chAr[lInd]
            if child is None:
                raise NodeError(self, localIndex)
            childIndex = rel_index ^ (localIndex << shift)
            child.DelNode(rel_level - 1, childIndex, relative=True)
            pass
        return
    
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
        
        if bot_level < self.levels:
            chAr = self.children
            if chAr is not None:
                for child in chAr.flat:
                    if child is not None:
                        for descendant in child.IterDepthFirst(bot_level,
                                                               min(top_level, self.levels - 1)):
                            yield descendant
                            continue
                        pass
                    continue
                pass
            pass
        if self.levels >= bot_level and self.levels <= top_level:
            yield self
            
    def __eq__(self, other):
        if not super(Node,self).__eq__(other):
            return False
                
        if self.children is None:
            if other.children is None:
                return True
            else:
                return False
            
        if other.children is None:
            return False
        
        return np.all(self.children == other.children)
        
    pass



class Tree(object):
    """An octree.
    """
    NodeClass = Node
    def __init__(self, maxLevels):
        """maxLevels is the maximum depth of the tree.
        """
        self.root = self.NodeClass(maxLevels, (0,0,0))
        self.levels = maxLevels
        return
    
    def GetNode(self, *args, **kwargs):
        return self.root.GetNode(*args, **kwargs)
    
    def IterDepthFirst(self, *args):
        return self.root.IterDepthFirst(*args)
    
    def DelNode(self, *args, **kwargs):
        return self.root.DelNode(*args, **kwargs)
    
    def __eq__(self, other):
        return self.root == other.root
    def __ne__(self, other):
        return not (self == other)

    pass

