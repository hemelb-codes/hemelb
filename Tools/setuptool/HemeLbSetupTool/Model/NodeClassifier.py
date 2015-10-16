import numpy as np
from HemeOct import SolidNode, FluidNode
import pdb

class NodeClassifier(object):
    
    def __call__(self, node, min_level=0):
        if node.levels < min_level:
            return
        
        if node.levels == 0:
            self._ClassifyLeafNode(node)
        else:
            if node.children is None:
                # We know nothing about our children, so all we know is that 
                # this block is all the same.
                return
            self._ClassifyChildren(node, min_level)
            self._ClassifyComplexNode(node)
            pass
        return
    
    def _ClassifyLeafNode(self, node):
        # We're at the bottom
        # Decide on fluidness
        inVotes = 0
        outVotes = 0
        for hits in node.intersections.itervalues():
            # A hit is a tuple (t, outward_flag, tri_ID)
            # Sort them based on the t, the distance along the link
            hits.sort(key=lambda x: x[0])
            if hits[0][1]:
                # The normal points out along the link, so we're inside this triangle
                inVotes += 1
            else:
                # We're outside this triangle
                outVotes += 1
                pass
                    
            if inVotes and outVotes:
                # Let's look at any non-unanimous links for now
                pdb.set_trace()
                pass
            
            # Count the votes and set the fluidness flag
            if inVotes > outVotes:
                node.fluid_count = 1
            else:
                node.fluid_count = 0
                pass
            # Set the edge flag
            node.edge = True
            continue
        return
    
    def _ClassifyChildren(self, node, min_level):
        for cn in node.children.flat:
            if cn is not None:
                self(cn, min_level)
                
    def _ClassifyComplexNode(self, node):
        know_child_state = np.zeros(node.shape, dtype=bool)
        for ijk, child in np.ndenumerate(node.children):
            if child is not None and hasattr(child, "fluid_count"):
                know_child_state[ijk] = True
                
        if not np.all(know_child_state):
            # Got to infer state of some children
            
            # Get the relative index of one we know.
            inds = know_child_state.nonzero()
            ijk = (inds[0][0], inds[1][0], inds[2][0])
            
#             found = False
#             for ijk, child in np.ndenumerate(node.children):
#                 if child is not None:
#                     if hasattr(child, "fluid_count"):
#                         found = True
#                         break
#                 continue
#             assert found, "Know nothing about node's children!"
            # Get the index, within the child node, of the cell that
            # is closest to the centre of the current node.
            corner_ijk = tuple(x^1 for x in ijk)
            # Now descend down to get the leaf node there,
            known_node = node.children[ijk]
            
            while known_node.levels > 0 and known_node.children is not None:
                known_node = known_node.children[corner_ijk]
                continue
                        
            # Because all unknown nodes have no boundary separating 
            # them from this cell, they are the same.
            if known_node.fluid_count > 0:
                UnknownNodeClass = FluidNode
            else:
                UnknownNodeClass = SolidNode
                pass
            
            for ijk in np.ndindex(*node.shape):
                if not know_child_state[ijk]:
                    halfsize = 2**(node.levels - 1)
                    offset = node.offset + halfsize * np.array(ijk)
                    node.children[ijk] = UnknownNodeClass(node.levels - 1, offset)
                    pass
                continue
            pass
        
        assert np.all(np.not_equal(node.children, None)), "Still missing some children!"
        
        # Now tot up the number of fluid sites contained in this node
        node.fluid_count = 0
        for ch in node.children.flat:
            node.fluid_count += ch.fluid_count
            continue
        
        # We're all solid so can delete the kids.
        if node.fluid_count == 0:
            node.children = None
            pass
        return