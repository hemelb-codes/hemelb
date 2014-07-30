import abc
import numpy as np
import quantities as pq

def UnpackArrayOfTriples(array):
    if not isinstance(array, np.ndarray):
        array = np.array(array)
        
    assert array.shape[-1] == 3
    return array[...,0], array[..., 1], array[..., 2]

class TransformerBase(object):
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def TransformPoint(self, pt):
        pass
    
    @abc.abstractmethod
    def GetTransformationMatrix(self, pt):
        pass

    @abc.abstractmethod
    def _GetInvInitArgs(self):
        pass
    
    def GetInverse(self):
        return self._Inverse(*self._GetInvInitArgs())
    
    def TransformVector(self, coords, vec):
        L_ij = self.GetTransformationMatrix(coords)
        ans = np.einsum('...ij,...j', L_ij, vec)
        
        if isinstance(vec, pq.Quantity):
            return ans * vec.units
        
        return ans
    pass
