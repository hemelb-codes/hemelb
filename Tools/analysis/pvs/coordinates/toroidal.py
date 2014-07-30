"""Define classes to convert between toroidal and Cartesian coordinate
systems.

"""
import numpy as np
from . import UnpackArrayOfTriples, TransformerBase

class ToroidalBase(TransformerBase):
    def __init__(self, a, c):
        TransformerBase.__init__(self)
        
        self.a = a
        self.c = c
        return
    
    def _GetInvInitArgs(self):
        return (self.a, self.c)

class ToroidalToCartesianTransformer(ToroidalBase):
    """A curvilinear coordinate system (a r, theta, phi) is chosen,
    where (a r, theta) are polar coordinates in the plane of cross
    section of the pipe, a is the radius of the cross-section and phi
    is the angle from the positive x-axis. Hence:

    x = (c + a r cos theta)cos phi
    y = -(c + a r cos theta)sin phi
    z = a r sin theta

    where c = RingRadius.
    """
    
    def TransformPoint(self, tCoords):
        r, theta, phi = UnpackArrayOfTriples(tCoords)
        cPlusarCosTheta = self.c + self.a * r * np.cos(theta)
        
        cCoords = np.empty(r.shape + (3,), dtype=r.dtype) * self.c.units
        cCoords[..., 0] = cPlusarCosTheta * np.cos(phi)
        cCoords[..., 1] = -cPlusarCosTheta * np.sin(phi)
        cCoords[..., 2] = self.a * r * np.sin(theta)
        return cCoords
    
    
    def GetTransformationMatrix(self, tCoords):
        r, theta, phi = UnpackArrayOfTriples(tCoords)
        sinTheta = np.sin(theta)
        cosTheta = np.cos(theta)
        sinPhi = np.sin(phi)
        cosPhi = np.cos(phi)

        L = np.empty(r.shape + (3,3), dtype=r.dtype)
        
        L[..., 0,0] = cosTheta * cosPhi
        L[..., 1,0] = -cosTheta * sinPhi
        L[..., 2,0] = sinTheta

        L[..., 0,1] = -sinTheta * cosPhi
        L[..., 1,1] = sinTheta * sinPhi
        L[..., 2,1] = cosTheta

        L[..., 0,2] = -sinPhi
        L[..., 1,2] = -cosPhi
        L[..., 2,2] = 0.
        
        return L
    pass
    
class CartesianToToroidalTransformer(ToroidalBase):
    
    def TransformPoint(self, cCoords):
        x, y, z = UnpackArrayOfTriples(cCoords)
        
        tCoords = np.empty(x.shape + (3,), dtype=x.dtype)
        rhoMinusC = np.sqrt(x**2 + y**2) - self.c
        tCoords[..., 0] = np.sqrt(z**2 + (rhoMinusC)**2) / self.a
        # Note z == rSinTheta
        tCoords[..., 1] = np.arctan2(z, rhoMinusC)        
        tCoords[..., 2] = np.arctan2(-y, x)
        return tCoords
    
    def GetTransformationMatrix(self, cCoords):
        # Get the toroidal coords
        tCoords = self.TransformPoint(cCoords)
        # Get the inverse's TM
        Linv = self.GetInverse().GetTransformationMatrix(tCoords)
        # Now transpose it along the LAST TWO dimensions.
        return np.einsum('...ji', Linv)
    
    pass

# Set up the inverses
CartesianToToroidalTransformer._Inverse = ToroidalToCartesianTransformer
ToroidalToCartesianTransformer._Inverse = CartesianToToroidalTransformer
