import numpy as N
from configobj import ConfigObj

def _floatarray(lst):
    return N.array([float(x) for x in lst])

def _intarray(lst):
    return N.array([float(x) for x in lst])

class Transformer(object):
    """Parses the coordinates output of the HemeLB segtool and allows
    simple conversions between coordinate systems.
    
    """
    
    mapper = {'voxel_size' : float,
              'res_factor' : int,
              'half_dim'   : _floatarray,
              'mesh_centre': _floatarray,
              'seed_pos'   : _floatarray,
              'seed_site'  : _intarray}
              
    def __init__(self, coordsFile):
        """trans = Transformer(coordinate filename)
        
        """
        c = ConfigObj(coordsFile)
        for k in self.mapper:
            setattr(self, k, self.mapper[k](c[k]))
            continue
        
        self.shift = N.zeros(3)
        self.shift = self.siteToMesh(self.seed_site) - self.seed_pos
        
        return
    
    def meshToSite(self, mesh):
        return ((mesh + self.half_dim) * self.res_factor / self.voxel_size).astype(int)
    
    def siteToMesh(self, site):
        return N.array(site, dtype=int) * (self.voxel_size/self.res_factor) - self.half_dim - self.shift
    
    def stlToMesh(self, stl):
        return stl - self.mesh_centre
    
    def meshToStl(self, mesh):
        return mesh + self.mesh_centre
    
    def stlToSite(self, stl):
        return self.meshToSite(self.stlToMesh(stl))
    
    def siteToStl(self, site):
        return self.meshToStl(self.siteToMesh(site))
    pass
