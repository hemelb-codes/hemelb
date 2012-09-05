import os.path
import numpy as np
from HemeLbSetupTool.Model.OutputGeneration import CylinderGenerator
from hemeTools.parsers.geometry.simple import ConfigLoader

class TestCylinderGenerator:
    def test_regression(self, tmpdir):
        rng = np.random.RandomState()
        
        basename = tmpdir.join('cyl')
        OutputGeometryFile = basename.strpath + '.gmy' 
        OutputXmlFile = basename.strpath + '.xml'
        VoxelSizeMetres = 0.1
        Axis = rng.normal(size=(3,))
        Axis /= np.sqrt(np.dot(Axis,Axis))
        
        LengthMetres = 1.32
        RadiusMetres = 0.74
        
        generator = CylinderGenerator(OutputGeometryFile, OutputXmlFile,
                                      VoxelSizeMetres, Axis, LengthMetres,
                                      RadiusMetres)
        generator.Execute()
        
        checker = CylinderTestingGmyParser(OutputGeometryFile, VoxelSizeMetres,
                                           Axis, LengthMetres, RadiusMetres)
        checker.Load()
        return
    
    pass

class CylinderTestingGmyParser(ConfigLoader):
    def __init__(self, filename, VoxelSize, Axis, Length, Radius):
        ConfigLoader.__init__(self, filename)
        self.VoxelSize = VoxelSize
        self.Axis = Axis
        self.Length = Length
        self.Radius = Radius
        return
    
    def IsInsideCylinder(self, x):
        xDOTn = np.dot(x, self.Axis)
        if xDOTn < -0.5 * self.Length or xDOTn > 0.5 * self.Length:
            return False
        perp = x - xDOTn * self.Axis
        if np.dot(perp, perp) > self.Radius**2:
            return False
        
        return True
        
    def OnEndHeader(self):
        assert self.Domain.VoxelSize == self.VoxelSize
        
    def OnEndSite(self, block, site):
        if site.IsFluid:
            assert self.IsInsideCylinder(site.Position)
        else:
            assert not self.IsInsideCylinder(site.Position)
        return
    
    def OnEndBlock(self, bIdx, bIjk):
        self.Domain.DeleteBlock(bIdx)
        return
    
    pass
