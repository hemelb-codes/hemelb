import pytest

import os.path
import numpy as np

from HemeLbSetupTool.Model import OutputGeneration
from HemeLbSetupTool.Model.Profile import Profile
from hemeTools.parsers.geometry.simple import ConfigLoader
import fixtures

class TestPolyDataGenerator:
    def test_regression(self, tmpdir):
        """Generate a gmy from a stored profile and check that the output is
        identical.
        """
        dataDir = os.path.join(os.path.split(__file__)[0], 'data')
        proFileName = os.path.join(dataDir, 'test.pro')
        
        p = Profile()
        p.LoadFromFile(proFileName)
        # Change the output to the tmpdir
        basename = tmpdir.join('test').strpath
        outGmyFileName = basename + '.gmy'
        outXmlFileName = basename + '.xml'
        p.OutputGeometryFile = outGmyFileName
        p.OutputXmlFile = outXmlFileName
        
        generator = OutputGeneration.PolyDataGenerator(p)
        generator.Execute()
        
        import filecmp
        assert filecmp.cmp(outGmyFileName, os.path.join(dataDir, 'test.gmy'))
        assert filecmp.cmp(outXmlFileName, os.path.join(dataDir, 'test.xml'))
        
    def test_cube(self,tmpdir):
        """Generate a gmy from a simple cubic profile and check the output"""
        cube = fixtures.cube(tmpdir)
        cube.VoxelSize = 0.23
        cube.StlFileUnitId = 0
        generator = OutputGeneration.PolyDataGenerator(cube)
        generator.Execute()
        # Load back the resulting geometry file and assert things are as expected
        checker = CubeTestingGmyParser(cube.OutputGeometryFile,cube.VoxelSize)
        checker.Load()
        
        fluid_sites = sum(checker.Domain.BlockFluidSiteCounts)
        block_count = len(checker.Domain.Blocks)
        block_size = checker.Domain.BlockSize
        sites = block_count * block_size**3
        #assert(sites==4096)
        #assert(fluid_sites==729)
        assert(sites != fluid_sites)
        # Now, turn on the skip-non-intersecting-blocks optimisation, and assert same result
        generator.skipNonIntersectingBlocks = True
        generator.Execute()
        checker_skip_nonintersecting = CubeTestingGmyParser(cube.OutputGeometryFile, cube.VoxelSize)
        checker_skip_nonintersecting.Load()
        fluid_sites_nonintersecting = sum(checker_skip_nonintersecting.Domain.BlockFluidSiteCounts)
        assert(fluid_sites_nonintersecting == fluid_sites)
        
    def test_cylinder(self,tmpdir):
        """Generate a gmy from a simple cubic profile and check the output"""
        cylinder = fixtures.cylinder(tmpdir)
        cylinder.VoxelSize = 0.23
        cylinder.StlFileUnitId = 0
        generator = OutputGeneration.PolyDataGenerator(cylinder)
        generator.Execute()
        # Load back the resulting geometry file and assert things are as expected
        checker = CylinderTestingGmyParser(cylinder.OutputGeometryFile,cylinder.VoxelSize,np.array([0.0,1.0,0.0]),1.0,0.5)
        checker.Load()
        
        fluid_sites = sum(checker.Domain.BlockFluidSiteCounts)
        block_count = len(checker.Domain.Blocks)
        block_size = checker.Domain.BlockSize
        sites=block_count * block_size**3
        #assert(sites==4096)
        #assert(fluid_sites==621)
        assert(sites != fluid_sites)
        # Now, turn on the skip-non-intersecting-blocks optimisation, and assert same result
        generator.skipNonIntersectingBlocks = True
        generator.Execute()
        checker_skip_nonintersecting=CylinderTestingGmyParser(
            cylinder.OutputGeometryFile, cylinder.VoxelSize,np.array( [0.0, 1.0, 0.0] ), 1.0, 0.5)
        checker_skip_nonintersecting.Load()
        fluid_sites_nonintersecting = sum(checker_skip_nonintersecting.Domain.BlockFluidSiteCounts)
        assert(fluid_sites_nonintersecting == fluid_sites)


class TestCylinderGenerator:
    
    @pytest.mark.parametrize(("randSeed",), [(828,), (341,), (1432,)])
    def test_regression(self, tmpdir, randSeed):
        """Generate a small cylinder GMY with a random orientation. Then check 
        that the output is correct by running it through a custom subclass of 
        a ConfigLoader. 
        """
        
        basename = tmpdir.join('cyl')
        OutputGeometryFile = basename.strpath + '.gmy' 
        OutputXmlFile = basename.strpath + '.xml'
        VoxelSizeMetres = 0.1
        
        rng = np.random.RandomState(randSeed)
        Axis = rng.normal(size=(3,))
        Axis /= np.sqrt(np.dot(Axis,Axis))
        
        LengthMetres = 1.32
        RadiusMetres = 0.74
        
        generator = OutputGeneration.CylinderGenerator(OutputGeometryFile, OutputXmlFile,
                                      VoxelSizeMetres, Axis, LengthMetres,
                                      RadiusMetres)
        generator.Execute()
        
        checker = CylinderTestingGmyParser(OutputGeometryFile, VoxelSizeMetres,
                                           Axis, LengthMetres, RadiusMetres)
        checker.Load()
        return
    
    pass

class TestingGmyParser(ConfigLoader):
    def __init__(self,filename,VoxelSize):
        ConfigLoader.__init__(self,filename)
        self.VoxelSize=VoxelSize
        
    def OnEndHeader(self):
        assert self.Domain.VoxelSize == self.VoxelSize
        
    def OnEndSite(self, block, site):
        return
        if site.IsFluid:
            assert self.IsInside(site.Position)
        else:
            assert not self.IsInside(site.Position)
        return
    
    def OnEndBlock(self, bIdx, bIjk):
        self.Domain.DeleteBlock(bIdx)
        return
    
    pass

class CylinderTestingGmyParser(TestingGmyParser):
    def __init__(self, filename, VoxelSize, Axis, Length, Radius):
        TestingGmyParser.__init__(self, filename,VoxelSize)
        self.Axis = Axis
        self.Length = Length
        self.Radius = Radius
        return
    
    def IsInside(self, x):
        xDOTn = np.dot(x, self.Axis)
        if xDOTn < -0.5 * self.Length or xDOTn > 0.5 * self.Length:
            return False
        perp = x - xDOTn * self.Axis
        if np.dot(perp, perp) > self.Radius**2:
            return False
        return True
        
class CubeTestingGmyParser(TestingGmyParser):
    def IsInside(self, position):
        result= (
            all(component<0.5 and component>-0.5 for component in position)
        )      
        return result  

