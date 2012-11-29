import os
from HemeLbSetupTool.Model.Profile import Profile

def fixture_path(name):
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here,name+".stl")

def temporary_profile(tmpdir,name):
    result=Profile()
    basename = tmpdir.join(name).strpath
    outGmyFileName = basename + '.gmy'
    outXmlFileName = basename + '.xml'
    result.OutputGeometryFile = outGmyFileName
    result.OutputXmlFile = outXmlFileName
    return result

def cube(tmpdir):
    # Create a profile which is a simple cube
    result=temporary_profile(tmpdir,'cube')
    result.StlFile=fixture_path('cube')
    return result