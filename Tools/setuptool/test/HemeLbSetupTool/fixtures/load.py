import os
from HemeLbSetupTool.Model.Profile import Profile

def fixture_path(name):
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here,name+".stl")

def cube():
    # Create a profile which is a simple cube
    result=Profile()
    result.StlFile=fixture_path('cube')