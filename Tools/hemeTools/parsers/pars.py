"""Parse the inlet/outlet parameters file that results from the HemeLB segtool.
"""
import os.path
from ..Plane import Plane, Pressure

def parse(filename):
    pars = file(filename)
    
    numIn = int(pars.readline())
    # get pressure info
    inPress = [Pressure(*lineToVec(pars.readline().strip()))
               for i in range(numIn)]
    
    numOut = int(pars.readline())
    # get pressure info
    outPress = [Pressure(*lineToVec(pars.readline().strip()))
                for i in range(numOut)]
    
    inNorms = [pars.readline().strip() for i in range(numIn)]
    outNorms = [pars.readline().strip() for i in range(numOut)]
    
    inCents = [pars.readline().strip() for i in range(numIn)]
    outCents = [pars.readline().strip() for i in range(numOut)]
    
    inlets = [Plane(normal=lineToVec(inNorms[i]),
                    centre=lineToVec(inCents[i]),
                    name='Inlet %d' % i,
                    pressure=inPress[i]) for i in range(numIn)]
    
    outlets = [Plane(normal=lineToVec(outNorms[i]),
                     centre=lineToVec(outCents[i]),
                     name='Outlet %d' % i,
                     pressure=outPress[i]) for i in range(numOut)]
    
    return inlets + outlets

def lineToVec(line):
    return [float(x) for x in line.split()]
