#!/usr/bin/env python
import sys
try:
    r = float(sys.argv[1])
    length = float(sys.argv[2])
    n = int(sys.argv[3])
    outfile = sys.argv[4]
except:
    print "Usage:\n%s radius(in mm) length(in mm) resolution outfile" % sys.argv[0]
    raise SystemExit

from TriangulatedCylinderSource import TriangulatedCylinderSource
from vtk import vtkSTLWriter

tcs = TriangulatedCylinderSource()
w = vtkSTLWriter()
w.SetInputConnection(tcs.GetOutputPort())

tcs.SetResolution(n)
tcs.SetHeight(length)
tcs.CappingOff()
w.SetFileName(outfile)

w.Write()

