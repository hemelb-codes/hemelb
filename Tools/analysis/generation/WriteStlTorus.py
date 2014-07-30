
import vtk
import numpy as np

def WriteStlTorus(c, a, n_c, n_a, filename):
    """VTK defines a torus by:

x(u,v) = (c+acosv)cosu
y(u,v) = (c+acosv)sinu
z(u,v) = a sinv

        ___-----___
     .--           --.
    /    _-------_    \
   /    /         \    \
  |    |           |    |
  <-2a->     +---c--->  |
  |    |           |    |
   \    \         /    /
    \    `-------'    /
     `--___     ___--'
           -----
  y
  ^
  |
 /--> x
/z

The radius from the centre to the middle of the ring of the torus
is c and a is the radius of the cross-section of the ring of the
torus where c, a > 0. Generally c > a, giving the usual torus
shape.

Then the torus is defined by: S(u,v) = (x(u,v), y(u,v), z(u,v)),
where 0 <= u <= 2 pi, 0 <= v <= 2 pi

Arguments are (in order):
c - radius from the centre to the middle of the ring
a - radius of the cross-section of the ring

n_c - number of sections around the torus
n_a - number of sections around the cross section

filename - name of file to write
"""
    torus = vtk.vtkParametricTorus()
    torus.SetRingRadius(c)                 # i.e. c
    torus.SetCrossSectionRadius(a)         # i.e. a

    torusSource = vtk.vtkParametricFunctionSource()
    torusSource.SetParametricFunction(torus)
    torusSource.SetUResolution(n_c)
    torusSource.SetVResolution(n_a)

    # This cleaner removes the degenerate points at the joins
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetTolerance(1e-6)
    cleaner.SetInputConnection(torusSource.GetOutputPort())

    writer = vtk.vtkSTLWriter()
    writer.SetInputConnection(cleaner.GetOutputPort())
    writer.SetFileName(filename)
    writer.Write()

if __name__ == "__main__":
    import sys
    try:
        args = iter(sys.argv[1:])
        c = float(args.next())
        a = float(args.next())
        n_c = int(args.next())
        n_a = int(args.next())
        filename = args.next()
        try:
            args.next()
        except StopIteration:
            # OK - args exhausted
            pass
        else:
            # Arguments left! Raise any old exception to be caught below
            raise ValueError
    except:
        print "Argument Error"
        print WriteStlTorus.__doc__
        raise SystemExit(1)
    else:
        WriteStlTorus(c, a, n_c, n_a, filename)
