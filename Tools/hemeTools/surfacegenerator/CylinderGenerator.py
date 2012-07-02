from TriangulatedCylinderSource import TriangulatedCylinderSource
from vtk import vtkSTLWriter

def CylinderGenerator(radius, length, resolution, outfile):
    """Write an STL file of a triangulated, uncapped, origin-centred cylinder
    of the specified radius, length and resolution to the specified file.
    """
    tcs = TriangulatedCylinderSource()
    w = vtkSTLWriter()
    w.SetInputConnection(tcs.GetOutputPort())

    tcs.SetRadius(radius)
    tcs.SetResolution(resolution)
    tcs.SetHeight(length)
    tcs.CappingOff()
    w.SetFileName(outfile)

    w.Write()

if __name__ == "__main__":
    import sys
    
    try:
        radius = float(sys.argv[1])
        length = float(sys.argv[2])
        resolution = int(sys.argv[3])
        outfile = sys.argv[4]
    except:
        print "Usage:\n%s radius(in mm) length(in mm) resolution outfile" % sys.argv[0]
        raise SystemExit(1)

    CylinderGenerator(radius, length, resolution, outfile)
