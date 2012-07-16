from TriangulatedCylinderSource import TriangulatedCylinderSource
from vtk import vtkSTLWriter
import numpy as np

def CylinderGenerator(radius, length, resolution, outfile, direction=(0,0,1)):
    """Write an STL file of a triangulated, uncapped, origin-centred cylinder
    of the specified radius, length and resolution to the specified file.
    """
    tcs = TriangulatedCylinderSource()
    w = vtkSTLWriter()
    w.SetInputConnection(tcs.GetOutputPort())

    tcs.SetRadius(radius)
    tcs.SetResolution(resolution)
    tcs.SetHeight(length)
    tcs.SetDirection(direction)
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
        direction = np.array([0., 0., 1.])
        
        try:
            rotx = np.deg2rad(float(sys.argv[5]))
            matx = np.array([[1.,           0.,           0.],
                             [0., np.cos(rotx),-np.sin(rotx)],
                             [0., np.sin(rotx), np.cos(rotx)]])
            
            roty = np.deg2rad(float(sys.argv[6]))    
            maty = np.array([[ np.cos(roty), 0., np.sin(roty)],
                             [           0., 1.,           0.],
                             [-np.sin(roty), 0., np.cos(roty)]])
            
            direction = maty.dot(matx.dot(direction))
            
        except IndexError:
            pass
        
    except:
        print "Usage:\n%s radius(in mm) length(in mm) resolution outfile [ rotation_about_x_axis rotation_about_y_axis ]" % sys.argv[0]
        raise SystemExit(1)

    CylinderGenerator(radius, length, resolution, outfile, direction)
