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
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('radius', type=float, help='cylinder radius / mm')
    p.add_argument('length', type=float, help='cylinder length / mm')
    p.add_argument('resolution', type=int, help='number of segments around circumference')
    p.add_argument('outfile', help='output file name')
    p.add_argument('rotx', type=float, help='rotation about x-axis in degrees', nargs='?', default=0.)
    p.add_argument('roty', type=float, help='rotation about y-axis in degrees', nargs='?', default=0.)

    args = p.parse_args()
        
    direction = np.array([0., 0., 1.])
        
    rotx = np.deg2rad(args.rotx)
    matx = np.array([[1.,           0.,           0.],
                     [0., np.cos(rotx),-np.sin(rotx)],
                     [0., np.sin(rotx), np.cos(rotx)]])
            
    roty = np.deg2rad(args.roty)
    maty = np.array([[ np.cos(roty), 0., np.sin(roty)],
                     [           0., 1.,           0.],
                     [-np.sin(roty), 0., np.cos(roty)]])
            
    direction = maty.dot(matx.dot(direction))
    
    CylinderGenerator(args.radius, args.length, args.resolution, args.outfile, direction)
