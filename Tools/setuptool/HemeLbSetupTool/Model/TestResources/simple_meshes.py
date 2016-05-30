# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 
import numpy as np
def mk_trivial():
    # Put the tri points on a square, in the (0,0,0) octant
    points = np.array([(1.2, 1.2, 1.2),
                       (1.2, 1.2, 2.2),
                       (1.2, 2.2, 1.2),
                       (1.2, 2.2, 2.2)], dtype=float)
    # Define the tris
    triangles = np.array([(0,1,2),
                          (2,1,3)], dtype=int)
    normals = np.array([[-1, 0, 0],
                        [-1, 0, 0]], dtype=float)
    return points, triangles, normals

def mk_trivial2():

    points = np.array([(1.2, 1.2, 1.2),
                       (1.2, 1.2, 2.2),
                       (1.2, 2.2, 1.2),
                       (1.2, 2.2, 2.2),
                       (1.2, 1.2, 3.2),
                       (1.2, 2.2, 3.2)], dtype=float)

    triangles = np.array([(0,1,2),
                          (2,1,3),(0,4,2)], dtype=int)
    normals = np.array([[-1, 0, 0],
                        [-1, 0, 0],[-1, 0, 0]], dtype=float)

    return points, triangles, normals

def mk_duct():
    from ..Iolets import Inlet, Outlet, Vector
    # This is a simple square duct with an inlet and outlet
    # Side length is 3 voxels, length is 12 voxels
    # Aligned along the line 2,2, 2.5 - 2,2,14.5
    points = np.array([(0.5, 0.5, 2.5),
                       (0.5, 3.5, 2.5),
                       (3.5, 3.5, 2.5),
                       (3.5, 0.5, 2.5),
                       
                       (0.5, 0.5, 14.5),
                       (0.5, 3.5, 14.5),
                       (3.5, 3.5, 14.5),
                       (3.5, 0.5, 14.5)],
                       
                       dtype=float)
    
    triangles = np.array([# -x
                          (0,4,1),
                          (4,5,1),
                          #+y
                          (1,5,2),
                          (5,6,2),
                          # +x
                          (2,6,3),
                          (6,7,3),
                          # -y
                          (3,7,0),
                          (7,4,0),
                          # inlet
                          (0,1,3),
                          (1,2,3),
                          # outlet
                          (4,5,7),
                          (5,6,7)
                          ],
                         dtype=int)
    normals = np.array([(-1,0,0),
                        (-1,0,0),
                        (0,+1,0),
                        (0,+1,0),
                        (+1,0,0),
                        (+1,0,0),
                        (0,-1,0),
                        (0,-1,0),
                        (0,0,-1),
                        (0,0,-1),
                        (0,0,+1),
                        (0,0,+1),],
                       dtype=float)
    labels = np.array([-1,-1, -1,-1, -1,-1, -1,-1, 0,0, 1,1])
    
    inlet = Inlet(Centre=Vector(2,2,2.5),
                  Normal=Vector(0,0,1),
                  Radius=2.5)
    outlet = Outlet(Centre=Vector(2,2,14.5),
                  Normal=Vector(0,0,-1),
                  Radius=2.5)
    return points, triangles, normals, labels, [inlet, outlet]
