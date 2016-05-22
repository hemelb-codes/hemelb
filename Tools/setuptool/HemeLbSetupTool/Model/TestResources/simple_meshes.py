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
