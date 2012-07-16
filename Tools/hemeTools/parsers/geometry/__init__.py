# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

"""Regarding indices, a few conventions:

1) Broadly there are two types of index:

 - Three-dimensional indices into a 3D array, these are typically a
   numpy.ndarray (with shape == (3,)) and are named with a suffix of
   'Idx'.

 - One-dimensional indices into a flattened array. These are just
   integers and have the suffix 'Ijk'.

2) Indices can refer to a number of things and have additional naming:

 -  b : Index of a block

 - sg : Index of a site in the whole domain (site global)

 - sl : Index of a site in the block (site local)
 
"""
import numpy as np

GeometryMagicNumber = 0x676d7904
MooreNeighbourhoodDirections = np.array(
    [[-1,-1,-1],
     [-1,-1, 0],
     [-1,-1,+1],
     [-1, 0,-1],
     [-1 , 0, 0],
     [-1, 0,+1],
     [-1,+1,-1],
     [-1,+1, 0],
     [-1,+1,+1],
     [ 0,-1,-1],
     [ 0,-1, 0],
     [ 0,-1,+1],
     [ 0, 0,-1],
     #[ 0, 0, 0], <= the null displacement is not part of the Moore N'hood
     [ 0, 0,+1],
     [ 0,+1,-1],
     [ 0,+1, 0],
     [ 0,+1,+1],
     [+1,-1,-1],
     [+1,-1, 0],
     [+1,-1,+1],
     [+1, 0,-1],
     [+1, 0, 0],
     [+1, 0,+1],
     [+1,+1,-1],
     [+1,+1, 0],
     [+1,+1,+1]]
    )
