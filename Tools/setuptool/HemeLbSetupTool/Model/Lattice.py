# This file is part of HemeLB and is CONFIDENTIAL. You may not work
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
#

import numpy as np

Neighbours = np.array((
                      (-1,-1,-1),
                      (-1,-1, 0),
                      (-1,-1,+1),
                      (-1, 0,-1),
                      (-1, 0, 0),
                      (-1, 0,+1),
                      (-1,+1,-1),
                      (-1,+1, 0),
                      (-1,+1,+1),
                      ( 0,-1,-1),
                      ( 0,-1, 0),
                      ( 0,-1,+1),
                      ( 0, 0,-1),
                      #( 0, 0, 0),
                      ( 0, 0,+1),
                      ( 0,+1,-1),
                      ( 0,+1, 0),
                      ( 0,+1,+1),
                      (+1,-1,-1),
                      (+1,-1, 0),
                      (+1,-1,+1),
                      (+1, 0,-1),
                      (+1, 0, 0),
                      (+1, 0,+1),
                      (+1,+1,-1),
                      (+1,+1, 0),
                      (+1,+1,+1)
                      ),
                     dtype=int)

Opposites = np.empty(26, dtype=int)
Opposites[:] = -1
for i in xrange(26):
    for j in xrange(i, 26):
        if np.all(Neighbours[i] == -Neighbours[j]):
            Opposites[i] = j
            Opposites[j] = i
            
for i in xrange(26):
    assert np.all(Neighbours[i] == -Neighbours[Opposites[i]])
