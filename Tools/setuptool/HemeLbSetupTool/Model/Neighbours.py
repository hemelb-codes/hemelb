"""HemeLB's ordering of the vectors linking a site to its 26 neighbours. 
"""
import numpy as np

neighbours = np.array(((-1,-1,-1),
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
                       # ( 0, 0, 0),
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
                       (+1,+1,+1)), dtype=int)
inverses = np.array((25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0))
assert np.all(neighbours + neighbours[inverses] == 0)

norm2 = np.sum(neighbours**2, axis=1)
norm = np.sqrt(norm2)