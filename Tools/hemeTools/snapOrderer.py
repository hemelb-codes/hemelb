"""Since the line-by-line ordering of snapshots can vary depending on
decomposition, this enforces a canonical ordering.

The ordering is based on the grid positions, z fastest varying.

"""

from .snapshots import HemeLbSnapshot
import numpy as N

def computeGridIdx(snap):
    return (snap.grid[:,0] * snap.bb_len[1] +
            snap.grid[:,1]) * snap.bb_len[2] + \
            snap.grid[:,2]

def order(snap):
    idx = computeGridIdx(snap)
    idxIdx = idx.argsort()
    
    new = snap.copy()

    for i, ind in enumerate(idxIdx):
        new[i] = snap[ind]
        continue
    
    return new
