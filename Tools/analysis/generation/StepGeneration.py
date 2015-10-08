#!/usr/bin/env python

def WriteStlStep(resolution, nh, nH, nw, nL1, nL2, outfile):
    """Generate STL for a backwards facing step, as in Armaly et
al. J. Fluid Mech. 127 p473 (1983). Briefly, flow is along x-axis,
step in y-direction (see Fig. 1)

    ___________________________________________
   /                                          /
  /w                                         /
 /                                          /__
--------------------------------------------  /
/ h                  / /       H             /
---------------------|/                     /
        L1           |----------------------
  y                              L2
  ^
  |
 /--> x
/z


We're going to tile the system with squares (which are then
subdivided to triangles), arguments are:

- the side length of a square (in mm)
- number of squares along h
- number of squares along H
- number of squares along w
- number of squares along L1
- number of squares along L2
- output filename

"""
    import numpy as np
    from vtk import VTK_TRIANGLE
    from enthought.tvtk.api import tvtk as vtk

    h = resolution * nh
    H = resolution * nH
    w = resolution * nw
    L1 = resolution * nL1
    L2 = resolution * nL2
    
    pd = vtk.PolyData()
    pd.points = vtk.Points()
    
    pointList = vtk.MergePoints()
    pointList.data_set = pd
    pointList.divisions = [(nL1+nL2)/4, nH/4, nw/4]
    pointList.init_point_insertion(pd.points, (0., nL1 + nL2, 0., nH, 0., nw))
    
    
    def appendPoint(pt):
        ptId = pointList.is_inserted_point(pt)
        if ptId >= 0:
            return ptId
        return pointList.insert_next_point(pt)
    
    class Quadder(object):
        def __init__(self, shifts):
            self.shifts = shifts
            self.quad = np.zeros(4, dtype=int)
            return
        
        def __call__(self, pt):
            for i in xrange(4):
                self.quad[i] = appendPoint(pt + self.shifts[i])
                continue
            return self.quad
        pass
    
    quadYp = Quadder([np.array([0., 0., 0.]),
                     np.array([0., 0., 1.]),
                     np.array([1., 0., 1.]),
                     np.array([1., 0., 0.])])
    quadYm = Quadder([np.array([0., 0., 0.]),
                     np.array([1., 0., 0.]),
                     np.array([1., 0., 1.]),
                     np.array([0., 0., 1.])])
    quadX = Quadder([np.array([0., 0., 0.]),
                     np.array([0.,-1., 0.]),
                     np.array([0.,-1., 1.]),
                     np.array([0., 0., 1.])])
    quadZp = Quadder([np.array([0., 0., 0.]),
                     np.array([1., 0., 0.]),
                     np.array([1., 1., 0.]),
                     np.array([0., 1., 0.])])
    quadZm = Quadder([np.array([0., 0., 0.]),
                     np.array([0., 1., 0.]),
                     np.array([1., 1., 0.]),
                     np.array([1., 0., 0.])])
    
    pt = np.zeros(3)

    nTris = 2 * (2*(nL1+nL2)*nw + (nH - nh)*nw + 2*(nL1*nh + nL2*nH))
    pd.allocate(nTris, nTris)
    
    def addQuad(quad):
        pd.insert_next_cell(VTK_TRIANGLE, quad[0:3])
        pd.insert_next_cell(VTK_TRIANGLE, quad[[0,2,3]])
        return
    
    # Top: x = 0..L1+L2; y = H; z = 0..w
    for pt[0] in xrange(0, nL1+nL2):
        pt[1] = nH
        for pt[2] in xrange(0, nw):
            addQuad(quadYp(pt))
            continue
        continue
    
    # Bottom: x = 0..L1; y = H-h; z = 0..w
    for pt[0] in xrange(0, nL1):
        pt[1] = nH-nh
        for pt[2] in xrange(0, nw):
            addQuad(quadYm(pt))
            continue
        continue
    
    # Step: x = L1; y =H-h .. 0; z = 0..w
    pt[0] = nL1
    for pt[1] in xrange(nH-nh, 0, -1):
        for pt[2] in xrange(0, nw):
            addQuad(quadX(pt))
            continue
        continue
    
    # Rest of bottom: x = L1..L1+L2; y = 0; z = 0..w
    for pt[0] in xrange(nL1, nL1+nL2):
        pt[1] = 0
        for pt[2] in xrange(0, nw):
            addQuad(quadYm(pt))
            continue
        continue

    # Side, first part: x = 0..L1; y = H-h..H; z = 0 or w
    for pt[0] in xrange(0, nL1):
        for pt[1] in xrange(nH-nh, nH):
            pt[2] = 0
            addQuad(quadZm(pt))
            pt[2] = nw
            addQuad(quadZp(pt))
            continue
        continue
    # Side, second part: x = L1..L1+L2; y = 0..H; z = 0 or w
    for pt[0] in xrange(nL1, nL1+nL2):
        for pt[1] in xrange(0, nH):
            pt[2] = 0
            addQuad(quadZm(pt))
            pt[2] = nw
            addQuad(quadZp(pt))            
            continue
        continue
    
    # Now scale by the resolution
    pointArray = pd.points.to_array()
    pointArray *= resolution
    
    stl = vtk.STLWriter()
    stl.file_name = outfile
    stl.input = pd
    stl.write()

    return

if __name__ == "__main__":
    import sys
    try:
        args = iter(sys.argv[1:])
        resolution = float(args.next())
        nh = int(args.next())
        nH = int(args.next())
        nw = int(args.next())
        nL1 = int(args.next())
        nL2 = int(args.next())
        outfile = args.next()
        try:
            args.next()
        except StopIteration:
            # OK - args exhausted
            pass
        else:
            # Arguments left!
            raise ValueError
    except:
        print WriteStlStep.__doc__
        raise SystemExit(1)
    else:
        WriteStlStep(resolution, nh, nH, nw, nL1, nL2, outfile)        
    
