import numpy as np
from scipy.sparse import lil_matrix
from vtk.util import numpy_support as convert
import fipy

import pdb

class PoiseuilleSolver(object):
    def __init__(self, polydata):
        self.mesh = FiPyTriangleMesher(polydata).GetMesh()
        
        self.speed = fipy.CellVariable(name="speed", mesh=self.mesh, value=0.)
        self.BCs = (fipy.FixedValue(faces=self.mesh.getExteriorFaces(), value=0),)
        
        self.equation = fipy.DiffusionTerm() + 1. == 0.
        return
    
    def Solve(self):
        self.equation.solve(var=self.speed, boundaryConditions=self.BCs)
        return self.speed
    
    pass

class FiPyTriangleMesher(object):
    def __init__(self, polyData):
        vertices3D = convert.vtk_to_numpy(polyData.GetPoints().GetData())
        # FiPy needs its coordinates in an array like:
        # [[x0, x1, ..., xn-1],
        #  [y0, y1, ..., yn-1]]
        # VTK gives us:
        # [[x0, y0, z0],
        #  [x1, y1, z1],
        #  ...,
        #  [xn-1, yn-1, zn-1]]
        self.vertices = vertices3D[:,0:2].transpose()

        # VTK gives us a flat array of cells like
        # [nPoints, pId0, pId1, pId2, ... nPoints,... ]
        cells = convert.vtk_to_numpy(polyData.GetPolys().GetData())
        self.nTris = nTris = polyData.GetPolys().GetNumberOfCells()
        cells.shape = (nTris, 4)
        
        # FiPy requires a list of the FACES (in 2D the lines) making up each 
        # cell. Make the array the max possible size for now.
        self.faces = np.zeros((2, 3*nTris), dtype=np.int)
        
        # Since these have to be unique, construct a sparse lookup matrix
         
        # IMPORTANT- since the default value of the matrix is zero, we're 
        # going to use 1-indexing for these so +1 when writing and -1 when 
        # reading.
        self.faceCache = lil_matrix((nTris,nTris), dtype=np.int)
        # Track the number of faces added.
        self.nFaces = 0
        
        # A cell is defined by the IDs of its faces
        self.cells = np.zeros((3, nTris), dtype=np.int)
        # Track the number added
        self.nCells = 0
        
        for tri in cells[:,1:]:
            self.AddTriangle(tri)
        
        self.faces = self.faces[:,:self.nFaces]
        # self.cells = self.cells.flatten()
        return
    
    def GetMesh(self):
        return fipy.meshes.numMesh.mesh2D.Mesh2D(self.vertices, self.faces, self.cells)
    
    def AddTriangle(self, triPtIds):
        for i in xrange(3):
            self.cells[i, self.nCells] = self.InsertUniqueFace(triPtIds[i],
                                                          triPtIds[(i+1) % 3])
        self.nCells += 1
        return
    
    def InsertUniqueFace(self, a, b):
        # Put them in ascending order
        if a > b:
            b, a = a, b
        
        faceId = self.faceCache[a,b] - 1 
        if faceId < 0:
            # It has not been inserted, so insert a face
            faceId = self.nFaces 
            self.nFaces += 1
            self.faceCache[a,b] = faceId + 1
            self.faces[0, faceId] = a
            self.faces[1, faceId] = b
            
        return faceId
    