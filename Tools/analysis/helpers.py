import math
import numpy as np

def space(grid_coords,voxel_size,origin):
  return grid_coords*voxel_size+origin 

def aligned_cylindrical(grid_coords,voxel_size,origin):
  cartesian=space(grid_coords,voxel_size,origin)
  x,y,z=cartesian.transpose()
  ans=np.empty_like(cartesian)
  ans[:,0]=np.sqrt(x**2+y**2)
  ans[:,1]=np.arctan2(y,x)
  ans[:,2]=z
  return ans

eval_helpers=globals()

# Computes the theoretical stress tensor for 3D Poiseuille flow in the z-direction at a point (radius, 0, z) for an arbitrary z plane.
#
# The tensor has the form:
#     [ -p  0  a ]
#     [  0 -p  0 ], with a = 0.5 * \Delta_p/\Delta_x * radius
#     [  a  0 -p ]
#
# The upper-triangular part of the tensor is stored row-wise in a 1D numpy.array and returned.
def poiseuille_stress_tensor(radius, pressure, pressure_diff_over_pipe_length):
  dv_z_over_dx = 0.5 * pressure_diff_over_pipe_length * abs(radius)
  return np.array([pressure,0.0,dv_z_over_dx,pressure,0.0,pressure])

# Computes the Frobenius norm of a symmetric matrix stored in compress format. symmetric_matrix is the upper-triangular part of
# the matrix stored row-wise in a 1D numpy.array
def frobenius_norm(symmetric_matrix):
  # This mask is used to exclude diagonal elements of the matrix from arithmetic operations
  off_diagonal_mask = mask=np.array([True,False,False,True,False,True]);
  # Square the matrix element-wise
  squared_symmetric_matrix = symmetric_matrix**2
  # Compute the Frobenius norm of symmetric_matrix. Off-diagonal matrix elements are added twice to account for the not stored strictly lower part of the matrix.
  norm = math.sqrt(sum(squared_symmetric_matrix) + np.ma.sum(np.ma.array(squared_symmetric_matrix, mask=off_diagonal_mask)))
  return norm
