import math
import numpy as np
import scipy.special

def space(grid_coords,voxel_size,origin):
  return [grid_coord*voxel_size+origin for grid_coord in grid_coords] 

def aligned_cylindrical(grid_coords,voxel_size,origin):
  cartesian=space(grid_coords,voxel_size,origin)
  return [[pow(pow(x,2)+pow(y,2),0.5),math.atan2(y,x),z] for x,y,z in cartesian]

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

def womersley_velocity(radial_position, time, womersley, pipe_radius, pressure_amplitude, pressure_period, density):
  print radial_position, time, womersley, pipe_radius, pressure_amplitude, pressure_period, density
  # Eq. (5.1) in Formaggia et al. "Cardiovascular Mathematics"
  bessel_numerator = scipy.special.jn(0, pow(1j, 3.0/2.0)*womersley*(radial_position/pipe_radius))
  bessel_denominator = scipy.special.jn(0, pow(1j, 3.0/2.0)*womersley)
  sol = (((pressure_amplitude * pressure_period) / (2j * np.pi * density))
         * (1 - bessel_numerator/bessel_denominator)
         * np.exp(2j * np.pi * (time/pressure_period))).real
  return sol
