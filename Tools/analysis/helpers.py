import math
import numpy as np
import scipy.special

def space(grid_coords,voxel_size,origin):
  return [grid_coord*voxel_size+origin for grid_coord in grid_coords] 

def aligned_cylindrical(grid_coords, voxel_size, origin, axis):
    cartesian=space(grid_coords,voxel_size,origin)

    # To get theta, we have to construct an 'up', i.e. theta = 0 direction.
    # We make this up = (0,1,0) - q * axis, with up.axis = 0.
    # I.e. it's a unit in the y-direction, moved along axis until we're perpendicular
    # to the axis.
    # Some straightforward maths gives q = axis_y, if axis is normalised.
    up = [-axis[1] * axis[0], 1 - axis[1] * axis[1], -axis[1] * axis[2]]
    up = [val / math.sqrt(sum([x*x for x in up])) for val in up]
    result =[]
    for cart in cartesian:
        aligned_z = sum([x*y for x,y in zip(cart,axis)])
        off_axis = [x - aligned_z * y for x,y in zip(cart,axis)]
        r = sum([x*x for x in off_axis]) ** 0.5
        theta = math.acos(
            sum([x*y for x,y in zip(up, off_axis)]) / r)
        result.append([r, theta, aligned_z])
    return result

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

def womersley_velocity(radial_position, time, womersley, pipe_radius, pipe_length, pressure_amplitude, pressure_period, density):
  # Eq. (5.1) in Formaggia et al. "Cardiovascular Mathematics"
  bessel_numerator = scipy.special.jn(0, pow(1j, 3.0/2.0)*womersley*(radial_position/pipe_radius))
  bessel_denominator = scipy.special.jn(0, pow(1j, 3.0/2.0)*womersley)
  sol = (((pressure_amplitude/pipe_length * pressure_period) / (2j * np.pi * density))
         * (1 - bessel_numerator/bessel_denominator)
         * np.exp(2j * np.pi * (time/pressure_period + 0.25))).real
  return sol

def vector_magnitude(vector):
  return np.linalg.norm(vector)

def has_converged(field_1, field_2, max_rel_epsilon, norm=lambda x,y: (x-y)/x):
  diff_field = [abs(norm(x,y)) for (x,y) in zip(field_1, field_2) ]
  index = diff_field.index(max(diff_field))
  return max(diff_field) < max_rel_epsilon

def map_velocity_pair(vel_1, vel_2):
  return sum((pow(x-y,2) for (x,y) in zip(vel_1, vel_2)))**0.5

def max_vector_magnitude(set):
  return max(map(vector_magnitude, set))

def rms(set):
  return (sum(x*x for x in set) / len(set)) ** 0.5

def ave(set):
  return sum(set) / float(len(set))

# Takes two sets of fields (which are themselves sets of points) zips them, maps them reduces them.
def zip_map_reduce(field_set_set_1, field_set_set_2, inner_map, inner_reduction, norm, reduction):
  zipped_together = zip(field_set_set_1, field_set_set_2)

  # Calculate the normalisation
  normalisation = [norm(set_1, set_2) for (set_1, set_2) in zipped_together]

  # Map each point pair to, e.g., the vector magnitude of the difference
  mapped = [ [inner_map(point_1, point_2) for (point_1, point_2) in zip(set_1, set_2)] for (set_1, set_2) in zipped_together ]

  # Reduce the sets to, e.g., their max / rms
  reductions = [inner_reduction(x) for x in mapped]

  # Take each reduced set, multiply by its normalisation then reduce that with the outer reduction function.
  result = reduction( [x*y for (x,y) in zip(reductions, normalisation)])
  return result

def rotate_to_axis(z_vel, axis):
  return [z_vel*x for x in axis]

class TBDAlgorithm(object):
  """
  This class implements the three-band decomposition (TBD) algorithm
  as described in Gizzi et al. 2011. A TBDAlgorithm object can
  be used to compute the TBD of a given wall shear stress signal and
  multiple risk factors (via multiple calls to compute_tbd_for_threshold()).
  """
  # Identifiers of the different bands in the algorithm
  negative_band, zero_band, positive_band, undefined_band = range(4)

  def __init__(self, signed_wss_magnitude):
    """
    Constructor.

    signed_wss_magnitude -- Iterable with the signed magnitude
    of the wall shear stress (wss) vector at a given surface point
    over time. The sign is given by the dot product of a given wss
    vector and the average of all of them.
    """
    self.signed_wss_magnitude = signed_wss_magnitude

  def compute_tbd_for_threshold(self, threshold):
    """
    Compute the TBD for a given threshold. Returns a numpy array with
    the number of negative, zero, and positive intervals resulting from
    the TBD analysis of the signal provided in the constructor and
    threshold.

    threshold -- risk factor of the TBD. In the same units of force
    per area as the signed_wss_magnitude argument of the constructor
    """
    self._init_tbd(threshold)
    
    for wss_signal_value in self.signed_wss_magnitude:
      self._update_tbd(wss_signal_value)

    return self._get_tbd_result()

  def _init_tbd(self, threshold):
    """
    Initialise the TBD for a given risk factor.

    threshold -- risk factor of the TBD. In the same units of force
    per area as the signed_wss_magnitude argument of the constructor
    """    
    self.threshold = threshold
    self.current_state = self.undefined_band
    self.num_intervals = np.array([0, 0, 0])

  def _update_tbd(self, new_wss):
    """
    Updates the TBD by processing a value of wss.

    new_wss -- signed magnitude of the wss vector at a given time
    """
    old_state = self.current_state
    self.current_state = self._get_band_for_value(new_wss)
    if (old_state != self.current_state):
      self.num_intervals[self.current_state] += 1

  def _get_band_for_value(self, wss_value):
    """
    Returns to which band a given value of wss belongs.

    wss_value -- signed magnitude of the wss vector at a given time
    """
    if (wss_value > self.threshold):
      return self.positive_band
    elif (wss_value < -self.threshold):
      return self.negative_band
    else:
      return self.zero_band

  def _get_tbd_result(self):
    """
    Returns a numpy array with the number of negative, zero, and positive
    intervals resulting from the TBD
    """
    return (self.threshold, self.num_intervals[0], self.num_intervals[1], self.num_intervals[2])

def compute_tbd(signed_wss_magnitude):
  """
  Computes the TBD decomposition of a signed wss signal. Returns a list of
  numpy arrays. The cardinality of the list is the number risk factors considered.
  Each numpy array has the number of negative, zero, and positive intervals
  resulting from the TBD.
  """
  tbd_list = []
  tbd_algorithm = TBDAlgorithm(signed_wss_magnitude)

  for risk_factor in np.arange(0.05, max(signed_wss_magnitude), 0.01):
    tbd = tbd_algorithm.compute_tbd_for_threshold(risk_factor)
    tbd_list.append(tbd)

  return tbd_list
