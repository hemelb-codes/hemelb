#ifndef HEMELB_CONSTANTS_H
#define HEMELB_CONSTANTS_H

#include <limits>

#include "mpiInclude.h"

namespace hemelb
{

  // Basic types for use in HemeLB. Any variable which scales as a function of the number of sites
  // can have type site_t, processors proc_t.
  // Any variable whose precision should roughly match that of the lattice sites' velocity
  // distributions can have type distribn_t.
  typedef long site_t;
  const MPI_Datatype site_mpi_t = MPI_LONG;
  typedef int proc_t;
  typedef double distribn_t;

  const unsigned int COLLISION_TYPES = 6;
  const double PI = 3.14159265358979323846264338327950288;
  const double DEG_TO_RAD = (PI / 180.0);
  // TODO this was used for a convergence test - we could reinstate that at some point.
  const double EPSILON = 1.0e-30;

  // TODO almost certainly filth.
  const distribn_t NO_VALUE = std::numeric_limits<distribn_t>::max();
  const int BIG_NUMBER2 = 1 << 30;
  const unsigned int BIG_NUMBER3 = 1U << 31U;

  const double REFERENCE_PRESSURE_mmHg = 80.0;
  const double mmHg_TO_PASCAL = 133.3223874;
  const double BLOOD_DENSITY_Kg_per_m3 = 1000.0;
  const double BLOOD_VISCOSITY_Pa_s = 0.004;
  const double PULSATILE_PERIOD_s = 60.0 / 70.0;
  const unsigned int VIS_FIELDS = 3;

  // These constants are used to pack a lot of stuff into a 32 bit int.
  // They are used in the setup tool and must be consistent.

  /* This is the number of boundary types. It was 4, but the
   * "CHARACTERISTIC_BOUNDARY" type is never used and I don't know what it is
   * meant to be. It is also not used in the setup tool, so we will drop it,
   * setting BOUNDARIES to 3
   */
  const unsigned int BOUNDARIES = 3U;
  const unsigned int INLET_BOUNDARY = 0U;
  const unsigned int OUTLET_BOUNDARY = 1U;
  const unsigned int WALL_BOUNDARY = 2U;
  // const unsigned int CHARACTERISTIC_BOUNDARY = 3U;

  const unsigned int SITE_TYPE_BITS = 2U;
  const unsigned int BOUNDARY_CONFIG_BITS = 14U;
  const unsigned int BOUNDARY_DIR_BITS = 4U;
  const unsigned int BOUNDARY_ID_BITS = 10U;

  const unsigned int BOUNDARY_CONFIG_SHIFT = 2U; // SITE_TYPE_BITS;
  const unsigned int BOUNDARY_DIR_SHIFT = 16U; // BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
  const unsigned int BOUNDARY_ID_SHIFT = 20U; // BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

  const unsigned int SITE_TYPE_MASK = ( (1U << 2U) - 1U); // ((1U << SITE_TYPE_BITS) - 1U);
  const unsigned int BOUNDARY_CONFIG_MASK = ( (1U << 14U) - 1U) << 2U; // ((1U << BOUNDARY_CONFIG_BITS) - 1U) << BOUNDARY_CONFIG_SHIFT;
  const unsigned int BOUNDARY_DIR_MASK = ( (1U << 4U) - 1U) << 16U; //((1U << BOUNDARY_DIR_BITS) - 1U)    << BOUNDARY_DIR_SHIFT;
  const unsigned int BOUNDARY_ID_MASK = ( (1U << 10U) - 1U) << 20U; // ((1U << BOUNDARY_ID_BITS) - 1U)     << BOUNDARY_ID_SHIFT
  const unsigned int PRESSURE_EDGE_MASK = 1U << 31U;

  const unsigned int FLUID = 1U;
  const unsigned int INLET = 2U;
  const unsigned int OUTLET = 4U;
  const unsigned int EDGE = 8U;

  // square of the speed of sound
  const double Cs2 = 1.0 / 3.0;
}

#endif //HEMELB_CONSTANTS_H
