#ifndef HEMELB_CONSTANTS_H
#define HEMELB_CONSTANTS_H

const unsigned int COLLISION_TYPES = 6;
const double PI = 3.14159265358979323846264338327950288;
const double DEG_TO_RAD = (PI / 180.0);
// TODO this was used for a convergence test - we could reinstate that at some point.
const double EPSILON = 1.0e-30;

// TODO almost certainly filth.
const double BIG_NUMBER = 1.0e+30;
const int BIG_NUMBER2 = 1 << 30;
const unsigned int BIG_NUMBER3 = 1U << 31U;

const double REFERENCE_PRESSURE_mmHg = 80.0;
const double mmHg_TO_PASCAL = 133.3223874;
const double BLOOD_DENSITY_Kg_per_m3 = 1000.0;
const double BLOOD_VISCOSITY_Pa_s = 0.004;
const double PULSATILE_PERIOD_s = 60.0 / 70.0;
const unsigned int VIS_FIELDS = 3;

// the constants needed to define the configuration of the lattice
// sites follow
const unsigned int BOUNDARIES = 4U;
const unsigned int INLET_BOUNDARY = 0U;
const unsigned int OUTLET_BOUNDARY = 1U;
const unsigned int WALL_BOUNDARY = 2U;
const unsigned int CHARACTERISTIC_BOUNDARY = 3U;

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

#endif //HEMELB_CONSTANTS_H
