#ifndef HEMELB_CONSTANTS_H
#define HEMELB_CONSTANTS_H

#define COLOURED_PIXELS_MAX    2048 * 2048

#define COLLISION_TYPES                6
#define NEIGHBOUR_PROCS_MAX            64
#define PI           3.14159265358979323846264338327950288
#define DEG_TO_RAD   (PI / 180.0)
#define EPSILON      1.0e-30

#define VON_MISES_STRESS   +1.0
#define SHEAR_STRESS       -1.0

#define STABLE                 1
#define UNSTABLE               0
#define STABLE_AND_CONVERGED   2

#define MACROSCOPIC_PARS   5
#define DENSITY            0
#define VELOCITY           1
#define STRESS             2

#define COMMS_LEVELS                   2

#define BCAST_FREQ   1

#define RECV_BUFFER_A   0
#define RECV_BUFFER_B   1

#define STEERABLE_PARAMETERS   20

#define REFERENCE_PRESSURE             80.0           // 80 mmHg
#define mmHg_TO_PASCAL                 133.3223874
#define BLOOD_DENSITY                  1000.0        // 1000 Kg m^(-3)
#define BLOOD_VISCOSITY                0.004         // 0.004 Pascal s
#define PULSATILE_PERIOD               0.857142857   // period of oscillation (in s) is
					             // chosen to be 1 min / 70
					             // beats per min
#define TOL                            1.0e-6

#define VIS_FIELDS                     3

// the constants needed to define the configuration of the lattice
// sites follow
const unsigned int SOLID_TYPE  = 0U;
const unsigned int FLUID_TYPE  = 1U;
const unsigned int INLET_TYPE  = 2U;
const unsigned int OUTLET_TYPE = 3U;

const unsigned int BOUNDARIES              = 4U;
const unsigned int INLET_BOUNDARY          = 0U;
const unsigned int OUTLET_BOUNDARY         = 1U;
const unsigned int WALL_BOUNDARY           = 2U;
const unsigned int CHARACTERISTIC_BOUNDARY = 3U;

const unsigned int SITE_TYPE_BITS       = 2U;
const unsigned int BOUNDARY_CONFIG_BITS = 14U;
const unsigned int BOUNDARY_DIR_BITS    = 4U;
const unsigned int BOUNDARY_ID_BITS     = 10U;

const unsigned int BOUNDARY_CONFIG_SHIFT = 2U;   // SITE_TYPE_BITS;
const unsigned int BOUNDARY_DIR_SHIFT    = 16U;  // BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
const unsigned int BOUNDARY_ID_SHIFT     = 20U;  // BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

const unsigned int SITE_TYPE_MASK       = ((1U <<  2U) - 1U);         // ((1U << SITE_TYPE_BITS) - 1U);
const unsigned int BOUNDARY_CONFIG_MASK = ((1U << 14U) - 1U) << 2U;   // ((1U << BOUNDARY_CONFIG_BITS) - 1U) << BOUNDARY_CONFIG_SHIFT;
const unsigned int BOUNDARY_DIR_MASK    = ((1U <<  4U) - 1U) << 16U;  //((1U << BOUNDARY_DIR_BITS) - 1U)    << BOUNDARY_DIR_SHIFT;
const unsigned int BOUNDARY_ID_MASK     = ((1U << 10U) - 1U) << 20U;  // ((1U << BOUNDARY_ID_BITS) - 1U)     << BOUNDARY_ID_SHIFT
const unsigned int PRESSURE_EDGE_MASK   = 1U << 31U;

const unsigned int FLUID  = 1U;
const unsigned int INLET  = 2U;
const unsigned int OUTLET = 4U;
const unsigned int EDGE   = 8U;

// square of the speed of sound
const double Cs2 = 1.0 / 3.0;

#endif //HEMELB_CONSTANTS_H
