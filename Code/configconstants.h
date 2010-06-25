#ifndef __configconstants_h_
#define __configconstants_h_

#define COLLISION_TYPES                6
#define NEIGHBOUR_PROCS_MAX            64

// the constants needed to define the configuration of the lattice
// sites follow

const unsigned int SOLID_TYPE  = 0U;
const unsigned int FLUID_TYPE  = 1U;
const unsigned int INLET_TYPE  = 2U;
const unsigned int OUTLET_TYPE = 3U;
const unsigned int NULL_TYPE   = 4U;

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

#endif //__configconstants_h_
