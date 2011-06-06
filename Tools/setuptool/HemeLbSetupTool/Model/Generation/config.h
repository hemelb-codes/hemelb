#ifndef HEMELBSETUPTOOL_CONFIG_H
#define HEMELBSETUPTOOL_CONFIG_H

// This is from the HemeLB tree
#include "constants.h"

namespace hemelb {
	enum SiteType {
		// This must be consistent with hemelb::geometry::LatticeData::SiteType
		SOLID_TYPE = 0U,
		FLUID_TYPE = 1U,
		INLET_TYPE = 2U,
		OUTLET_TYPE = 3U
	};
}
#endif // HEMELBSETUPTOOL_CONFIG_H
