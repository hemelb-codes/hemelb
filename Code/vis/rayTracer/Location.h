#ifndef HEMELB_VIS_RAYTRACER_LOCATION_H
#define HEMELB_VIS_RAYTRACER_LOCATION_H


#include "constants.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class Location
      {
      public:
	site_t i, j, k;

	Location();

	Location(site_t iI, site_t iJ, site_t iK);

	Location(site_t iX);
	
	Location operator+(const Location right);

	Location operator*(const site_t multiplier);

	static void UpdateMinLocation(Location& io_store_location, Location& i_compare_location);
      
	static void UpdateMaxLocation(Location& io_store_location, const Location& i_compare_location);
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_LOCATION_H
