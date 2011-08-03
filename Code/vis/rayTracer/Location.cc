#include "util/utilityFunctions.h"
#include "vis/rayTracer/Location.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      Location::Location() { }

      Location::Location(site_t iI, site_t iJ, site_t iK)
      {
	i = iI;
	j = iJ;
	k = iK;
      }

      Location::Location(site_t iX)
      {
	i=iX;
	j=iX;
	k=iX;
      }

      Location Location::operator+(const Location right)
      {
	return Location(i + right.i,
			j + right.j,
			k + right.k);
      }

      Location Location::operator*(const site_t multiplier)
      {
	return Location(i * multiplier,
			j * multiplier,
			k * multiplier);
      }

      void Location::UpdateMinLocation(Location& io_store_location, Location& i_compare_location)
      {
	io_store_location.i = 
	  util::NumericalFunctions::min
	  (io_store_location.i, 
	   i_compare_location.i);

	io_store_location.j = 
	  util::NumericalFunctions::min
	  (io_store_location.j, 
	   i_compare_location.j);

	io_store_location.k = 
	  util::NumericalFunctions::min
	  (io_store_location.k, 
	   i_compare_location.k);
      }

      void Location::UpdateMaxLocation(Location& io_store_location, const Location& i_compare_location)
      {
	io_store_location.i = 
	  util::NumericalFunctions::max
	  (io_store_location.i, 
	   i_compare_location.i);

	io_store_location.j = 
	  util::NumericalFunctions::max
	  (io_store_location.j, 
	   i_compare_location.j);

	io_store_location.k = 
	  util::NumericalFunctions::max
	  (io_store_location.k, 
	   i_compare_location.k);
      }
    
    }
  }
}
