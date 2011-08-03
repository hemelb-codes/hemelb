#ifndef HEMELB_VIS_RAYTRACER_LOCATION_H
#define HEMELB_VIS_RAYTRACER_LOCATION_H

#include "constants.h"
#include "util/utilityFunctions.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      template <class T = site_t> 
	class Location
	{
	public:
	T i, j, k;

	Location() {};

	Location(T iI, T iJ, T iK) :
	i(iI), j(iJ), k(iK)
	{}
      
	Location(T iX) :
	i(iX), j(iX), k(iX)
	{}

	template < class OldTypeT >
	Location<T>(const Location<OldTypeT> & iOldLocation)
	{
	  i = static_cast<T>(iOldLocation.i);
	  j = static_cast<T>(iOldLocation.j);
	  k = static_cast<T>(iOldLocation.k);
	}
	
	Location<T> operator+(const Location<T> right)
	{
	  return Location(i + right.i,
			  j + right.j,
			  k + right.k);
	}

	Location<T> operator-(const Location<T> right)
	{
	  return Location(i - right.i,
			  j - right.j,
			  k - right.k);
	}

	template < class MultiplierT >
	Location<T> operator*(const MultiplierT multiplier)
	{
	  return Location(i * multiplier,
			  j * multiplier,
			  k * multiplier);
	}

	static Location<T> MaxLimit() 
	{
	  return Location(std::numeric_limits<T>::max());
	}

	static Location<T> MinLimit()
	{
	  return Location(std::numeric_limits<T>::min());
	}

	static void UpdateMinLocation(Location& io_store_location, const Location& i_compare_location)
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
      
	static void UpdateMaxLocation(Location& io_store_location, const Location& i_compare_location)
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
	
	};
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_LOCATION_H
