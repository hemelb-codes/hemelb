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
      //Location is essentially a 3D vector, storing the 
      //x, y and z co-ordinate in the templated numeric type
      //Other methods are defined for convenience
      template <class T = site_t> 
	class Location
	{
	public:
	T x, y, z;

	Location() {};

	Location(T iX, T iY, T iZ) :
	x(iX), y(iY), z(iZ)
	{}
      
	Location(T iX) :
	x(iX), y(iX), z(iX)
	{}

	//Copy constructor - can be used to perform type converstion 
	template < class OldTypeT >
	Location<T>(const Location<OldTypeT> & iOldLocation)
	{
	  x = static_cast<T>(iOldLocation.x);
	  y = static_cast<T>(iOldLocation.y);
	  z = static_cast<T>(iOldLocation.z);
	}
	
	//Vector addition
	Location<T> operator+(const Location<T> right)
	{
	  return Location(x + right.x,
			  y + right.y,
			  z + right.z);
	}

	//Vector subraction
	Location<T> operator-(const Location<T> right)
	{
	  return Location(x - right.x,
			  y - right.y,
			  z - right.z);
	}
	
	//Scalar multiplication
	template < class MultiplierT >
	Location<T> operator*(const MultiplierT multiplier)
	{
	  return Location(x * multiplier,
			  y * multiplier,
			  z * multiplier);
	}

	static Location<T> MaxLimit() 
	{
	  return Location(std::numeric_limits<T>::max());
	}

	static Location<T> MinLimit()
	{
	  return Location(std::numeric_limits<T>::min());
	}

	//Updates the Location in the first Location paramter with the smallest of each
	//of the x, y and z co-ordinates independently of both Locations
	static void UpdateMinLocation(Location<T>& io_store_location, const Location<T>& i_compare_location)
	{
	  io_store_location.x = 
	  util::NumericalFunctions::min
	  (io_store_location.x, 
	   i_compare_location.x);

	  io_store_location.y = 
	  util::NumericalFunctions::min
	  (io_store_location.y, 
	   i_compare_location.y);

	  io_store_location.z = 
	  util::NumericalFunctions::min
	  (io_store_location.z, 
	   i_compare_location.z);
	}
      
	//Updates the Location in the first Location paramter with the largest of each
	//of the x, y and z co-ordinates independently of both Locations
	static void UpdateMaxLocation(Location<T>& io_store_location, const Location<T>& i_compare_location)
	{
	  io_store_location.x = 
	  util::NumericalFunctions::max
	  (io_store_location.x, 
	   i_compare_location.x);

	  io_store_location.y = 
	  util::NumericalFunctions::max
	  (io_store_location.y, 
	   i_compare_location.y);

	  io_store_location.z = 
	  util::NumericalFunctions::max
	  (io_store_location.z, 
	   i_compare_location.z);
	}
	
	};
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_LOCATION_H
