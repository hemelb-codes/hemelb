#ifndef HEMELB_VIS_RAYTRACER_XYCOORDINATE_H
#define HEMELB_VIS_RAYTRACER_XYCOORDINATE_H

#include "constants.h"
#include "util/utilityFunctions.h"


namespace hemelb
{
  namespace vis
  {
    template <class T> 
      class XYCoordinates
      {
      public:
      T x, y;

      XYCoordinates() {};

      XYCoordinates(T iX, T iY) :
      x(iX), y(iY)
      {}
      
      XYCoordinates(T iN) :
      x(iN), y(iN)
      {}

      //Copy constructor - can be used to perform type converstion 
      template < class OldTypeT >
      XYCoordinates<T>(const XYCoordinates<OldTypeT> & iOldXYCoordinates)
      {
	x = static_cast<T>(iOldXYCoordinates.x);
	y = static_cast<T>(iOldXYCoordinates.y);
      }

      //Equality
      bool operator==(const XYCoordinates<T> right) 
      {
	if(x != right.x) { return false; }
	if(y != right.y) { return false; }
	return true;
      }
      
	
      //Vector addition
      XYCoordinates<T> operator+(const XYCoordinates<T> right) const
      {
	return XYCoordinates(x + right.x,
			   y + right.y);
      }

      //Vector addition
      void operator+=(const XYCoordinates<T> right) 
      {
	x += right.x;
	y += right.y;
      }
      

      //Vector subtraction
      XYCoordinates<T> operator-(const XYCoordinates<T> right) const
      {
	return XYCoordinates(x - right.x,
			   y - right.y);
      }
	
      //Scalar multiplication
      template < class MultiplierT >
      XYCoordinates<T> operator*(const MultiplierT multiplier) const
      {
	return XYCoordinates(x * multiplier,
			   y * multiplier);
      }

      static XYCoordinates<T> MaxLimit() 
      {
	return XYCoordinates(std::numeric_limits<T>::max());
      }

      static XYCoordinates<T> MinLimit()
      {
	return XYCoordinates(std::numeric_limits<T>::min());
      }

      //Updates the XYCoordinates in the first XYCoordinates paramter with the smallest of each
      //of the x and y co-ordinatess independently of both XYCoordinates
      static void UpdateMinXYCoordinates(XYCoordinates<T>& io_store_location, const XYCoordinates<T>& i_compare_location)
      {
	io_store_location.x = 
	util::NumericalFunctions::min
	(io_store_location.x, 
	 i_compare_location.x);

	io_store_location.y = 
	util::NumericalFunctions::min
	(io_store_location.y, 
	 i_compare_location.y);
      } 
      
      //Updates the XYCoordinates in the first XYCoordinate paramter with the largest of each
      //of the x and y co-ordinates independently of both XYCoordinates
      static void UpdateMaxXYCoordinates(XYCoordinates<T>& io_store_location, const XYCoordinates<T>& i_compare_location)
      {
	io_store_location.x = 
	util::NumericalFunctions::max
	(io_store_location.x, 
	 i_compare_location.x);

	io_store_location.y = 
	util::NumericalFunctions::max
	(io_store_location.y, 
	 i_compare_location.y);
      }
	
      };

    template <class T>
    XYCoordinates<T> operator+(const XYCoordinates<T> left, const XYCoordinates <T> right) 
    {
      return left + right;
    }
   
    template <class T>
    XYCoordinates<T> operator-(const XYCoordinates<T> left, const XYCoordinates <T> right)
    {
      return left - right;
    }
    
    template <class T>
    XYCoordinates<T> operator*(const T left, const XYCoordinates <T> right)
    {
      return right*left;
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_XYCOORDINATES_H
