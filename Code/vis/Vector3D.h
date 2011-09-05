#ifndef HEMELB_VIS_RAYTRACER_VECTOR3D_H
#define HEMELB_VIS_RAYTRACER_VECTOR3D_H

#include "constants.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace vis
  {
    // Vector3D represents a 3D vector, storing the
    // x, y and z co-ordinate in the templated numeric type
    // Other methods are defined for convenience
    template<class T = site_t>
    class Vector3D
    {
      public:
        T x, y, z;

        Vector3D()
        {
        }
        ;

        Vector3D(T iX, T iY, T iZ) :
            x(iX), y(iY), z(iZ)
        {
        }

        Vector3D(T iX) :
            x(iX), y(iX), z(iX)
        {
        }

        //Copy constructor - can be used to perform type converstion
        template<class OldTypeT>
        Vector3D<T>(const Vector3D<OldTypeT> & iOldVector3D)
        {
          x = static_cast<T>(iOldVector3D.x);
          y = static_cast<T>(iOldVector3D.y);
          z = static_cast<T>(iOldVector3D.z);
        }

        //Equality
        bool operator==(const Vector3D<T> right)
        {
          if (x != right.x)
          {
            return false;
          }
          if (y != right.y)
          {
            return false;
          }
          if (z != right.z)
          {
            return false;
          }

          return true;
        }

        //Vector addition
        Vector3D<T> operator+(const Vector3D<T> right) const
        {
          return Vector3D(x + right.x, y + right.y, z + right.z);
        }

        //Vector addition
        void operator+=(const Vector3D<T> right)
        {
          x += right.x;
          y += right.y;
          z += right.z;
        }

        //Vector subraction
        Vector3D<T> operator-(const Vector3D<T> right) const
        {
          return Vector3D(x - right.x, y - right.y, z - right.z);
        }

        //Scalar multiplication
        template<class MultiplierT>
        Vector3D<T> operator*(const MultiplierT multiplier) const
        {
          return Vector3D(x * multiplier, y * multiplier, z * multiplier);
        }

        static Vector3D<T> MaxLimit()
        {
          return Vector3D(std::numeric_limits<T>::max());
        }

        static Vector3D<T> MinLimit()
        {
          return Vector3D(std::numeric_limits<T>::min());
        }

        //Updates the Vector3D in the first Vector3D paramter with the smallest of each
        //of the x, y and z co-ordinates independently of both Vector3Ds
        static void UpdateMinVector3D(Vector3D<T>& io_store_location,
                                      const Vector3D<T>& i_compare_location)
        {
          io_store_location.x = util::NumericalFunctions::min(io_store_location.x,
                                                              i_compare_location.x);

          io_store_location.y = util::NumericalFunctions::min(io_store_location.y,
                                                              i_compare_location.y);

          io_store_location.z = util::NumericalFunctions::min(io_store_location.z,
                                                              i_compare_location.z);
        }

        //Updates the Vector3D in the first Vector3D paramter with the largest of each
        //of the x, y and z co-ordinates independently of both Vector3Ds
        static void UpdateMaxVector3D(Vector3D<T>& io_store_location,
                                      const Vector3D<T>& i_compare_location)
        {
          io_store_location.x = util::NumericalFunctions::max(io_store_location.x,
                                                              i_compare_location.x);

          io_store_location.y = util::NumericalFunctions::max(io_store_location.y,
                                                              i_compare_location.y);

          io_store_location.z = util::NumericalFunctions::max(io_store_location.z,
                                                              i_compare_location.z);
        }

    };

    template<class T>
    Vector3D<T> operator+(const Vector3D<T> left, const Vector3D<T> right)
    {
      return left + right;
    }

    template<class T>
    Vector3D<T> operator-(const Vector3D<T> left, const Vector3D<T> right)
    {
      return left - right;
    }

    template<class T>
    Vector3D<T> operator*(const T left, const Vector3D<T> right)
    {
      return right * left;
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_VECTOR3D_H
