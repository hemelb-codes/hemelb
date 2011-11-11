#ifndef HEMELB_UTIL_VECTOR3D_H
#define HEMELB_UTIL_VECTOR3D_H

#include "constants.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace util
  {
    namespace Direction
    {
      enum Direction
      {
        X,
        Y,
        Z
      };
    }

    //Vector3D is essentially a 3D vector, storing the 
    //x, y and z co-ordinate in the templated numeric type
    //Other methods are defined for convenience
    template<class T>
    class Vector3D
    {
      public:
        T x, y, z;

        Vector3D()
        {
        }

        Vector3D(const T iX, const T iY, const T iZ) :
          x(iX), y(iY), z(iZ)
        {
        }

        Vector3D(const T iX) :
          x(iX), y(iX), z(iX)
        {
        }

        T GetByDirection(Direction::Direction lDirection)
        {
          switch (lDirection)
          {
            case Direction::X:
              return x;
              break;

            case Direction::Y:
              return y;
              break;

              // It can only be Direction::Z, but having the default case prevents a compiler
              // warning (that we might not return anything).
            default:
              return z;
              break;
          }
        }

        //Copy constructor - can be used to perform type conversion
        template<class OldTypeT>
        Vector3D<T> (const Vector3D<OldTypeT> & iOldVector3D)
        {
          x = (T) (iOldVector3D.x);
          y = (T) (iOldVector3D.y);
          z = (T) (iOldVector3D.z);
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
        Vector3D& operator+=(const Vector3D<T> right)
        {
          x += right.x;
          y += right.y;
          z += right.z;

          return *this;
        }

        //Normalisation
        void Normalise()
        {
          T magnitudeSquared = x * x + y * y + z * z;

          if (magnitudeSquared == 0)
          {
            return;
          }

          T lInverseMagnitude = 1.0F / sqrtf(magnitudeSquared);

          x *= lInverseMagnitude;
          y *= lInverseMagnitude;
          z *= lInverseMagnitude;
        }

        //Dot product
        T DotProduct(Vector3D<T> iVector) const
        {
          return (x * iVector.x + y * iVector.y + z * iVector.z);
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

        //Updates the Vector3D in with the smallest of each of the x, y and z
        //co-ordinates independently of both Vector3Ds
        void UpdatePointwiseMin(const Vector3D<T>& iCompareVector)
        {
          x = util::NumericalFunctions::min(x, iCompareVector.x);

          y = util::NumericalFunctions::min(y, iCompareVector.y);

          z = util::NumericalFunctions::min(z, iCompareVector.z);
        }

        //Updates the Vector3D with the largest of each of the x, y and z
        //co-ordinates independently of both Vector3Ds
        void UpdatePointwiseMax(const Vector3D<T>& iCompareVector)
        {
          x = util::NumericalFunctions::max(x, iCompareVector.x);

          y = util::NumericalFunctions::max(y, iCompareVector.y);

          z = util::NumericalFunctions::max(z, iCompareVector.z);
        }

        static Vector3D<T> MaxLimit()
        {
          return Vector3D(std::numeric_limits<T>::max());
        }

        static Vector3D<T> MinLimit()
        {
          return Vector3D(std::numeric_limits<T>::min());
        }

    };

  }
}

#endif // HEMELB_UTIL_VECTOR3D_H
