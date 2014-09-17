// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UTIL_VECTOR3D_H
#define HEMELB_UTIL_VECTOR3D_H

#include <cmath>
#include <iterator>
#include <ostream>
#include <algorithm>
#include <limits>
#include "util/static_assert.h"
#include "util/utilityFunctions.h"
#include "util/Vector3DArithmeticTraits.h"

namespace hemelb
{
  namespace util
  {
    namespace Direction
    {
      /**
       * The Cartesian directions.
       */
      enum Direction
      {
        X, //!< X
        Y, //!< Y
        Z
      //!< Z
      };
    }

    /**
     * Base class for Vector3D that handles IndexError.
     *
     * Note that we DO NOT initialise the static member variable
     * HandlerFunction* handler in Vector3D.cc, as you might expect!
     * The reason is that we want to use different handlers within HemeLB
     * and the setup tool. For HemeLB proper, this initialisation is done in
     * Vector3DHemeLb.cc
     */
    class Vector3DBase
    {
      public:
        /**
         * Typedef for the error handler function.
         */
        typedef void (HandlerFunction)(int direction);

        /**
         * Set the function used to handle errors, overriding the default.
         *
         * @param A pointer to the function to use to deal with an index
         * error. Must raise an exception or terminate the simulation or
         * similar.
         */
        static void SetIndexErrorHandler(HandlerFunction* func);
      protected:
        /**
         * Called by subclasses to handle an index error. Wraps the real
         * function used.
         *
         * @param The (out of range) index supplied
         */
        static HandlerFunction HandleIndexError;
      private:
        /**
         * The pointer to the real handler function to use.
         */
        static HandlerFunction* handler;
    };

    // Forward definition of template class fully defined below.
    template<typename T> class Vector3DIterator;

    /**
     * Three dimensional vector class template.
     *
     * It exposes its members publicly for use and overloads arithmetic
     * operators as well as providing a number of useful methods. Other
     * methods are defined for convenience
     */
    template<class T>
    class Vector3D : public Vector3DBase
    {
      public:
        /**
         * The type of the element
         */
        typedef T value_type;

        /**
         * An iterator over elements of this Vector3D
         */
        typedef Vector3DIterator<T> iterator;

        /**
         * Iterator pointing to first element of the Vector3D.
         * @return
         */
        iterator begin()
        {
          return iterator(*this, 0U);
        }

        /**
         * Iterator pointing just past the last element of the Vector3D.
         * @return
         */
        iterator end()
        {
          return iterator(*this, 3U);
        }

        // Be friends with all other Vec3 instantiations.
        // template<class > friend class Vector3D;

        /**
         * x, y and z components.
         */
        T x, y, z;

        /**
         * Default constructor, instantiates with zeros.
         */
        Vector3D() :
            x(0), y(0), z(0)
        {
        }

        /**
         * Constructor accepting elements
         * @param x-component
         * @param y-component
         * @param z-component
         */
        Vector3D(const T iX, const T iY, const T iZ) :
            x(iX), y(iY), z(iZ)
        {
        }

        /**
         * Constructor filling in each component with the supplied argument.
         * @param used for all components
         */
        Vector3D(const T iX) :
            x(iX), y(iX), z(iX)
        {
        }

        /**
         * Get a component by Direction
         * @param lDirection
         * @return The component
         */
        T GetByDirection(Direction::Direction lDirection)
        {
          switch (lDirection)
          {
            case Direction::X:
              return x;

            case Direction::Y:
              return y;

            case Direction::Z:
              return z;

            default:
              HandleIndexError(lDirection);
              // HandleIndexError will either throw or bring the simulation
              // down. We include this return statement to suppress warnings!
              return x;
          }
        }

        /**
         * Get a component by index (x = 0, y = 1, z = 2).
         * @param index
         * @return component
         */
        T& operator[](int index)
        {
          switch (index)
          {
            case 0:
              return x;

            case 1:
              return y;

            case 2:
              return z;

            default:
              HandleIndexError(index);
              // HandleIndexError will either throw or bring the simulation
              // down. We include this return statement to suppress warnings!
              return x;
          }
        }

        /**
         * Get a component by index (x = 0, y = 1, z = 2).
         * @param index
         * @return component
         */
        const T& operator[](int index) const
        {
          switch (index)
          {
            case 0:
              return x;

            case 1:
              return y;

            case 2:
              return z;

            default:
              HandleIndexError(index);
              // HandleIndexError will either throw or bring the simulation
              // down. We include this return statement to suppress warnings!
              return x;
          }
        }

        /**
         * Copy constructor. Can perform type conversion.
         * @param source
         */
        template<class OldTypeT>
        Vector3D(const Vector3D<OldTypeT> & iOldVector3D)
        {
          x = (T) (iOldVector3D.x);
          y = (T) (iOldVector3D.y);
          z = (T) (iOldVector3D.z);
        }

        /**
         * Equality
         */
        bool operator==(const Vector3D<T> right) const
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

        /**
         * In-place normalisation (i.e. this will be a unit vector).
         * @return reference to this
         */
        Vector3D& Normalise()
        {
          (*this) /= GetMagnitude();
          return *this;
        }

        /**
         * Compute the unit vector that points in this direction.
         * @return the unit vector
         */
        Vector3D GetNormalised() const
        {
          Vector3D normed(*this);
          normed.Normalise();
          return normed;
        }

        /**
         * Dot product between this vector and another
         * @param other
         * @return
         */
        T Dot(const Vector3D& otherVector) const
        {
          return (x * otherVector.x + y * otherVector.y + z * otherVector.z);
        }
        /**
         * Dot product between two Vector3Ds.
         * @param V1
         * @param V2
         * @return
         */
        static T Dot(const Vector3D &V1, const Vector3D &V2)
        {
          return V1.Dot(V2);
        }

        /**
         * Cross product between two Vector3D.
         * @param V1
         * @param V2
         * @return
         */
        static Vector3D Cross(const Vector3D& V1, const Vector3D& V2)
        {
          return Vector3D(V1.y * V2.z - V1.z * V2.y,
                          V1.z * V2.x - V1.x * V2.z,
                          V1.x * V2.y - V1.y * V2.x);

        }

        /**
         * Cross product between this vector and another.
         * @param other
         * @return
         */
        Vector3D Cross(const Vector3D& other) const
        {
          return Cross(*this, other);
        }

        /**
         * Compute the magnitude squared of the vector
         * @return magnitude**2
         */
        T GetMagnitudeSquared() const
        {
          return this->Dot(*this);
        }

        /**
         * Compute the magnitude of the vector
         * @return magnitude
         */
        T GetMagnitude() const
        {
          HEMELB_STATIC_ASSERT(std::numeric_limits<T>::is_specialized);
          HEMELB_STATIC_ASSERT(!std::numeric_limits<T>::is_integer);
          return std::sqrt(GetMagnitudeSquared());
        }

        /**
         * Vector addition
         * @param right
         * @return this + right
         */
        Vector3D operator+(const Vector3D right) const
        {
          return Vector3D(x + right.x, y + right.y, z + right.z);
        }

        /**
         * In-place vector addition
         * @param right
         * @return the updated vector
         */
        Vector3D& operator+=(const Vector3D<T> right)
        {
          x += right.x;
          y += right.y;
          z += right.z;

          return *this;
        }

        /**
         * Vector unary negation
         * @return -this
         */
        Vector3D operator-() const
        {
          return Vector3D(-x, -y, -z);
        }

        /**
         * Vector subtraction
         * @param right
         * @return this - right
         */
        Vector3D operator-(const Vector3D right) const
        {
          return Vector3D(x - right.x, y - right.y, z - right.z);
        }

        /**
         * In-place vector subtraction
         * @param right
         * @return this - right
         */
        Vector3D& operator-=(const Vector3D right)
        {
          x -= right.x;
          y -= right.y;
          z -= right.z;
          return *this;
        }

        /**
         * Multiplication by a scalar
         * @param multiplier
         * @return
         */
        template<class MultiplierT>
        Vector3D<typename Vector3DArithmeticTraits<T, MultiplierT>::operatorReturnType> operator*(const MultiplierT multiplier) const
        {
          return Vector3D<typename Vector3DArithmeticTraits<T, MultiplierT>::operatorReturnType>(x * multiplier,
                                                                                                 y * multiplier,
                                                                                                 z * multiplier);
        }

        /**
         * In-place multiplication by a scalar
         * @param multiplier
         * @return
         */
        template<class MultiplierT>
        Vector3D& operator*=(const MultiplierT multiplier)
        {
          x *= multiplier;
          y *= multiplier;
          z *= multiplier;
          return *this;
        }

        /**
         * Division by a scalar
         * @param divisor
         * @return
         */
        template<class DivisorT>
        Vector3D<typename Vector3DArithmeticTraits<T, DivisorT>::operatorReturnType> operator/(const DivisorT divisor) const
        {
          return Vector3D<typename Vector3DArithmeticTraits<T, DivisorT>::operatorReturnType>(x / divisor,
                                                                                              y / divisor,
                                                                                              z / divisor);
        }

        /**
         * In-place divison by a scalar
         * @param divisor
         * @return
         */
        template<class DivisorT>
        Vector3D& operator/=(const DivisorT divisor)
        {
          x /= divisor;
          y /= divisor;
          z /= divisor;
          return *this;
        }

        /**
         * Scalar modulus
         * @param divisor
         */
        template<class ModuloT>
        Vector3D operator%(const ModuloT divisor) const
        {
          return Vector3D(x % divisor, y % divisor, z % divisor);
        }

        /**
         * In-place scalar modulus
         * @param divisor
         */
        template<class ModuloT>
        Vector3D& operator%=(const ModuloT divisor)
        {
          x %= divisor;
          y %= divisor;
          z %= divisor;
          return *this;
        }

        /**
         * Point-wise multiplication
         */
        Vector3D PointwiseMultiplication(const Vector3D& rightArgument) const
        {
          return Vector3D(x * rightArgument.x, y * rightArgument.y, z * rightArgument.z);
        }

        /**
         * Point-wise division
         */
        Vector3D PointwiseDivision(const Vector3D& rightArgument) const
        {
          return Vector3D(x / rightArgument.x, y / rightArgument.y, z / rightArgument.z);
        }

        /**
         * Updates the each component of this Vector3D with the smaller of the
         * corresponding component in this and the other Vector3D.
         *
         * @param vector to compare against
         */
        void UpdatePointwiseMin(const Vector3D& iCompareVector)
        {
          x = std::min(x, iCompareVector.x);

          y = std::min(y, iCompareVector.y);

          z = std::min(z, iCompareVector.z);
        }

        /**
         * Updates the each component of this Vector3D with the larger of the
         * corresponding component in this and the other Vector3D.
         *
         * @param vector to compare against
         */
        void UpdatePointwiseMax(const Vector3D& iCompareVector)
        {
          x = std::max(x, iCompareVector.x);

          y = std::max(y, iCompareVector.y);

          z = std::max(z, iCompareVector.z);
        }

        /**
         * Return whether the vector is within the bounds specified.
         * @param vMin
         * @param vMax
         * @return
         */
        bool IsInRange(const Vector3D& vMin, const Vector3D& vMax) const
        {
          return NumericalFunctions::IsInRange(x, vMin.x, vMax.x)
              && NumericalFunctions::IsInRange(y, vMin.y, vMax.y)
              && NumericalFunctions::IsInRange(z, vMin.z, vMax.z);
        }

        /**
         * Vector filled with the maximum value for the element type.
         * @return
         */
        static Vector3D MaxLimit()
        {
          return Vector3D(std::numeric_limits<T>::max());
        }

        /**
         * Vector filled with the minimum value for the element type.
         * @return
         */
        static Vector3D MinLimit()
        {
          return Vector3D(std::numeric_limits<T>::min());
        }

        /**
         * Factory for Vector3Ds of ones.
         * @return
         */
        static Vector3D Ones()
        {
          return Vector3D(1);
        }

        /**
         * Factory for Vector3Ds of zeros.
         * @return
         */
        static Vector3D Zero()
        {
          return Vector3D(0);
        }
    };

    /**
     * Template class for iterators over the elements of a Vector3D instantiation.
     */
    template<typename T>
    class Vector3DIterator : public std::iterator<std::forward_iterator_tag, T>
    {
      public:
        /**
         * The type of vector over which we will iterate.
         */
        typedef Vector3D<T> vector;
      protected:
        vector* vec; //!< The Vector3D
        unsigned int i; //!< Current position in the vector

      public:
        /**
         * Default constructor
         */
        Vector3DIterator() :
            vec(NULL), i(0)
        {
        }

        /**
         * Construct an iterator over the given vector, starting at the given
         * element
         * @param vector
         * @param element index
         */
        Vector3DIterator(vector& vec, unsigned int i = 0) :
            vec(&vec), i(i)
        {
        }

        /**
         * Copy constructor.
         * @param other
         */
        Vector3DIterator(const Vector3DIterator& other) :
            vec(other.vec), i(other.i)
        {
        }

        /**
         * Assignment
         * @param other
         * @return
         */
        Vector3DIterator& operator=(const Vector3DIterator& other)
        {
          if (this == &other)
          {
            return (*this);
          }
          this->vec = other.vec;
          this->i = other.i;

          return (*this);
        }

        /**
         * Advance to next element.
         * @return the iterator advanced to the next position.
         */
        Vector3DIterator& operator++()
        {
          this->i++;
          return *this;
        }

        /**
         * Test for equality.
         * @param other
         * @return
         */
        bool operator==(const Vector3DIterator& other) const
        {
          return (this->vec == other.vec) && (this->i == other.i);
        }

        /**
         * Test for inequality.
         * @param other
         * @return
         */
        bool operator!=(const Vector3DIterator& other) const
        {
          return ! (*this == other);
        }

        /**
         * Dereference
         * @return element at the current position.
         */
        T& operator*()
        {
          return (*this->vec)[this->i];
        }

        /**
         * Deference
         * @return pointer to element at the current position.
         */
        T* operator->()
        {
          return & (* (*this));
        }
    };

    template<typename T>
    std::ostream& operator<<(std::ostream& o, Vector3D<T> const& v3)
    {
      return o << "(" << v3.x << "," << v3.y << "," << v3.z << ")";
    }

    namespace detail
    {
      bool CheckNextChar(std::istream& i, char c);
    }

    template<typename T>
    std::istream& operator>>(std::istream& i, Vector3D<T>& v3)
    {
      if (detail::CheckNextChar(i, '('))
      {
        i >> v3.x;
        if (detail::CheckNextChar(i, ','))
        {
          i >> v3.y;
          if (detail::CheckNextChar(i, ','))
          {
            i >> v3.z;
            detail::CheckNextChar(i, ')');
          }
        }
      }
      return i;
    }
  }
}

#endif // HEMELB_UTIL_VECTOR3D_H
