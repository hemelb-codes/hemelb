// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_VECTOR3D_H
#define HEMELB_UTIL_VECTOR3D_H

#include <cmath>
#include <iterator>
#include <ostream>
#include <algorithm>
#include <limits>
#include "util/utilityFunctions.h"

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
      static_assert(std::is_arithmetic<T>::value, "Vector3D only allowed with arithmetic types");
      public:
        // The type of the element
        using value_type = T;

        // An iterator over elements of this Vector3D
        using iterator = Vector3DIterator<T>;

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

        // Default constructor. It's a value type so this has undefined value when default constructed.
        constexpr Vector3D() = default;

        /**
         * Constructor accepting elements
         * @param x-component
         * @param y-component
         * @param z-component
         */
        constexpr Vector3D(const T iX, const T iY, const T iZ) :
            x(iX), y(iY), z(iZ)
        {
        }

        /**
         * Constructor filling in each component with the supplied argument.
         * @param used for all components
         */
        explicit constexpr Vector3D(const T iX) :
            x(iX), y(iX), z(iX)
        {
        }

        // Copy constructor from Vector3D<this type>
        constexpr Vector3D(const Vector3D&) = default;

        // Move constructor from Vector3D<this type>
        //
        // Note that this doesn't gain us anything, unless T is
        // usefully moveable. We also don't provide a converting move
        // c'tor as that will need a copy even if one of the types
        // involved is moveable.
        constexpr Vector3D(Vector3D&&) = default;

        // Constructor converting elements from another type
        template<class U>
        explicit constexpr Vector3D(const Vector3D<U>& i) :
	    x(T(i.x)), y(T(i.y)), z(T(i.z))
        {
        }

        // Go the other way and create a copy as the supplied type
        template <class U>
        Vector3D<U> as() const {
          return Vector3D<U>{U(x), U(y), U(z)};
        }

        // Copy and move assignment are both explicitly defaulted
        constexpr Vector3D& operator=(const Vector3D&) = default;
        constexpr Vector3D& operator=(Vector3D&&) = default;

        /**
         * Get a component by Direction
         * @param lDirection
         * @return The component
         */
        T& GetByDirection(Direction::Direction lDirection)
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
         * Equality
         */
        constexpr bool operator==(const Vector3D& right) const
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
        constexpr Vector3D GetNormalised() const
        {
          Vector3D normed(*this);
          normed.Normalise();
          return normed;
        }

        // Helper traits
        // Note: standard defines additive (+, -) and multiplicative (*, /, %) categories of arithmetic operators
        // 
        // First, result type of doing T OPERATOR decay(scalar)
        template <typename U>
	using add_result_t = decltype(std::declval<T>() + std::declval<typename std::decay<U>::type>());
        template <typename U>
	using mul_result_t = decltype(std::declval<T>() * std::declval<typename std::decay<U>::type>());
        // Second, predicate whether that type is the same as T
        template <typename U>
	static constexpr bool add_result_same_v = std::is_same<T, add_result_t<U>>::value;
        template <typename U>
	static constexpr bool mul_result_same_v = std::is_same<T, mul_result_t<U>>::value;

        /**
         * Dot product between this vector and another
         * @param other
         * @return
         */
        template <typename U, typename RES = mul_result_t<U>>
        constexpr RES Dot(const Vector3D<U>& otherVector) const
        {
          return (x * otherVector.x + y * otherVector.y + z * otherVector.z);
        }
        /**
         * Dot product between two Vector3Ds.
         * @param V1
         * @param V2
         * @return
         */
        template <typename U, typename RES = mul_result_t<U>>
        constexpr static RES Dot(const Vector3D &V1, const Vector3D<U> &V2)
        {
          return V1.Dot(V2);
        }

        /**
         * Cross product between two Vector3D.
         * @param V1
         * @param V2
         * @return
         */
        template <typename U, typename RES = mul_result_t<U>>
        constexpr static Vector3D<RES> Cross(const Vector3D& V1, const Vector3D<U>& V2)
        {
          return Vector3D<RES>(V1.y * V2.z - V1.z * V2.y,
			       V1.z * V2.x - V1.x * V2.z,
			       V1.x * V2.y - V1.y * V2.x);
        }

        /**
         * Cross product between this vector and another.
         * @param other
         * @return
         */
        template <typename U, typename RES = mul_result_t<U>>
	constexpr Vector3D<RES> Cross(const Vector3D<U>& other) const
        {
          return Cross<U, RES>(*this, other);
        }

        template<class OTHER> Vector3D<OTHER> cast() const
        {
          return Vector3D<OTHER>(*this);
        }

        /**
         * Compute the magnitude squared of the vector
         * @return magnitude**2
         */
        constexpr T GetMagnitudeSquared() const
        {
          return this->Dot(*this);
        }

        /**
         * Compute the magnitude of the vector
         * @return magnitude
         */
        T GetMagnitude() const
        {
          static_assert(!std::is_integral<T>::value, "Vector3D::GetMagnitude only makes sense for floating types");
          return std::sqrt(GetMagnitudeSquared());
        }

        /**
         * Vector addition
         * @param right
         * @return this + right
         */
        template <typename U, typename RES = add_result_t<U>>
        constexpr Vector3D<RES> operator+(const Vector3D<U>& right) const
        {
          return Vector3D<RES>(x + right.x, y + right.y, z + right.z);
        }

        /**
         * In-place vector addition
         * @param right
         * @return the updated vector
         */
        template <typename U, typename = typename std::enable_if<add_result_same_v<U>>::type>
        constexpr Vector3D& operator+=(const Vector3D<U>& right)
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
        constexpr Vector3D operator-() const
        {
          return Vector3D(-x, -y, -z);
        }

        /**
         * Vector subtraction
         * @param right
         * @return this - right
         */
        template <typename U, typename RES = add_result_t<U>>
        constexpr Vector3D<RES> operator-(const Vector3D<U>& right) const
        {
          return Vector3D<RES>(x - right.x, y - right.y, z - right.z);
        }

        /**
         * In-place vector subtraction
         * @param right
         * @return this - right
         */
        template <typename U, typename = typename std::enable_if<add_result_same_v<U>>::type>
        constexpr Vector3D& operator-=(const Vector3D<U>& right)
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
        template<class U,
		 typename = typename std::enable_if<std::is_arithmetic<U>::value>::type,
		 typename RES = mul_result_t<U>>
        friend constexpr Vector3D<RES> operator*(const Vector3D& lhs, const U& rhs)
        {
          return Vector3D<RES>{
	    lhs.x * rhs,
	    lhs.y * rhs,
	    lhs.z * rhs
          };
        }

        // For scalar * vector just swap and call above
        template<class U,
		 typename = typename std::enable_if<std::is_arithmetic<U>::value>::type,
		 typename RES = mul_result_t<U>>
        friend constexpr Vector3D<RES> operator*(const U& lhs, const Vector3D& rhs)
        {
          return rhs*lhs;
        }

        /**
         * In-place multiplication by a scalar
         * @param multiplier
         * @return
         */
        template<class U, typename = typename std::enable_if<mul_result_same_v<U>>::type>
        constexpr Vector3D& operator*=(const U& multiplier)
        {
          x *= multiplier;
          y *= multiplier;
          z *= multiplier;
          return *this;
        }

        // Division by a scalar
        // Return type is a new Vector3D of the type of (T / U)
        template<class U,
		 typename = typename std::enable_if<std::is_arithmetic<U>::value>::type,
		 typename RES = mul_result_t<U>>
        constexpr Vector3D<RES> operator/(const U& divisor) const
        {
          return Vector3D<RES>{x / divisor,
			       y / divisor,
			       z / divisor};
        }

        // In-place divison by a scalar
        // Only enabled if the type of T/U == T
        // Returns the updated object
        template<class U, typename = typename std::enable_if<mul_result_same_v<U>>::type>
        constexpr Vector3D& operator/=(const U& divisor)
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
        template<class U, typename RES = mul_result_t<U>>
        constexpr Vector3D<RES> operator%(const U& divisor) const
        {
          return Vector3D<RES>{x % divisor, y % divisor, z % divisor};
        }

        /**
         * In-place scalar modulus
         * @param divisor
         */
        template<class U, typename = typename std::enable_if<mul_result_same_v<U>>::type>
        constexpr Vector3D& operator%=(const U& divisor)
        {
          x %= divisor;
          y %= divisor;
          z %= divisor;
          return *this;
        }

        /**
         * Point-wise multiplication
         */
        constexpr Vector3D PointwiseMultiplication(const Vector3D& rightArgument) const
        {
          return Vector3D(x * rightArgument.x, y * rightArgument.y, z * rightArgument.z);
        }

        /**
         * Point-wise division
         */
        constexpr Vector3D PointwiseDivision(const Vector3D& rightArgument) const
        {
          return Vector3D(x / rightArgument.x, y / rightArgument.y, z / rightArgument.z);
        }

        /**
         * Updates the each component of this Vector3D with the smaller of the
         * corresponding component in this and the other Vector3D.
         *
         * @param vector to compare against
         */
        constexpr void UpdatePointwiseMin(const Vector3D& iCompareVector)
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
        constexpr void UpdatePointwiseMax(const Vector3D& iCompareVector)
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
        constexpr bool IsInRange(const Vector3D& vMin, const Vector3D& vMax) const
        {
          return NumericalFunctions::IsInRange(x, vMin.x, vMax.x)
              && NumericalFunctions::IsInRange(y, vMin.y, vMax.y)
              && NumericalFunctions::IsInRange(z, vMin.z, vMax.z);
        }

        /**
         * Vector filled with the maximum value for the element type.
         * @return
         */
        static constexpr Vector3D MaxLimit()
        {
          return Vector3D(std::numeric_limits<T>::max());
        }

        /**
         * Vector filled with the minimum value for the element type.
         * @return
         */
        static constexpr Vector3D MinLimit()
        {
          return Vector3D(std::numeric_limits<T>::min());
        }

        /**
         * Factory for Vector3Ds of ones.
         * @return
         */
        static constexpr Vector3D Ones()
        {
          return Vector3D(1);
        }

        /**
         * Factory for Vector3Ds of zeros.
         * @return
         */
        static constexpr Vector3D Zero()
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
            vec(nullptr), i(0)
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

#ifdef HEMELB_CODE
#include "net/MpiDataType.h"

namespace hemelb
{
  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<float> >::RegisterMpiDataType();
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<int64_t> >::RegisterMpiDataType();
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<double> >::RegisterMpiDataType();
  }
}
#endif

#endif // HEMELB_UTIL_VECTOR3D_H
