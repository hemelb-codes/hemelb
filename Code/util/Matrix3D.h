// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UTIL_MATRIX3D_H
#define HEMELB_UTIL_MATRIX3D_H

#include <cassert>
#include "units.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace util
  {

    class Matrix3D
    {
      public:
        /**
         * Convenience accessor.
         *
         * @param row
         * @return
         */
        distribn_t* operator [](const unsigned int row);
        distribn_t const * operator [](const unsigned int row) const;

        /**
         * Multiplies all the entries of the matrix by a given value
         *
         * @param value multiplier value
         */
        void operator*=(distribn_t value);

        /**
         * Adds a given value to all the entries along the diagonal:
         *    matrix = matrix + value*I
         * where I is the identity matrix.
         *
         * @param value value to be added to the matrix diagonal
         */
        void addDiagonal(distribn_t value);

        /**
         * Computes the matrix-vector product:
         *    result =  matrix*multiplier
         *
         * @param multiplier vector to be multiplied by the matrix
         * @param result matrix-vector product result.
         */
        void timesVector(const util::Vector3D<double>& multiplier,
                         util::Vector3D<double>& result) const;

        /**
         * Computes the product of a matrix and a scalar.
         *
         * @param scalarValue
         * @return matrix-scalar product
         */
        Matrix3D operator*(distribn_t scalarValue) const;

        Matrix3D operator*(Matrix3D const &a) const;
        util::Vector3D<double> operator*(util::Vector3D<double>  const &a) const;

      private:

        //! Internal data representation
        distribn_t matrix[3][3];

    };
    std::ostream& operator<<(std::ostream& stream, Matrix3D const &a);

    //! Rotation matrix for a to b
    //! Rotation matrix maps one basis onto another:
    //!    - a0 (normalised a) maps to b0 (normalised b)
    //!    - a0.Cross(b0) maps to a0.Cross(b0)
    //!    - a0.Cross(a0, Cross(b0)) maps to b0.Cross(a0.Cross(b0))
    util::Matrix3D rotationMatrix(LatticePosition const& a, LatticePosition const& b);

    //! Rotation matrix around axis by theta
    util::Matrix3D rotationMatrix(LatticePosition const& a, Dimensionless const& b);

  } /* namespace util */
} /* namespace hemelb */
#endif /* HEMELB_UTIL_MATRIX3D_H */
