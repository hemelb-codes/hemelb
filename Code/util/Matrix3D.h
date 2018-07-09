
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
        void timesVector(const util::Vector3D<double>& multiplier, util::Vector3D<double>& result) const;

        /**
         * Computes the product of a matrix and a scalar.
         *
         * @param scalarValue
         * @return matrix-scalar product
         */
        Matrix3D operator*(distribn_t scalarValue) const;

      private:

        //! Internal data representation
        distribn_t matrix[3][3];

    };
  } /* namespace util */
} /* namespace hemelb */
#endif /* HEMELB_UTIL_MATRIX3D_H */
