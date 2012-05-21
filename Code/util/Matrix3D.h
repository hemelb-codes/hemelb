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
        distribn_t matrix[3][3];

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
         * @param result matrix-vector product result. result is assumed initialised to 0.
         */
        void timesVector(const util::Vector3D<double>& multiplier, util::Vector3D<double>& result);
    };
  } /* namespace util */
} /* namespace hemelb */
#endif /* HEMELB_UTIL_MATRIX3D_H */
