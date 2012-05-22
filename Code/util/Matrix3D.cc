/*
 * Matrix3D.cc
 *
 *  Created on: May 21, 2012
 *      Author: mobernabeu
 */

#include "util/Matrix3D.h"

namespace hemelb
{
  namespace util
  {

    distribn_t* Matrix3D::operator [](const unsigned int row)
    {
      return matrix[row];
    }

    /**
     * Multiplies all the entries of the matrix by a given value
     *
     * @param value multiplier value
     */
    void Matrix3D::operator*=(distribn_t value)
    {
      for (unsigned row = 0; row < 3; row++)
      {
        for (unsigned column = 0; column < 3; column++)
        {
          matrix[row][column] *= value;
        }
      }
    }

    /**
     * Adds a given value to all the entries along the diagonal:
     *    matrix = matrix + value*I
     * where I is the identity matrix.
     *
     * @param value value to be added to the matrix diagonal
     */
    void Matrix3D::addDiagonal(distribn_t value)
    {
      for (unsigned row = 0; row < 3; row++)
      {
        matrix[row][row] += value;
      }
    }

    /**
     * Computes the matrix-vector product:
     *    result =  matrix*multiplier
     *
     * @param multiplier vector to be multiplied by the matrix
     * @param result matrix-vector product result. result is assumed initialised to 0.
     */
    void Matrix3D::timesVector(const util::Vector3D<double>& multiplier, util::Vector3D<double>& result)
    {
      assert(result==0.0);

      for (unsigned row = 0; row < 3; row++)
      {
        for (unsigned column = 0; column < 3; column++)
        {
          result[row] += matrix[row][column] * multiplier[column];
        }
      }
    }

  } /* namespace util */
} /* namespace hemelb */
