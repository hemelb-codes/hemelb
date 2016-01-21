
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

    void Matrix3D::addDiagonal(distribn_t value)
    {
      for (unsigned row = 0; row < 3; row++)
      {
        matrix[row][row] += value;
      }
    }

    void Matrix3D::timesVector(const util::Vector3D<double>& multiplier, util::Vector3D<double>& result) const
    {
      for (unsigned row = 0; row < 3; row++)
      {
        result[row] = 0.0;
        for (unsigned column = 0; column < 3; column++)
        {
          result[row] += matrix[row][column] * multiplier[column];
        }
      }
    }

    Matrix3D Matrix3D::operator*(distribn_t scalarValue) const
    {
      Matrix3D returnMatrix;
      for (unsigned rowIndex = 0; rowIndex < 3; ++rowIndex)
      {
        for (unsigned columnIndex = 0; columnIndex < 3; ++columnIndex)
        {
          returnMatrix[rowIndex][columnIndex] = scalarValue * matrix[rowIndex][columnIndex];
        }
      }
      return returnMatrix;
    }


  } /* namespace util */
} /* namespace hemelb */
