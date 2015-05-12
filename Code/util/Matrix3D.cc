// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
    distribn_t const * Matrix3D::operator [](const unsigned int row) const
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

    void Matrix3D::timesVector(const util::Vector3D<double>& multiplier,
                               util::Vector3D<double>& result) const
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

    Matrix3D Matrix3D::operator*(Matrix3D const &a) const
    {
      Matrix3D result;
      for(size_t i(0); i < 3; ++i)
      {
        for(size_t j(0); j < 3; ++j)
        {
          result.matrix[i][j] = 0;
          for(size_t k(0); k < 3; ++k)
          {
            result.matrix[i][j] += matrix[i][k] * a.matrix[k][j];
          }
        }
      }
      return result;
    }

    util::Vector3D<double> Matrix3D::operator*(util::Vector3D<double> const &a) const
    {
      util::Vector3D<double> result;
      timesVector(a, result);
      return result;
    }

    std::ostream& operator<<(std::ostream& stream, Matrix3D const &a)
    {
      return stream << a[0][0] << " " << a[0][1] << " " << a[0][2] << "\n"
             << a[1][0] << " " << a[1][1] << " " << a[1][2] << "\n"
             << a[2][0] << " " << a[2][1] << " " << a[2][2];
    }

  } /* namespace util */
} /* namespace hemelb */
