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

    util::Matrix3D rotationMatrix(LatticePosition const& a, LatticePosition const& b)
    {
      util::Matrix3D result;
      LatticePosition const a0 = a.GetNormalised();
      LatticePosition const b0 = b.GetNormalised();
      // Special case where a0 == b0: Rotation is identity
      // Also works for a = 0 and b = 0... But that's just bad input
      if (a0.Cross(b0).GetMagnitude() < 1e-8)
      {
        for (size_t i(0); i < 3; ++i)
        {
          for (size_t j(0); j < 3; ++j)
          {
            result[i][j] = i == j ?
              1 :
              0;
          }
        }
        return result;
      }

      LatticePosition const u = (a0.Cross(b0)).GetNormalised();
      LatticePosition const a1 = a0.Cross(u).GetNormalised();
      LatticePosition const b1 = b0.Cross(u).GetNormalised();

      // HemeLB doesn't need linear algebra, so it doesn't have it - cos that would be
      // complicated - which in turn implies it doesn't need it. Instead, do lets use a square
      // wheel.
      util::Matrix3D A, B;
      A[0][0] = a0[0];
      A[0][1] = a0[1];
      A[0][2] = a0[2];
      A[1][0] = u[0];
      A[1][1] = u[1];
      A[1][2] = u[2];
      A[2][0] = a1[0];
      A[2][1] = a1[1];
      A[2][2] = a1[2];
      B[0][0] = b0[0];
      B[0][1] = u[0];
      B[0][2] = b1[0];
      B[1][0] = b0[1];
      B[1][1] = u[1];
      B[1][2] = b1[1];
      B[2][0] = b0[2];
      B[2][1] = u[2];
      B[2][2] = b1[2];
      return B * A;
    }

    util::Matrix3D rotationMatrix(LatticePosition const& axis, Dimensionless const& theta)
    {
      auto const u = axis.GetNormalised();
      Matrix3D result;
      result[0][0] = std::cos(theta) + u.x * u.x * (1 - std::cos(theta));
      result[1][0] = u.y * u.x * (1 - std::cos(theta)) + u.z * std::sin(theta);
      result[2][0] = u.z * u.x * (1 - std::cos(theta)) - u.y * std::sin(theta);

      result[0][1] = u.y * u.x * (1 - std::cos(theta)) - u.z * std::sin(theta);
      result[1][1] = std::cos(theta) + u.y * u.y * (1 - std::cos(theta));
      result[2][1] = u.z * u.y * (1 - std::cos(theta)) + u.x * std::sin(theta);

      result[0][2] = u.x * u.z * (1 - std::cos(theta)) + u.y * std::sin(theta);
      result[1][2] = u.z * u.y * (1 - std::cos(theta)) - u.x * std::sin(theta);
      result[2][2] = std::cos(theta) + u.z * u.z * (1 - std::cos(theta));

      return result;
    }

  } /* namespace util */
} /* namespace hemelb */
