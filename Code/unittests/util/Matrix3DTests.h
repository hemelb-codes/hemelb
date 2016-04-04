
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_UTIL_MATRIX3DTESTS_H
#define HEMELB_UNITTESTS_UTIL_MATRIX3DTESTS_H

#include <cppunit/TestFixture.h>
#include "util/Matrix3D.h"

namespace hemelb
{
  namespace unittests
  {
    namespace util
    {
      using namespace hemelb::util;

      class Matrix3DTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(Matrix3DTests);
          CPPUNIT_TEST(TestMatrix3D);
          CPPUNIT_TEST(TestMatrix3DScalarProduct);
          CPPUNIT_TEST_SUITE_END();

          // Returns the entries of the following 3X3 matrix:
          //     [1 2 3]
          //     [4 5 6]
          //     [7 8 9]
          double ComputeMatrixElement(unsigned row, unsigned column)
          {
            CPPUNIT_ASSERT(row<3);
            CPPUNIT_ASSERT(column<3);
            return row * 3.0 + (column + 1.0);
          }

        public:

          void setUp()
          {
            // Initialise the matrix with some data
            for (unsigned row = 0; row < 3; row++)
            {
              for (unsigned column = 0; column < 3; column++)
              {
                matrix[row][column] = ComputeMatrixElement(row, column);
              }
            }
          }

          void TestMatrix3D()
          {
            // Test operator*=
            matrix *= 2;
            for (unsigned row = 0; row < 3; row++)
            {
              for (unsigned column = 0; column < 3; column++)
              {
                CPPUNIT_ASSERT_EQUAL(matrix[row][column], 2.0 * ComputeMatrixElement(row, column));
              }
            }

            // Test addDiagonal
            matrix.addDiagonal(3);
            for (unsigned diag = 0; diag < 3; diag++)
            {
              CPPUNIT_ASSERT_EQUAL(matrix[diag][diag], 3.0 + 2.0 * ComputeMatrixElement(diag, diag));
            }

            // Test the matrix-vector product
            Vector3D<double> multiplierVector(1, 2, 3);
            Vector3D<double> resultVector;
            matrix.timesVector(multiplierVector, resultVector);
            CPPUNIT_ASSERT_EQUAL(resultVector[0], 31.0);
            CPPUNIT_ASSERT_EQUAL(resultVector[1], 70.0);
            CPPUNIT_ASSERT_EQUAL(resultVector[2], 109.0);
          }

          void TestMatrix3DScalarProduct()
          {
            // Test operator*
            Matrix3D halfMatrix = matrix * 0.5;
            for (unsigned row = 0; row < 3; row++)
            {
              for (unsigned column = 0; column < 3; column++)
              {
                CPPUNIT_ASSERT_EQUAL(halfMatrix[row][column], 0.5 * ComputeMatrixElement(row, column));
              }
            }
          }

        private:
          Matrix3D matrix;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(Matrix3DTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_UTIL_MATRIX3DTESTS_H */
