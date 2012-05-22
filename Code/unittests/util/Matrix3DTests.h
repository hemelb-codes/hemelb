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
          CPPUNIT_TEST(TestMatrix3D);CPPUNIT_TEST_SUITE_END();

          double ComputeMatrixElement(unsigned row, unsigned column)
          {
            return row*3.0 + (column + 1.0);
          }

        public:
          void TestMatrix3D()
          {
            Matrix3D matrix;

            for (unsigned row = 0; row < 3; row++)
            {
              for (unsigned column = 0; column < 3; column++)
              {
                matrix[row][column] = ComputeMatrixElement(row, column);
              }
            }

            matrix *= 2;

            for (unsigned row = 0; row < 3; row++)
            {
              for (unsigned column = 0; column < 3; column++)
              {
                CPPUNIT_ASSERT_EQUAL(matrix[row][column], 2.0 * ComputeMatrixElement(row, column));
              }
            }

            matrix.addDiagonal(3);

            for (unsigned diag = 0; diag < 3; diag++)
            {
              CPPUNIT_ASSERT_EQUAL(matrix[diag][diag], 3.0 + 2.0 * ComputeMatrixElement(diag, diag));
            }

            Vector3D<double> multiplierVector(1,2,3);
            Vector3D<double> resultVector;
            matrix.timesVector(multiplierVector, resultVector);

            CPPUNIT_ASSERT_EQUAL(resultVector[0], 31.0);
            CPPUNIT_ASSERT_EQUAL(resultVector[1], 70.0);
            CPPUNIT_ASSERT_EQUAL(resultVector[2], 109.0);
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(Matrix3DTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_UTIL_MATRIX3DTESTS_H */
