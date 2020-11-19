// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "util/Matrix3D.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::util;

    namespace {
      // Returns the entries of the following 3X3 matrix:
      //     [1 2 3]
      //     [4 5 6]
      //     [7 8 9]
      double ComputeMatrixElement(unsigned row, unsigned column)
      {
	REQUIRE(row<3);
	REQUIRE(column<3);
	return row * 3.0 + (column + 1.0);
      }
    }

    TEST_CASE("Matrix3DTests") {
      // Initialise the matrix with some data
      Matrix3D matrix;
      for (unsigned row = 0; row < 3; row++) {
	for (unsigned column = 0; column < 3; column++) {
	  matrix[row][column] = ComputeMatrixElement(row, column);
	}
      }

      SECTION("Test operator*=") {
	matrix *= 2;
	for (unsigned row = 0; row < 3; row++) {
	  for (unsigned column = 0; column < 3; column++) {
	    REQUIRE(matrix[row][column] == 2.0 * ComputeMatrixElement(row, column));
	  }
	}
      }

      SECTION("Test addDiagonal") {
	matrix.addDiagonal(3);
	for (unsigned diag = 0; diag < 3; diag++) {
	  REQUIRE(matrix[diag][diag] == 3.0 + ComputeMatrixElement(diag, diag));
	}
      }

      SECTION("Test the matrix-vector product") {
	Vector3D<double> multiplierVector(1, 2, 3);
	Vector3D<double> resultVector;
	matrix.timesVector(multiplierVector, resultVector);
	REQUIRE(resultVector[0] == 14.0);
	REQUIRE(resultVector[1] == 32.0);
	REQUIRE(resultVector[2] == 50.0);
      }

      SECTION("Test Matrix3DScalarProduct operator*") {
	auto halfMatrix = matrix * 0.5;
	for (unsigned row = 0; row < 3; row++) {
	  for (unsigned column = 0; column < 3; column++) {
	    REQUIRE(halfMatrix[row][column] == 0.5 * ComputeMatrixElement(row, column));
	  }
	}
      }

    }

  }
}
