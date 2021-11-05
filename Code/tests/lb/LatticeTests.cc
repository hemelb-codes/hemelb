// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "lb/lattices/D3Q15.h"
#include "lb/lattices/D3Q19.h"
#include "lb/lattices/D3Q27.h"

#include "tests/lb/LbTestsHelper.h"

namespace hemelb
{
  namespace tests
  {

    Approx apprx(double x) {
      return Approx(x).margin(1e-10);
    }

    TEMPLATE_TEST_CASE("Lattices work as expected", "[lb]",
		       lb::lattices::D3Q15, lb::lattices::D3Q19, lb::lattices::D3Q27) {
      using LatticeType = TestType;
      /*
	Interface being tested:

	static const unsigned int NUMVECTORS
	static const int CX[NUMVECTORS];
	static const int CY[NUMVECTORS];
	static const int CZ[NUMVECTORS];

	Require that each vector element is in {0,1,-1} and that each vector is unique.
      */

      for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
	{
	  REQUIRE(LatticeType::CX[direction] >= -1);
	  REQUIRE(LatticeType::CY[direction] >= -1);
	  REQUIRE(LatticeType::CZ[direction] >= -1);
	  REQUIRE(LatticeType::CX[direction] <= 1);
	  REQUIRE(LatticeType::CY[direction] <= 1);
	  REQUIRE(LatticeType::CZ[direction] <= 1);

	  REQUIRE(LatticeType::CX[direction] >= -1);
	  REQUIRE(LatticeType::CYD[direction] >= -1);
	  REQUIRE(LatticeType::CZD[direction] >= -1);
	  REQUIRE(LatticeType::CXD[direction] <= 1);
	  REQUIRE(LatticeType::CYD[direction] <= 1);
	  REQUIRE(LatticeType::CZD[direction] <= 1);

              
	  for (Direction otherDirection = 0; otherDirection < LatticeType::NUMVECTORS; ++otherDirection)
	    {
	      if (otherDirection == direction)
                {
                  continue;
                }

	      REQUIRE((LatticeType::CX[direction] != LatticeType::CX[otherDirection]
		       || LatticeType::CY[direction] != LatticeType::CY[otherDirection]
		       || LatticeType::CZ[direction] != LatticeType::CZ[otherDirection]));
                
	      REQUIRE((LatticeType::CXD[direction] != LatticeType::CXD[otherDirection]
		       || LatticeType::CYD[direction] != LatticeType::CYD[otherDirection]
		       || LatticeType::CZD[direction] != LatticeType::CZD[otherDirection]));
	    }
	}

      /*
	static const int INVERSEDIRECTIONS[NUMVECTORS];

	Require that inverses are in pairs (this inherently checks for uniqueness, obvs)
      */
      for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
	{
	  Direction inverse = LatticeType::INVERSEDIRECTIONS[direction];

	  REQUIRE(direction == LatticeType::INVERSEDIRECTIONS[inverse]);
	}

      /*
	static void CalculateDensityAndMomentum(const distribn_t f[],
	distribn_t &density,
	distribn_t &v_x,
	distribn_t &v_y,
	distribn_t &v_z);
      */

      distribn_t f_data[LatticeType::NUMVECTORS];

      // The 3 here is essentially a seed that relates to the magnitude of the density.
      LbTestsHelper::InitialiseAnisotropicTestData<LatticeType>(3, f_data);

      distribn_t density, momentum[3], expectedDensity;
      util::Vector3D<distribn_t> expectedMomentum;
      LatticeType::CalculateDensityAndMomentum(f_data, density, momentum[0], momentum[1], momentum[2]);

      LbTestsHelper::CalculateRhoMomentum<LatticeType>(f_data, expectedDensity, expectedMomentum);
            
      REQUIRE(apprx(expectedDensity) == density);
            
      REQUIRE(apprx(expectedMomentum[0]) == momentum[0]);
      REQUIRE(apprx(expectedMomentum[1]) == momentum[1]);
      REQUIRE(apprx(expectedMomentum[2]) == momentum[2]);

      /*
	static void CalculateFeq(const distribn_t &density,
	const distribn_t &xMomentum,
	const distribn_t &yMomentum,
	const distribn_t &zMomentum,
	distribn_t f_eq[]);

	static void CalculateEntropicFeq(const distribn_t &density,
	const distribn_t &v_x,
	const distribn_t &v_y,
	const distribn_t &v_z,
	distribn_t f_eq[]);
      */
      distribn_t equilibriumF[LatticeType::NUMVECTORS], expectedEquilibriumF[LatticeType::NUMVECTORS],
	equilibriumEntropicFAnsumali[LatticeType::NUMVECTORS],
	expectedEquilibriumEntropicFAnsumali[LatticeType::NUMVECTORS],
	equilibriumEntropicFChikatamarla[LatticeType::NUMVECTORS];

      // These values chosen as they're pairwise coprime. Probably doesn't matter.
      distribn_t targetDensity = 0.95;
      util::Vector3D<distribn_t> targetH{ 0.002, 0.003, 0.004 };

      LatticeType::CalculateFeq(targetDensity, targetH[0], targetH[1], targetH[2], equilibriumF);
      LatticeType::CalculateEntropicFeqAnsumali(targetDensity,
						targetH[0],
						targetH[1],
						targetH[2],
						equilibriumEntropicFAnsumali);
      LatticeType::CalculateEntropicFeqChik(targetDensity,
					    targetH[0],
					    targetH[1],
					    targetH[2],
					    equilibriumEntropicFChikatamarla);

      LbTestsHelper::CalculateLBGKEqmF<LatticeType>(targetDensity,
						    targetH,
						    expectedEquilibriumF);

      LbTestsHelper::CalculateAnsumaliEntropicEqmF<LatticeType>(targetDensity,
								targetH,
								expectedEquilibriumEntropicFAnsumali);

      for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
	{
	  REQUIRE(apprx(expectedEquilibriumF[direction]) == equilibriumF[direction]);
	  REQUIRE(apprx(expectedEquilibriumEntropicFAnsumali[direction]) ==
				       equilibriumEntropicFAnsumali[direction]);
	}

      /*
       * It's also the case that these should be invertible (i.e. the density and velocity should be what we started with).
       */
      distribn_t entropicCalculatedDensityAnsumali, entropicCalculatedMomentumAnsumali[3],
	entropicCalculatedDensityChikatamarla, entropicCalculatedMomentumChikatamarla[3], calculatedDensity,
	calculatedMomentum[3];

      LatticeType::CalculateDensityAndMomentum(equilibriumF,
					       calculatedDensity,
					       calculatedMomentum[0],
					       calculatedMomentum[1],
					       calculatedMomentum[2]);

      LatticeType::CalculateDensityAndMomentum(equilibriumEntropicFAnsumali,
					       entropicCalculatedDensityAnsumali,
					       entropicCalculatedMomentumAnsumali[0],
					       entropicCalculatedMomentumAnsumali[1],
					       entropicCalculatedMomentumAnsumali[2]);

      LatticeType::CalculateDensityAndMomentum(equilibriumEntropicFChikatamarla,
					       entropicCalculatedDensityChikatamarla,
					       entropicCalculatedMomentumChikatamarla[0],
					       entropicCalculatedMomentumChikatamarla[1],
					       entropicCalculatedMomentumChikatamarla[2]);

      REQUIRE(apprx(calculatedDensity) == targetDensity);
      REQUIRE(apprx(entropicCalculatedDensityAnsumali) == targetDensity);
      REQUIRE(apprx(entropicCalculatedDensityChikatamarla) == targetDensity);

      for (Direction direction = 0; direction < 3; direction++)
	{
	  REQUIRE(apprx(calculatedMomentum[direction]) == targetH[direction]);
	  REQUIRE(apprx(entropicCalculatedMomentumAnsumali[direction]) == targetH[direction]);
	  REQUIRE(apprx(entropicCalculatedMomentumChikatamarla[direction]) ==
				       targetH[direction]);
	}

      /*
       * static LatticeInfo* GetLatticeInfo();
       *
       * This must have the same values as the corresponding properties on the main lattice object.
       */
      auto& latticeInfo = LatticeType::GetLatticeInfo();

      REQUIRE(latticeInfo.GetNumVectors() == LatticeType::NUMVECTORS);

      for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
	{
	  REQUIRE(latticeInfo.GetInverseIndex(direction) == LatticeType::INVERSEDIRECTIONS[direction]);

	  const util::Vector3D<int>& velocityVector = latticeInfo.GetVector(direction);

	  for (Direction index = 0; index < 3; index++)
	    {
	      REQUIRE(velocityVector[index] == LatticeType::discreteVelocityVectors[index][direction]);
	    }
	}

      /**
       * @todo: Currently untested
       * * CalculateDensityMomentumFEq (it should just call other functions that *are* tested)
       * * CalculateEntropicDensityMomentumFEq (as above)
       * * CalculateVonMisesStress (probably needs a manually calculated test)
       * * CalculateShearStress (as above)
       * * CalculatePiTensor(as above)
       * * CalculateShearRate (as above)
       *
       */

      /*
	inline static void CalculateForceActingOnAPoint(const distribn_t density,
	const distribn_t tau,
	const distribn_t fNonEquilibrium[],
	const util::Vector3D<double>& wallNormal,
	util::Vector3D<double>& forceActing)
      */
      {
	// Test 1: non equilibrium distribution function equals to 0 and wall normal (1,0,0). Result should be (pressure, 0, 0)
	std::vector<distribn_t> nonEquilibriumF(LatticeType::NUMVECTORS, 0.0);
	LatticeDensity density = 3.0;
	util::Vector3D<Dimensionless> wallNormal(1, 0, 0);
	util::Vector3D<LatticeStress> traction;
	LatticeType::CalculateTractionOnAPoint(density, 1.0, nonEquilibriumF.data(), wallNormal, traction);

	REQUIRE(traction[0] == (density - 1) * Cs2);
	REQUIRE(traction[1] == 0.0);
	REQUIRE(traction[2] == 0.0);

	util::Vector3D<LatticeStress> tangentialComponentTraction;
	LatticeType::CalculateTangentialProjectionTraction(density,
							   1.0,
							   nonEquilibriumF.data(),
							   wallNormal,
							   tangentialComponentTraction);

	REQUIRE(tangentialComponentTraction[0] == 0.0);
	REQUIRE(tangentialComponentTraction[1] == 0.0);
	REQUIRE(tangentialComponentTraction[2] == 0.0);

      }

      /*
	inline static void CalculateStressTensor(const distribn_t density,
	const distribn_t tau,
	const distribn_t fNonEquilibrium[],
	util::Matrix3D& stressTensor)
      */
      {
	/*
	 * Test 1: non equilibrium distribution function equals to 0, except for directions (1,0,0) and (0,1,0) where
	 * it is equal to 1 and directions (1,1,1) and (-1,-1,-1) where it is equal to 2 and -2. Resulting tensor should
	 * be diagonal with (p + 1 - 1/(2*tau)) in the first two entries of the diagonal and p in the last one, where p
	 * is the pressure. The contribution of directions (1,1,1) and (-1,-1,-1) should cancel out.
	 */
	std::vector<distribn_t> nonEquilibriumF(LatticeType::NUMVECTORS, 0.0);
	nonEquilibriumF[1] = 1.0;
	nonEquilibriumF[3] = 1.0;
	nonEquilibriumF[7] = 2.0;
	nonEquilibriumF[8] = -2.0;

	LatticeDensity density = 6.0;
	distribn_t tau = 0.75;

	util::Vector3D<Dimensionless> wallNormal(1.0, 0.0, 0.0);
	util::Matrix3D stressTensor;
	LatticeType::CalculateStressTensor(density, tau, nonEquilibriumF.data(), stressTensor);

	for (unsigned rowIndex = 0; rowIndex < 3; ++rowIndex)
	  {
	    for (unsigned columnIndex = 0; columnIndex < 3; ++columnIndex)
	      {
		if (rowIndex == columnIndex)
                  {
                    if (rowIndex != 2)
		      {
			REQUIRE( (density - 1) * Cs2 + 1 - 1 / (2 * tau)
				 == stressTensor[rowIndex][columnIndex]);
		      }
                    else
		      {
			REQUIRE( (density - 1) * Cs2 == stressTensor[rowIndex][columnIndex]);
		      }
                  }
		else
                  {
                    REQUIRE(0.0 == stressTensor[rowIndex][columnIndex]);
                  }
	      }
	  }
      }
    }

  }
}

