
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_LATTICETESTS_H
#define HEMELB_UNITTESTS_LBTESTS_LATTICETESTS_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include "lb/lattices/D3Q15.h"
#include "lb/lattices/D3Q19.h"
#include "lb/lattices/D3Q27.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      class LatticeTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (LatticeTests);
          CPPUNIT_TEST (TestD3Q15);
          CPPUNIT_TEST (TestD3Q19);
          CPPUNIT_TEST (TestD3Q27);CPPUNIT_TEST_SUITE_END();

        public:

          LatticeTests() :
              epsilon(1e-10)
          {
          }
          // had to move this here for compilation portability, see http://stackoverflow.com/questions/370283/why-cant-i-have-a-non-integral-static-const-member-in-a-class
          // wasn't compiling under cray compilers)

          void setUp()
          {
          }

          void tearDown()
          {

          }

          void TestD3Q15()
          {
            TestLattice<lb::lattices::D3Q15>();
          }

          void TestD3Q19()
          {
            TestLattice<lb::lattices::D3Q19>();
          }

          void TestD3Q27()
          {
            TestLattice<lb::lattices::D3Q27>();
          }

        private:
          template<class LatticeType>
          void TestLattice()
          {
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
              CPPUNIT_ASSERT(LatticeType::CX[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CY[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CZ[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CX[direction] <= 1);
              CPPUNIT_ASSERT(LatticeType::CY[direction] <= 1);
              CPPUNIT_ASSERT(LatticeType::CZ[direction] <= 1);

              CPPUNIT_ASSERT(LatticeType::CX[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CYD[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CZD[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CXD[direction] <= 1);
              CPPUNIT_ASSERT(LatticeType::CYD[direction] <= 1);
              CPPUNIT_ASSERT(LatticeType::CZD[direction] <= 1);

              
              for (Direction otherDirection = 0; otherDirection < LatticeType::NUMVECTORS; ++otherDirection)
              {
                if (otherDirection == direction)
                {
                  continue;
                }

                CPPUNIT_ASSERT(LatticeType::CX[direction] != LatticeType::CX[otherDirection]
                    || LatticeType::CY[direction] != LatticeType::CY[otherDirection]
                    || LatticeType::CZ[direction] != LatticeType::CZ[otherDirection]);
                
                CPPUNIT_ASSERT(LatticeType::CXD[direction] != LatticeType::CXD[otherDirection]
                    || LatticeType::CYD[direction] != LatticeType::CYD[otherDirection]
                    || LatticeType::CZD[direction] != LatticeType::CZD[otherDirection]);
              }
            }

            /*
             static const int INVERSEDIRECTIONS[NUMVECTORS];

             Require that inverses are in pairs (this inherently checks for uniqueness, obvs)
             */
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              Direction inverse = LatticeType::INVERSEDIRECTIONS[direction];

              CPPUNIT_ASSERT(direction == LatticeType::INVERSEDIRECTIONS[inverse]);
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

            distribn_t density, momentum[3], expectedDensity, expectedMomentum[3];
            LatticeType::CalculateDensityAndMomentum(f_data, density, momentum[0], momentum[1], momentum[2]);

            LbTestsHelper::CalculateRhoMomentum<LatticeType>(f_data, expectedDensity, expectedMomentum);
            
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedDensity, density, epsilon);
            
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedMomentum[0], momentum[0], epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedMomentum[1], momentum[1], epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedMomentum[2], momentum[2], epsilon);

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
            distribn_t targetDensity = 0.95, targetH[3] = { 0.002, 0.003, 0.004 };

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
                                                          targetH[0],
                                                          targetH[1],
                                                          targetH[2],
                                                          expectedEquilibriumF);

            LbTestsHelper::CalculateAnsumaliEntropicEqmF<LatticeType>(targetDensity,
                                                                      targetH[0],
                                                                      targetH[1],
                                                                      targetH[2],
                                                                      expectedEquilibriumEntropicFAnsumali);

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedEquilibriumF[direction], equilibriumF[direction], epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedEquilibriumEntropicFAnsumali[direction],
                                           equilibriumEntropicFAnsumali[direction],
                                           epsilon);
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

            CPPUNIT_ASSERT_DOUBLES_EQUAL(calculatedDensity, targetDensity, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedDensityAnsumali, targetDensity, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedDensityChikatamarla, targetDensity, epsilon);

            for (Direction direction = 0; direction < 3; direction++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(calculatedMomentum[direction], targetH[direction], epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedMomentumAnsumali[direction], targetH[direction], epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedMomentumChikatamarla[direction],
                                           targetH[direction],
                                           epsilon);
            }

            /*
             * static LatticeInfo* GetLatticeInfo();
             *
             * This must have the same values as the corresponding properties on the main lattice object.
             */
            lb::lattices::LatticeInfo& latticeInfo = LatticeType::GetLatticeInfo();

            CPPUNIT_ASSERT(latticeInfo.GetNumVectors() == LatticeType::NUMVECTORS);

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              CPPUNIT_ASSERT(latticeInfo.GetInverseIndex(direction) == LatticeType::INVERSEDIRECTIONS[direction]);

              const util::Vector3D<int>& velocityVector = latticeInfo.GetVector(direction);

              for (Direction index = 0; index < 3; index++)
              {
                CPPUNIT_ASSERT(velocityVector[index] == LatticeType::discreteVelocityVectors[index][direction]);
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

              CPPUNIT_ASSERT_EQUAL(traction[0], (density - 1) * Cs2);
              CPPUNIT_ASSERT_EQUAL(traction[1], 0.0);
              CPPUNIT_ASSERT_EQUAL(traction[2], 0.0);

              util::Vector3D<LatticeStress> tangentialComponentTraction;
              LatticeType::CalculateTangentialProjectionTraction(density,
                                                                 1.0,
                                                                 nonEquilibriumF.data(),
                                                                 wallNormal,
                                                                 tangentialComponentTraction);

              CPPUNIT_ASSERT_EQUAL(tangentialComponentTraction[0], 0.0);
              CPPUNIT_ASSERT_EQUAL(tangentialComponentTraction[1], 0.0);
              CPPUNIT_ASSERT_EQUAL(tangentialComponentTraction[2], 0.0);

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
                      CPPUNIT_ASSERT_EQUAL( (density - 1) * Cs2 + 1 - 1 / (2 * tau),
                                           stressTensor[rowIndex][columnIndex]);
                    }
                    else
                    {
                      CPPUNIT_ASSERT_EQUAL( (density - 1) * Cs2, stressTensor[rowIndex][columnIndex]);
                    }
                  }
                  else
                  {
                    CPPUNIT_ASSERT_EQUAL(0.0, stressTensor[rowIndex][columnIndex]);
                  }
                }
              }
            }

          }
          const double epsilon;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (LatticeTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_LATTICETESTS_H */
