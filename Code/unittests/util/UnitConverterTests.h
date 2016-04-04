
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_UNITTESTS_UTIL_UNITCONVERTERTESTS_H
#define HEMELB_UNITTESTS_UTIL_UNITCONVERTERTESTS_H

#include <cppunit/TestFixture.h>
#include "util/UnitConverter.h"
#include "util/Matrix3D.h"
#include "constants.h"
#include "lb/lattices/Lattice.h"
#include "lb/lattices/D3Q15.h"

namespace hemelb
{
  namespace unittests
  {
    namespace util
    {
      using namespace hemelb::util;

      class UnitConverterTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(UnitConverterTests);
          CPPUNIT_TEST(TestSimpleStressTensor);
          CPPUNIT_TEST(TestSimpleTractionVector);CPPUNIT_TEST_SUITE_END();

        public:

          void setUp()
          {
            unitConverter = new util::UnitConverter(1., 1., Vector3D<double>(0.));
            pressMmHg = 81.0;
            densityLatt = unitConverter->ConvertPressureToLatticeUnits(pressMmHg) / Cs2;
            tau = 0.5;
            epsilon = 1e-9;
          }

          void tearDown()
          {
            delete unitConverter;
          }

          void TestPressure()
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(pressMmHg,
                                         unitConverter->ConvertPressureToPhysicalUnits(densityLatt * Cs2),
                                         epsilon);
          }

          void TestSimpleStressTensor()
          {
            std::vector<distribn_t> fNonEquilibrium(lb::lattices::D3Q15::NUMVECTORS, 0.0);

            util::Matrix3D stressTensor;
            lb::lattices::D3Q15::CalculateStressTensor(densityLatt, tau, fNonEquilibrium.data(), stressTensor);
            util::Matrix3D stressTensorPhys = unitConverter->ConvertFullStressTensorToPhysicalUnits(stressTensor);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(pressMmHg, stressTensorPhys[0][0] / mmHg_TO_PASCAL, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(pressMmHg, stressTensorPhys[1][1] / mmHg_TO_PASCAL, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(pressMmHg, stressTensorPhys[2][2] / mmHg_TO_PASCAL, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL((distribn_t) 0., stressTensorPhys[1][0], epsilon);
          }

          void TestSimpleTractionVector()
          {
            std::vector<distribn_t> fNonEquilibrium(lb::lattices::D3Q15::NUMVECTORS, 0.0);
            util::Vector3D<Dimensionless> wallNormal(0.0);
            wallNormal[0] = 1.0;

            util::Vector3D<LatticeStress> traction;
            lb::lattices::D3Q15::CalculateTractionOnAPoint(densityLatt,
                                                           tau,
                                                           fNonEquilibrium.data(),
                                                           wallNormal,
                                                           traction);
            util::Vector3D<PhysicalStress> tractionPhys = unitConverter->ConvertTractionToPhysicalUnits(traction,
                                                                                                        wallNormal);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(pressMmHg, tractionPhys[0] / mmHg_TO_PASCAL, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tractionPhys[1] / mmHg_TO_PASCAL, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tractionPhys[2] / mmHg_TO_PASCAL, epsilon);

          }

        private:
          util::UnitConverter* unitConverter;
          PhysicalPressure pressMmHg;
          LatticeDensity densityLatt;
          distribn_t tau;
          double epsilon;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(UnitConverterTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_UTIL_UNITCONVERTERTESTS_H */
