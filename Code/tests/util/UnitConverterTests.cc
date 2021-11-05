// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "util/UnitConverter.h"
#include "util/Matrix3D.h"
#include "constants.h"
#include "lb/lattices/Lattice.h"
#include "lb/lattices/D3Q15.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::util;

    TEST_CASE("UnitConverterTests") {
      auto unitConverter = util::UnitConverter(1., 1., Vector3D<double>(0.), 1000.0, 0.0);
      PhysicalPressure pressMmHg = 81.0;
      LatticeDensity densityLatt = unitConverter.ConvertPressureToLatticeUnits(pressMmHg) / Cs2;
      distribn_t tau = 0.5;
      const double epsilon = 1e-9;
      auto apprx = [&](double x) {
	return Approx(x).margin(epsilon);
      };

      SECTION("TestPressure") {
	REQUIRE(apprx(pressMmHg) ==
		unitConverter.ConvertPressureToPhysicalUnits(densityLatt * Cs2));
      }

      SECTION("TestSimpleStressTensor") {
	std::vector<distribn_t> fNonEquilibrium(lb::lattices::D3Q15::NUMVECTORS, 0.0);

	util::Matrix3D stressTensor;
	lb::lattices::D3Q15::CalculateStressTensor(densityLatt, tau, fNonEquilibrium.data(), stressTensor);
	util::Matrix3D stressTensorPhys = unitConverter.ConvertFullStressTensorToPhysicalUnits(stressTensor);

	REQUIRE(apprx(pressMmHg) == stressTensorPhys[0][0] / mmHg_TO_PASCAL);
	REQUIRE(apprx(pressMmHg) == stressTensorPhys[1][1] / mmHg_TO_PASCAL);
	REQUIRE(apprx(pressMmHg) == stressTensorPhys[2][2] / mmHg_TO_PASCAL);
	REQUIRE(apprx((distribn_t) 0.) == stressTensorPhys[1][0]);
      }

      SECTION("TestSimpleTractionVector") {
	std::vector<distribn_t> fNonEquilibrium(lb::lattices::D3Q15::NUMVECTORS, 0.0);
	util::Vector3D<Dimensionless> wallNormal(0.0);
	wallNormal[0] = 1.0;

	util::Vector3D<LatticeStress> traction;
	lb::lattices::D3Q15::CalculateTractionOnAPoint(densityLatt,
						       tau,
						       fNonEquilibrium.data(),
						       wallNormal,
						       traction);
	util::Vector3D<PhysicalStress> tractionPhys = unitConverter.ConvertTractionToPhysicalUnits(traction,
												   wallNormal);

	REQUIRE(apprx(pressMmHg) == tractionPhys[0] / mmHg_TO_PASCAL);
	REQUIRE(apprx(0.0) == tractionPhys[1] / mmHg_TO_PASCAL);
	REQUIRE(apprx(0.0) == tractionPhys[2] / mmHg_TO_PASCAL);
      }

    }
  }
}
