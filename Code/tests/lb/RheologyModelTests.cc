// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.


#include <catch2/catch.hpp>
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace tests
  {
    using namespace lb::kernels::rheologyModels;

    // Tests for the functionality of the non-Newtonian rheology models.
    //
    // The idea here is that I ran the models for a wide range of
    // shear-rates, plotted the results and compared against the
    // literature to a good agreement. Here we just test each of the
    // models against a few known values . It should be enough to
    // detect bugs introduced in any of the components involved.
    template<class RHEO>
    void CompareModelAgainsHardcodedValues(RHEO&& rheo,
					   const distribn_t density,
					   const std::vector<distribn_t>& shearRates,
					   const std::vector<distribn_t>& truthViscs,
					   const std::string& modelName) {
      REQUIRE(shearRates.size() == truthViscs.size());

      for (auto i = 0U; i < shearRates.size(); ++i) {
	distribn_t viscosity = rheo.CalculateViscosityForShearRate(shearRates[i],
								   density);
	INFO("Wrong " << modelName << " viscosity for shear rate " << shearRates[i]);
	CHECK(Approx(truthViscs[i]).margin(1e-10) == viscosity);
      }
    }

    TEST_CASE("RheologyModelTests") {
      const distribn_t density = 1.0;
      const std::vector<distribn_t> shearRates{1e-4, 1e-2, 1, 1e2, 1e4};
      // Ensure we are using coarsegrained blood parameters
      const lb::LbmParameters lbp{1.0, 1.0, DEFAULT_FLUID_DENSITY_Kg_per_m3, 0.004};
      lb::kernels::InitParams ip;
      ip.lbmParams = &lbp;

      SECTION("CarreauYasuda") {
	auto cy = CarreauYasudaRheologyModelHumanFit{ip};
	const std::vector<distribn_t> carreauViscosities{
	  0.1579855, 0.1283352, 2.597279780e-02, 4.282559500e-03, 3.521182515e-03
        };
	CompareModelAgainsHardcodedValues(cy,
					  density,
					  shearRates,
					  carreauViscosities,
					  "CarreauYasuda");
      }

      SECTION("Casson") {
	auto casson = CassonRheologyModel{ip};
	const std::vector<distribn_t> cassonViscosities{
	  0.16, 0.16, 6.185169e-02, 5.5308969e-03, 3.241821969e-03};
	CompareModelAgainsHardcodedValues(casson,
					  density,
					  shearRates,
					  cassonViscosities,
					  "Casson");
      }

      SECTION("TruncatedPowerLaw") {
	auto tpl = TruncatedPowerLawRheologyModel{ip};
	const std::vector<distribn_t> powerLawViscosities{
	  0.16, 4e-02, 4e-03, 3.500030078e-03, 3.500030078e-03};
	CompareModelAgainsHardcodedValues(tpl,
					  density,
					  shearRates,
					  powerLawViscosities,
					  "TruncatedPowerLaw");
      }
    }
  }
}

