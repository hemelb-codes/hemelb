
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.


#include <catch2/catch.hpp>
#include "lb/kernels/rheologyModels/RheologyModels.h"

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
    template<class RHEOLOGY_MODEL>
    void CompareModelAgainsHardcodedValues(const distribn_t density,
					   const std::vector<distribn_t>& shearRates,
					   const std::vector<distribn_t>& truthViscs,
					   const std::string& modelName) {
      REQUIRE(shearRates.size() == truthViscs.size());

      for (auto i = 0U; i < shearRates.size(); ++i) {
	distribn_t viscosity = RHEOLOGY_MODEL::CalculateViscosityForShearRate(shearRates[i],
									      density);
	INFO("Wrong " << modelName << " viscosity for shear rate " << shearRates[i]);
	REQUIRE(Approx(truthViscs[i]).margin(1e-10) == viscosity);
      }
    }

    TEST_CASE("RheologyModelTests") {
      const distribn_t density = 1.0;
      const std::vector<distribn_t> shearRates{1e-4, 1e-2, 1, 1e2, 1e4};

      SECTION("CarreauYasuda") {
	const std::vector<distribn_t> carreauViscosities{
	  0.0001579855, 0.0001283352, 2.597279780e-05, 4.282559500e-06, 3.521182515e-06
        };
	CompareModelAgainsHardcodedValues<CarreauYasudaRheologyModelHumanFit>(density,
									      shearRates,
									      carreauViscosities,
									      "CarreauYasuda");
      }

      SECTION("Casson") {
	const std::vector<distribn_t> cassonViscosities{
	  0.00016, 0.00016, 6.185169e-05, 5.5308969e-06, 3.241821969e-06};
	CompareModelAgainsHardcodedValues<CassonRheologyModel>(density,
							       shearRates,
							       cassonViscosities,
							       "Casson");
      }

      SECTION("TruncatedPowerLaw") {
	const std::vector<distribn_t> powerLawViscosities{
	  0.00016, 4e-05, 4e-06, 3.500030078e-06, 3.500030078e-06};
	CompareModelAgainsHardcodedValues<TruncatedPowerLawRheologyModel>(density,
									  shearRates,
									  powerLawViscosities,
									  "TruncatedPowerLaw");
      }
    }
  }
}

