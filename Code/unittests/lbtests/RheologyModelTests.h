
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_RHEOLOGYMODELTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_RHEOLOGYMODELTESTS_H

#include <cppunit/TestFixture.h>
#include "lb/kernels/rheologyModels/RheologyModels.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {

      /**
       * Class containing tests for the functionality of the non-Newtonian rheology models.
       *
       * The idea here is that I ran the models for a wide range of shear-rates, plotted the
       * results and compared against the literature to a good agreement. Here we just test
       * each of the models against a few known values . It should be enough to detect bugs
       * introduced in any of the components involved.
       */
      class RheologyModelTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(RheologyModelTests);
          CPPUNIT_TEST(TestRheologyModels);CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            density = 1.0;

            shearRates.push_back(1e-4);
            shearRates.push_back(1e-2);
            shearRates.push_back(1);
            shearRates.push_back(1e2);
            shearRates.push_back(1e4);

            carreauViscosities.push_back(0.0001579855);
            carreauViscosities.push_back(0.0001283352);
            carreauViscosities.push_back(2.597279780e-05);
            carreauViscosities.push_back(4.282559500e-06);
            carreauViscosities.push_back(3.521182515e-06);

            cassonViscosities.push_back(0.00016);
            cassonViscosities.push_back(0.00016);
            cassonViscosities.push_back(6.185169e-05);
            cassonViscosities.push_back(5.5308969e-06);
            cassonViscosities.push_back(3.241821969e-06);

            powerLawViscosities.push_back(0.00016);
            powerLawViscosities.push_back(4e-05);
            powerLawViscosities.push_back(4e-06);
            powerLawViscosities.push_back(3.500030078e-06);
            powerLawViscosities.push_back(3.500030078e-06);
          }

          void tearDown()
          {
          }

          template<class RHEOLOGY_MODEL>
          void CompareModelAgainsHardcodedValues(const std::vector<distribn_t>& shearRates,
                                                 const std::vector<distribn_t>& truthViscosities,
                                                 const std::string& modelName) const
          {
            CPPUNIT_ASSERT_EQUAL(shearRates.size(), truthViscosities.size());

            std::vector<distribn_t>::const_iterator shearRate;
            std::vector<distribn_t>::const_iterator truthVis;
            for (shearRate = shearRates.begin(), truthVis = truthViscosities.begin();
                shearRate != shearRates.end(); ++shearRate, ++truthVis)
            {
              distribn_t viscosity = RHEOLOGY_MODEL::CalculateViscosityForShearRate(*shearRate,
                                                                                    density);

              std::stringstream message;
              message << "Wrong " << modelName << " viscosity for shear rate " << *shearRate;
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(message.str(), *truthVis, viscosity, 1e-10);
            }
          }

          void TestRheologyModels()
          {
            CompareModelAgainsHardcodedValues<
                lb::kernels::rheologyModels::CarreauYasudaRheologyModelHumanFit>(shearRates,
                                                                                 carreauViscosities,
                                                                                 "CarreauYasuda");
            CompareModelAgainsHardcodedValues<lb::kernels::rheologyModels::CassonRheologyModel>(shearRates,
                                                                                                cassonViscosities,
                                                                                                "Casson");
            CompareModelAgainsHardcodedValues<
                lb::kernels::rheologyModels::TruncatedPowerLawRheologyModel>(shearRates,
                                                                             powerLawViscosities,
                                                                             "TruncatedPowerLaw");
          }

        private:

          distribn_t density;
          std::vector<distribn_t> shearRates;
          std::vector<distribn_t> carreauViscosities;
          std::vector<distribn_t> cassonViscosities;
          std::vector<distribn_t> powerLawViscosities;

      };
      CPPUNIT_TEST_SUITE_REGISTRATION(RheologyModelTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_RHEOLOGYMODELTESTS_H */
