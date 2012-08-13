#ifndef HEMELB_UNITTESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H
#include <cppunit/TestFixture.h>
#include "configuration/SimConfig.h"
#include "lb/collisions/Collisions.h"
#include "lb/SimulationState.h"
#include "topology/NetworkTopology.h"
#include "unittests/FourCubeLatticeData.h"
#include "unittests/OneInOneOutSimConfig.h"
#include "unittests/helpers/FolderTestFixture.h"

#include <iostream>
namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class FourCubeBasedTestFixture : public helpers::FolderTestFixture
      {


        public:
          FourCubeBasedTestFixture():initParams(){}
          void setUp()
          {
            // Initialise the network topology (necessary for using the inlets and oulets.
            int args = 1;
            char** argv = NULL;
            bool success;
            topology::NetworkTopology::Instance()->Init(args, argv, &success);
            latDat = FourCubeLatticeData::Create();

            simConfig = new OneInOneOutSimConfig();
            simState = new hemelb::lb::SimulationState(simConfig->GetTimeStepLength(),simConfig->GetTotalTimeSteps());
            lbmParams = new lb::LbmParameters(simState->GetTimeStepLength(),
                                              latDat->GetVoxelSize());
            unitConverter = new util::UnitConverter(lbmParams, simState, latDat->GetVoxelSize());

            initParams.maximumDensityGradient = 0.5;
            initParams.latDat = latDat;
            initParams.siteCount = initParams.latDat->GetLocalFluidSiteCount();
            initParams.lbmParams=lbmParams;
            numSites = initParams.latDat->GetLocalFluidSiteCount();
          }

          void tearDown()
          {
            delete latDat;
            delete lbmParams;
            delete simState;
            delete simConfig;
            delete unitConverter;
          }

        protected:
          FourCubeLatticeData* latDat;
          lb::kernels::InitParams initParams;
          site_t numSites;
          lb::LbmParameters* lbmParams;
          configuration::SimConfig* simConfig;
          lb::SimulationState* simState;
          util::UnitConverter* unitConverter;
        private:

      };
    }
  }
}
#endif // ONCE
