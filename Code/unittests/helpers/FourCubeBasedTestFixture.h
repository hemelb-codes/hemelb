
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H
#include <cppunit/TestFixture.h>
#include "configuration/SimConfig.h"
#include "lb/collisions/Collisions.h"
#include "lb/SimulationState.h"
#include "net/IOCommunicator.h"
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
          FourCubeBasedTestFixture() :
            initParams()
          {
          }
          void setUp()
          {
            helpers::FolderTestFixture::setUp();

            latDat = FourCubeLatticeData::Create(Comms());
            std::string path("");
            simConfig = new OneInOneOutSimConfig(path);
            simState = new hemelb::lb::SimulationState(simConfig->GetTimeStepLength(),
                                                       simConfig->GetTotalTimeSteps());
            lbmParams = new lb::LbmParameters(simState->GetTimeStepLength(),
                                              simConfig->GetVoxelSize());
            unitConverter = &simConfig->GetUnitConverter();

            initParams.latDat = latDat;
            initParams.siteCount = initParams.latDat->GetLocalFluidSiteCount();
            initParams.lbmParams = lbmParams;
            numSites = initParams.latDat->GetLocalFluidSiteCount();
          }

          void tearDown()
          {
            delete latDat;
            delete lbmParams;
            delete simState;
            delete simConfig;
            helpers::FolderTestFixture::tearDown();
          }

        protected:
          FourCubeLatticeData* latDat;
          lb::kernels::InitParams initParams;
          site_t numSites;
          lb::LbmParameters* lbmParams;
          configuration::SimConfig* simConfig;
          lb::SimulationState* simState;
          const util::UnitConverter* unitConverter;
        private:

      };
    }
  }
}
#endif // ONCE
