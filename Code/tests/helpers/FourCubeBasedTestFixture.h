// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H
#define HEMELB_TESTS_HELPERS_FOURCUBEBASEDTESTFIXTURE_H

#include <iostream>
#include <memory>

#include <catch2/catch.hpp>

#include "configuration/SimConfig.h"
#include "lb/collisions/Collisions.h"
#include "lb/SimulationState.h"
#include "net/IOCommunicator.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/FolderTestFixture.h"


namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      class FourCubeBasedTestFixtureBase : public FolderTestFixture {

      public:
	FourCubeBasedTestFixtureBase(int cubesize);
	~FourCubeBasedTestFixtureBase();

      protected:
	FourCubeLatticeData* latDat;
	lb::kernels::InitParams initParams;
	site_t numSites;
	lb::LbmParameters* lbmParams;
	configuration::SimConfig* simConfig;
	std::unique_ptr<lb::SimulationState> simState;
	const util::UnitConverter* unitConverter;
	int cubeSize;
	int cubeSizeWithHalo;
      private:
	std::string path;
      };

      template <int CUBESIZE = 4>
      class FourCubeBasedTestFixture : public FourCubeBasedTestFixtureBase {
      public:
	FourCubeBasedTestFixture() : FourCubeBasedTestFixtureBase(CUBESIZE) {
	}
      };
    }
  }
}
#endif // ONCE
