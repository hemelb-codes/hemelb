//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_GUOFORCING_H
#define HEMELB_UNITTESTS_GUOFORCING_H

#include <cppunit/TestFixture.h>
#include "lb/lattices/D3Q15.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/Comparisons.h"

namespace hemelb { namespace unittests {

  class GuoForcingTests
        : public helpers::FourCubeBasedTestFixture {
    CPPUNIT_TEST_SUITE(GuoForcingTests);
      CPPUNIT_TEST(testHydroVarsGetForce);
    CPPUNIT_TEST_SUITE_END();
    public:
      void setUp() {
        FourCubeBasedTestFixture::setUp();
        // Set forces equal to site index
        size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
        for(size_t site(0); site < nFluidSites; ++site)
          latDat->GetSite(site).SetForce( GetMeAForce(site) );
      }

      void testHydroVarsGetForce() {
        size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
        typedef lb::kernels::GuoForcingLBGK<lb::lattices::D3Q15> t_Kernel;
        typedef t_Kernel::KHydroVars t_HydroVars;
        for(size_t i(1); i < nFluidSites; i <<= 1)
          CPPUNIT_ASSERT(helpers::is_zero(
            t_HydroVars(latDat->GetSite(i)).force - GetMeAForce(i)
          ));
      }

    protected:
      LatticeForceVector GetMeAForce(size_t site) {
        return LatticeForceVector(
            LatticeForce(site),
            LatticeForce(site) + 0.1,
            LatticeForce(site) + 0.3
        );
      }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION(GuoForcingTests);
}}

#endif
