// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_WALLNODEDNCTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_WALLNODEDNCTESTS_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "unittests/FourCubeLatticeData.h"
#include "unittests/helpers/HasCommsTestFixture.h"
#include "redblood/WallCellPairIterator.h"
#include "lb/lattices/D3Q15.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class WallNodeDnCTests : public helpers::HasCommsTestFixture
      {
          CPPUNIT_TEST_SUITE (WallNodeDnCTests);
          CPPUNIT_TEST (testWallNodeDnC);CPPUNIT_TEST_SUITE_END();

          LatticeDistance const cutoff = 3.0;
          LatticeDistance const halo = 0.5;
          typedef lb::lattices::D3Q15 Lattice;

        public:
          void setUp()
          {
            latticeData.reset(FourCubeLatticeData::Create(Comms(), 27 + 2));
            for (site_t i(0); i < latticeData->GetLocalFluidSiteCount(); ++i)
            {
              auto const site = latticeData->GetSite(i);
              if (not site.IsWall())
              {
                continue;
              }
              for (Direction d(0); d < Lattice::NUMVECTORS; ++d)
              {
                if (site.HasWall(d))
                {
                  latticeData->SetBoundaryDistance(i, d, 0.5);
                }
              }
            }
          }
          void testWallNodeDnC()
          {
            using namespace hemelb::redblood;
            auto const dnc = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);

            // Checking the middle of the world first: No wall nodes
            auto const center = dnc.equal_range(LatticePosition(13, 13, 13));
            CPPUNIT_ASSERT(center.first == center.second);

            // Checking that upper side has 3*3*5 wall sites,
            // 3*3 because that's the size of a DnC box
            // 5 because there are five possible directions in D3Q15 pointing towards the wall
            auto upper = dnc.equal_range(LatticePosition(0, 3.5 * 3, 3.5 * 3));
            CPPUNIT_ASSERT_EQUAL(45l, std::distance(upper.first, upper.second));
            for (; upper.first != upper.second; ++upper.first)
            {
              CPPUNIT_ASSERT(upper.first->second.nearBorder bitand size_t(Borders::CENTER));
              CPPUNIT_ASSERT_EQUAL( (upper.first->second.nearBorder bitand size_t(Borders::BOTTOM))
                                       != 0,
                                   upper.first->second.node.x < halo);
            }
          }

        private:
          std::unique_ptr<unittests::FourCubeLatticeData> latticeData;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (WallNodeDnCTests);
    }
  }
}

#endif  // ONCE
