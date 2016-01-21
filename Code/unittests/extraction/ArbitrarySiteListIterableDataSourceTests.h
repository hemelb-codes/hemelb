
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCETESTS_H
#define HEMELB_UNITTESTS_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCETESTS_H

#include <string>
#include <cstdio>

#include <cppunit/TestFixture.h>
#include "extraction/ArbitrarySiteListIterableDataSource.h"
#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/MockNetHelper.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace unittests
  {
    namespace extraction
    {
      class ArbitrarySiteListIterableDataSourceTests : public FourCubeBasedTestFixture,
                                                       public MockNetHelper
      {
          CPPUNIT_TEST_SUITE (ArbitrarySiteListIterableDataSourceTests);
          CPPUNIT_TEST (TestManagerAddition);
          CPPUNIT_TEST (TestVelocityExtraction);CPPUNIT_TEST_SUITE_END();

        public:

          void setUp()
          {
            FourCubeBasedTestFixture::setUp();
            data = &latDat->GetNeighbouringData();
            MockNetHelper::setUp(1, 0);
            manager = new hemelb::geometry::neighbouring::NeighbouringDataManager(*latDat, *data, *netMock);
            TestableIterableDataSource = new hemelb::extraction::ArbitrarySiteListIterableDataSource();
          }

          void tearDown()
          {
            FourCubeBasedTestFixture::tearDown();
            MockNetHelper::tearDown();
            delete manager;
            delete TestableIterableDataSource;
          }

          void TestManagerAddition()
          {
            TestableIterableDataSource->SetManager(manager);
          }

          void TestVelocityExtraction()
          {
            TestableIterableDataSource->SetManager(manager);
            util::Vector3D<float> normal = util::Vector3D<float>(1.0,1.0,1.0);
            float velocity = TestableIterableDataSource->GetVelocityRelativeToNormal(manager, normal);
            std::cout << "Obtained velocity is: " << velocity << std::endl;
          }

        protected:
          hemelb::extraction::ArbitrarySiteListIterableDataSource *TestableIterableDataSource;
          hemelb::geometry::neighbouring::NeighbouringDataManager *manager;
          hemelb::geometry::neighbouring::NeighbouringLatticeData *data;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (ArbitrarySiteListIterableDataSourceTests);
    }
  }
}

#endif /* ARBITRARYSITELISTITERABLEDATASOURCETESTS_H */
