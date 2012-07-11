/*
 * ArbitrarySiteListIterableDataSourceTest.h
 *
 *  Created on: 3 Jul 2012
 *      Author: derek
 */

#ifndef HEMELB_UNITTESTS_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCETEST_H_
#define HEMELB_UNITTESTS_EXTRACTION_ARBITRARYSITELISTITERABLEDATASOURCETEST_H_

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

#endif /* ARBITRARYSITELISTITERABLEDATASOURCETEST_H_ */
