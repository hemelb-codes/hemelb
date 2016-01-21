
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGLATTICEDATATESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGLATTICEDATATESTS_H
#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "geometry/neighbouring/NeighbouringSite.h"
#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/CppUnitCompareSiteData.h"
namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {
      namespace neighbouring
      {
        using namespace hemelb::geometry::neighbouring;
        class NeighbouringLatticeDataTests : public FourCubeBasedTestFixture
        {
            CPPUNIT_TEST_SUITE (NeighbouringLatticeDataTests);
            CPPUNIT_TEST (TestConstruct);
            CPPUNIT_TEST (TestInsertAndRetrieveSiteData);
            CPPUNIT_TEST (TestInsertAndRetrieveDistance);
            CPPUNIT_TEST (TestInsertAndRetrieveNormal);
            CPPUNIT_TEST (TestInsertAndRetrieveDistributions);
            CPPUNIT_TEST (TestNeighbouringSite);

            CPPUNIT_TEST_SUITE_END();

          public:
            NeighbouringLatticeDataTests() :
                FourCubeBasedTestFixture(),data(NULL), exampleSite(NULL), dummyId(54)
            {
            }

            void setUp()
            {
              FourCubeBasedTestFixture::setUp();
              data = &latDat->GetNeighbouringData();
              exampleSite = new Site<LatticeData>(latDat->GetSite(24));
            }

            void tearDown()
            {
              delete exampleSite;
              FourCubeBasedTestFixture::tearDown();
            }

            void TestConstruct()
            {
              // PASS -- just verify setUp and tearDown
            }

            void TestInsertAndRetrieveSiteData()
            {
              data->GetSiteData(dummyId) = exampleSite->GetSiteData();
              CPPUNIT_ASSERT_EQUAL(exampleSite->GetSiteData(), data->GetSiteData(dummyId));
            }

            void TestInsertAndRetrieveDistance()
            {
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
              {
                data->GetCutDistances(dummyId)[direction]=exampleSite->GetWallDistance < lb::lattices::D3Q15 > (direction + 1);
              }

              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
              {
                CPPUNIT_ASSERT_EQUAL(exampleSite->GetWallDistance < lb::lattices::D3Q15 > (direction + 1),
                                     data->GetCutDistance<lb::lattices::D3Q15>(dummyId, direction + 1));
              }
            }

            void TestInsertAndRetrieveNormal()
            {
              data->GetNormalToWall(dummyId) = exampleSite->GetWallNormal();
              CPPUNIT_ASSERT_EQUAL(exampleSite->GetWallNormal(), data->GetNormalToWall(dummyId));
            }

            void TestInsertAndRetrieveDistributions()
            {
              std::vector<distribn_t> distribution;
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
              {
                distribution.push_back(exampleSite->GetFOld<lb::lattices::D3Q15>()[direction]);
              }

              data->GetDistribution(dummyId) = distribution;

              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
              {
                CPPUNIT_ASSERT_EQUAL(exampleSite->GetFOld<lb::lattices::D3Q15>()[direction],
                                     data->GetFOld(dummyId * lb::lattices::D3Q15::NUMVECTORS)[direction]);
              }
            }

            void TestNeighbouringSite()
            {

              std::vector<distribn_t> distances;
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
              {
                distances.push_back(exampleSite->GetWallDistance < lb::lattices::D3Q15 > (direction + 1));
              }

              std::vector<distribn_t> distribution;
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
              {
                distribution.push_back(exampleSite->GetFOld<lb::lattices::D3Q15>()[direction]);
              }
              data->SaveSite(dummyId,
                             distribution,
                             distances,
                             exampleSite->GetWallNormal(),
                             exampleSite->GetSiteData());

              NeighbouringSite neighbouringSite = data->GetSite(dummyId);

              CPPUNIT_ASSERT_EQUAL(exampleSite->GetSiteData(), neighbouringSite.GetSiteData());
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
              {
                CPPUNIT_ASSERT_EQUAL(exampleSite->GetWallDistance < lb::lattices::D3Q15 > (direction + 1),
                                     neighbouringSite.GetWallDistance<lb::lattices::D3Q15>(direction + 1));
              }
              CPPUNIT_ASSERT_EQUAL(exampleSite->GetWallNormal(), neighbouringSite.GetWallNormal());
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
              {
                CPPUNIT_ASSERT_EQUAL(exampleSite->GetFOld<lb::lattices::D3Q15>()[direction],
                                     neighbouringSite.GetFOld<lb::lattices::D3Q15>()[direction]);
              }
            }

          private:
            NeighbouringLatticeData *data;
            Site<LatticeData> *exampleSite;
            site_t dummyId;
        };
        // CPPUNIT USES LINENUMBER TO REGISTER MACRO
        // EXTRA LINE
        // EXTRA LINE
        CPPUNIT_TEST_SUITE_REGISTRATION (NeighbouringLatticeDataTests);
      }
    }
  }
}

#endif
