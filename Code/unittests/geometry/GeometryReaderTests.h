// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#include "geometry/LatticeData.h"
#include <cppunit/TestFixture.h>
#include "lb/lattices/D3Q15.h"
#include "resources/Resource.h"
#include "unittests/FourCubeLatticeData.h"
#include "unittests/helpers/FolderTestFixture.h"
namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {
      using namespace hemelb::geometry;
      using namespace resources;

      class GeometryReaderTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE ( GeometryReaderTests);
          CPPUNIT_TEST ( TestRead);
          CPPUNIT_TEST ( TestSameAsFourCube);CPPUNIT_TEST_SUITE_END();

        public:

          GeometryReaderTests() :
            comms(*net::IOCommunicator::Instance()), timings(comms)
          {
          }

          void setUp()
          {
            reader = new GeometryReader(false,
                                        hemelb::lb::lattices::D3Q15::GetLatticeInfo(),
                                        timings, comms);
            lattice = NULL;
            fourCube = FourCubeLatticeData::Create();
            FolderTestFixture::setUp();
            CopyResourceToTempdir("four_cube.xml");
            CopyResourceToTempdir("four_cube.gmy");
            simConfig = NULL;
            simConfig = configuration::SimConfig::New("four_cube.xml");
          }

          void tearDown()
          {
            FolderTestFixture::tearDown();
            delete reader;
            delete lattice;
            delete fourCube;
            delete simConfig;
          }

          void TestRead()
          {
            reader->LoadAndDecompose(simConfig->GetDataFilePath());
          }

          void TestSameAsFourCube()
          {
            Geometry readResult = reader->LoadAndDecompose(simConfig->GetDataFilePath());

            for (site_t i = 1; i < 5; i++)
            {
              for (site_t j = 1; j < 5; j++)
              {
                bool isWallSite = (i == 1 || i == 4 || j == 1 || j == 4);

                for (site_t k = 1; k < 5; k++)
                {
                  //std::cout << i << "," << j << "," << k << " > " << std::setbase(8) << fourCube->GetSiteData(i*16+j*4+k) << " : " << globalLattice->GetSiteData(i,j,k) << std::endl;
                  util::Vector3D<site_t> location(i, j, k);
                  site_t siteIndex =
                      fourCube->GetGlobalNoncontiguousSiteIdFromGlobalCoords(location);

                  hemelb::geometry::SiteData siteData(readResult.Blocks[0].Sites[siteIndex]);
                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetSiteData(),
                                       siteData);
                  //                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetSiteData().GetOtherRawData(),
                  //                                       siteData.GetOtherRawData());
                  //
                  //                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetSiteData().GetWallIntersectionData(),
                  //                                       siteData.GetWallIntersectionData());

                  CPPUNIT_ASSERT_EQUAL(isWallSite,
                                       readResult.Blocks[0].Sites[siteIndex].wallNormalAvailable);

                  if (isWallSite)
                  {
                    /// @todo: #597 use CPPUNIT_ASSERT_EQUAL directly (having trouble with Vector3D templated over different types at the minute)
                    /// CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetWallNormal(), readResult.Blocks[0].Sites[siteIndex].wallNormal);
                    bool sameNormal =
                        (fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetWallNormal()
                            == readResult.Blocks[0].Sites[siteIndex].wallNormal);
                    CPPUNIT_ASSERT(sameNormal);
                  }
                }
              }
            }

          }

        private:
          net::IOCommunicator& comms;
          GeometryReader *reader;
          LatticeData* lattice;
          configuration::SimConfig * simConfig;
          reporting::Timers timings;
          hemelb::geometry::LatticeData *fourCube;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION ( GeometryReaderTests);
    }
  }
}
#endif // ONCE
