
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#include "geometry/LatticeData.h"
#include <cppunit/TestFixture.h>
#include "lb/lattices/D3Q15.h"
#include "resources/Resource.h"
#include "unittests/FourCubeLatticeData.h"
#include "unittests/helpers/FolderTestFixture.h"
#include "unittests/helpers/LaddFail.h"

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

          GeometryReaderTests()
          {
          }

          void setUp()
          {
            FolderTestFixture::setUp();
            timings = new reporting::Timers(Comms());
            reader = new GeometryReader(hemelb::lb::lattices::D3Q15::GetLatticeInfo(),
                                        *timings, Comms());
            lattice = NULL;
            fourCube = FourCubeLatticeData::Create(Comms());
            CopyResourceToTempdir("four_cube.xml");
            CopyResourceToTempdir("four_cube.gmy");
            simConfig = NULL;
            simConfig = configuration::SimConfig::New("four_cube.xml");
          }

          void tearDown()
          {
            FolderTestFixture::tearDown();
            delete timings;
            delete reader;
            delete lattice;
            delete fourCube;
            delete simConfig;
          }

          void TestRead()
          {
            LADD_FAIL();
            reader->LoadAndDecompose(simConfig->GetDataFilePath());
          }

          void TestSameAsFourCube()
          {
            LADD_FAIL();
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
          GeometryReader *reader;
          LatticeData* lattice;
          configuration::SimConfig * simConfig;
          reporting::Timers* timings;
          hemelb::geometry::LatticeData *fourCube;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION ( GeometryReaderTests);
    }
  }
}
#endif // ONCE
