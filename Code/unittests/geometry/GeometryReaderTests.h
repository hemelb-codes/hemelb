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
      // open the protected classes for testing
      class TestableLatticeData : public LatticeData
      {
        public:
          class GeometryReader : public hemelb::geometry::GeometryReader
          {
            public:
              GeometryReader(const bool reserveSteeringCore,
                             const hemelb::lb::lattices::LatticeInfo& latticeInfo,
                             reporting::Timers &timings) :
                hemelb::geometry::GeometryReader(reserveSteeringCore, latticeInfo, timings)
              {
              }
          };
      };
      class GeometryReaderTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE( GeometryReaderTests);
          CPPUNIT_TEST( TestRead);
          CPPUNIT_TEST( TestSameAsFourCube);CPPUNIT_TEST_SUITE_END();

        public:

          GeometryReaderTests() :
            timings()
          {
          }

          void setUp()
          {
            reader = new TestableLatticeData::GeometryReader(false,
                                                             hemelb::lb::lattices::D3Q15::GetLatticeInfo(),
                                                             timings);
            lattice = NULL;
            bool dummy;
            topology::NetworkTopology::Instance()->Init(0, NULL, &dummy);
            fourCube = FourCubeLatticeData::Create();
            FolderTestFixture::setUp();
            CopyResourceToTempdir("four_cube.xml");
            CopyResourceToTempdir("four_cube.gmy");
            simConfig = configuration::SimConfig::Load("four_cube.xml");
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

            site_t siteIndex = 0;
            for (site_t i = 0; i < 4; i++)
            {
              for (site_t j = 0; j < 4; j++)
              {
                for (site_t k = 0; k < 4; k++)
                {
                  //std::cout << i << "," << j << "," << k << " > " << std::setbase(8) << fourCube->GetSiteData(i*16+j*4+k) << " : " << globalLattice->GetSiteData(i,j,k) << std::endl;
                  util::Vector3D<site_t> location(i, j, k);

                  hemelb::geometry::SiteData siteData(readResult.Blocks[0].Sites[siteIndex]);

                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetSiteData().GetOtherRawData(),
                                       siteData.GetOtherRawData());

                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(location)).GetSiteData().GetIntersectionData(),
                                       siteData.GetIntersectionData());

                  siteIndex++;
                }
              }
            }

          }

        private:
          TestableLatticeData::GeometryReader *reader;
          TestableLatticeData* lattice;
          configuration::SimConfig * simConfig;
          reporting::Timers timings;
          hemelb::geometry::LatticeData *fourCube;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION( GeometryReaderTests);
    }
  }
}
#endif // ONCE
