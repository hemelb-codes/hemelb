#ifndef HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#include "geometry/LatticeData.h"
#include <cppunit/TestFixture.h>
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
                             hemelb::geometry::GeometryReadResult& readResult) :
                  hemelb::geometry::GeometryReader(reserveSteeringCore, readResult)
              {
              }
          };
      };
      class GeometryReaderTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE(GeometryReaderTests);
          CPPUNIT_TEST(TestRead);
          CPPUNIT_TEST(TestSameAsFourCube);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            readResult = new GeometryReadResult();
            reader = new TestableLatticeData::GeometryReader(false, *readResult);
            lattice = NULL;
            bool dummy;
            topology::NetworkTopology::Instance()->Init(0, NULL, &dummy);
            fourCube = FourCubeLatticeData::Create();
            FolderTestFixture::setUp();
            CopyResourceToTempdir("four_cube.xml");
            CopyResourceToTempdir("four_cube.dat");
            simConfig = configuration::SimConfig::Load("four_cube.xml");
          }

          void tearDown()
          {
            FolderTestFixture::tearDown();
            delete reader;
            delete lattice;
            delete fourCube;
            delete simConfig;
            delete readResult;
          }

          void TestRead()
          {
            reader->LoadAndDecompose(simConfig->DataFilePath, timings);
          }

          void TestSameAsFourCube()
          {
            reader->LoadAndDecompose(simConfig->DataFilePath, timings);

            site_t siteIndex = 0;
            for (site_t i = 0; i < 4; i++)
            {
              for (site_t j = 0; j < 4; j++)
              {
                for (site_t k = 0; k < 4; k++)
                {
                  //std::cout << i << "," << j << "," << k << " > " << std::setbase(8) << fourCube->GetSiteData(i*16+j*4+k) << " : " << globalLattice->GetSiteData(i,j,k) << std::endl;
                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSite(fourCube->GetContiguousSiteId(i, j, k)).GetSiteData().GetRawValue(),
                                       readResult->Blocks[0].Sites[siteIndex].siteData.GetRawValue());

                  siteIndex++;
                }
              }
            }

          }

        private:
          TestableLatticeData::GeometryReader *reader;
          GeometryReadResult* readResult;
          TestableLatticeData* lattice;
          configuration::SimConfig * simConfig;
          reporting::Timers timings;
          hemelb::geometry::LatticeData *fourCube;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(GeometryReaderTests);
    }
  }
}
#endif // ONCE
