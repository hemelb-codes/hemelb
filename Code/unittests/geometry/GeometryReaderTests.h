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
          class GeometryReader : public LatticeData::GeometryReader
          {
            public:
              GeometryReader(const bool reserveSteeringCore) :
                LatticeData::GeometryReader(reserveSteeringCore)
              {
              }
          };

          typedef LatticeData::GlobalLatticeData TestableGlobalLatticeData;
      };
      class GeometryReaderTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE( GeometryReaderTests);
          CPPUNIT_TEST( TestRead);
          CPPUNIT_TEST( TestSameAsFourCube);CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            reader = new TestableLatticeData::GeometryReader(false);
            globalLattice = NULL;
            bool dummy;
            topology::NetworkTopology::Instance()->Init(0, NULL, &dummy);
            fourCube = FourCubeLatticeData::CubeGlobalLatticeData::Create();
            FolderTestFixture::setUp();
            CopyResourceToTempdir("four_cube.xml");
            CopyResourceToTempdir("four_cube.dat");
            simConfig = configuration::SimConfig::Load("four_cube.xml");
          }

          void tearDown()
          {
            FolderTestFixture::tearDown();
            delete reader;
            delete globalLattice;
            delete fourCube;
            delete simConfig;
          }

          void TestRead()
          {
            globalLattice = reader->LoadAndDecompose(simConfig->DataFilePath, timings);
          }

          void TestSameAsFourCube()
          {
            globalLattice = reader->LoadAndDecompose(simConfig->DataFilePath, timings);
            for (site_t i = 0; i < 4; i++)
            {
              for (site_t j = 0; j < 4; j++)
              {
                for (site_t k = 0; k < 4; k++)
                {
                  //std::cout << i << "," << j << "," << k << " > " << std::setbase(8) << fourCube->GetSiteData(i*16+j*4+k) << " : " << globalLattice->GetSiteData(i,j,k) << std::endl;
                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSiteData(i, j, k),
                                       globalLattice->GetSiteData(i, j, k));
                }
              }
            }

          }

        private:
          TestableLatticeData::GeometryReader *reader;
          TestableLatticeData::TestableGlobalLatticeData* globalLattice;
          configuration::SimConfig * simConfig;
          reporting::Timers timings;
          FourCubeLatticeData::CubeGlobalLatticeData *fourCube;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION( GeometryReaderTests);
    }
  }
}
#endif // ONCE
