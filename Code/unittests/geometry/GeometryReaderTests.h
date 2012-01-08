#ifndef HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#include "geometry/LatticeData.h"
#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "unittests/FourCubeLatticeData.h"
namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {
      using namespace hemelb::geometry;
      using namespace resources;

      class GeometryReaderTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE( GeometryReaderTests);
          CPPUNIT_TEST( TestRead);
          CPPUNIT_TEST( TestSameAsFourCube);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            reader = new hemelb::geometry::GeometryReader(false, readResult);
            params = new lb::LbmParameters(1000, 0.1);
            bool dummy;
            topology::NetworkTopology::Instance()->Init(0, NULL, &dummy);
            fourCube = FourCubeLatticeData::Create();
            simConfig = configuration::SimConfig::Load(Resource("four_cube.xml").Path().c_str());
          }

          void tearDown()
          {
            delete reader;
            delete params;
            delete fourCube;
            delete simConfig;
          }

          void TestRead()
          {
            reader->LoadAndDecompose(simConfig->DataFilePath, params, timings);
          }

          void TestSameAsFourCube()
          {
            reader->LoadAndDecompose(simConfig->DataFilePath, params, timings);
            site_t siteIndex = 0;
            for (site_t i = 0; i < 4; i++)
            {
              for (site_t j = 0; j < 4; j++)
              {
                for (site_t k = 0; k < 4; k++)
                {
                  //std::cout << i << "," << j << "," << k << " > " << std::setbase(8) << fourCube->GetSiteData(i*16+j*4+k) << " : " << globalLattice->GetSiteData(i,j,k) << std::endl;
                  CPPUNIT_ASSERT_EQUAL(fourCube->GetSiteData(siteIndex).GetRawValue(),
                                       readResult.Blocks[0].Sites[siteIndex].siteData.GetRawValue());

                  siteIndex++;
                }
              }
            }

          }

        private:
          GeometryReader *reader;
          GeometryReadResult readResult;
          configuration::SimConfig * simConfig;
          reporting::Timers timings;
          lb::LbmParameters *params;
          hemelb::geometry::LatticeData *fourCube;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION( GeometryReaderTests);
    }
  }
}
#endif // ONCE
