#ifndef HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_GEOMETRYREADERTESTS_H
#include "geometry/LatticeData.h"
#include <cppunit/TestFixture.h>
#include "unittests/resources/Resource.h"
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
          class GlobalLatticeData : public LatticeData::GlobalLatticeData
          {
          };
      };
      class GeometryReaderTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(GeometryReaderTests);
          CPPUNIT_TEST(TestConstruct);
          CPPUNIT_TEST(TestRead);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            reader = new TestableLatticeData::GeometryReader(false);
            globalLattice= new TestableLatticeData::GlobalLatticeData();
            params=new lb::LbmParameters(1000, 0.1);
            simConfig = configuration::SimConfig::Load( Resource("config.xml").Path().c_str());
          }

          void tearDown()
          {
            delete reader;
          }

          void TestConstruct()
          {
            CPPUNIT_ASSERT(reader); // No smoke from construction
          }

          void TestRead()
          {
            reader->LoadAndDecompose(globalLattice, params, simConfig, timings);
          }

        private:
          TestableLatticeData::GeometryReader *reader;
          TestableLatticeData::GlobalLatticeData *globalLattice;
          configuration::SimConfig * simConfig;
          reporting::Timers timings;
          lb::LbmParameters *params;

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(GeometryReaderTests);
    }
  }
}
#endif // ONCE
