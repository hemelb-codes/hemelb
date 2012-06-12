#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITETESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGSITETESTS_H
namespace unittests
{
  namespace hemelb
  {
    namespace geometry
    {
      namespace neighbouring
      {
        using namespace hemelb::geometry::neighbouring;
        class NeighbouringSiteTests : public CppUnit::TestFixture
        {
            CPPUNIT_TEST_SUITE (NeighbouringSiteTests);
            CPPUNIT_TEST (TestConstruct);

            CPPUNIT_TEST_SUITE_END();

          public:
            NeighbouringSiteTests()
            {
            }

            void setUp()
            {

            }

            void tearDown()
            {
            }

            void TestConstruct()
            {
              // PASS -- just verify setUp and tearDown
            }

          private:

        };

        // CPPUNIT USES LINENUMBER GRR
        // EXTRA LINE
        CPPUNIT_TEST_SUITE_REGISTRATION (NeighbouringSiteTests);
      }
    }
  }
}

#endif
