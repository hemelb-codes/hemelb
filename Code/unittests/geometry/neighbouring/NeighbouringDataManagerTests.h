#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
namespace unittests
{
  namespace hemelb
  {
    namespace geometry
    {
      namespace neighbouring
      {
        using namespace hemelb::geometry::neighbouring;
        class NeighbouringDataManagerTests : public CppUnit::TestFixture
        {
            CPPUNIT_TEST_SUITE (NeighbouringDataManagerTests);
            CPPUNIT_TEST (TestConstruct);

            CPPUNIT_TEST_SUITE_END();

          public:
            NeighbouringDataManagerTests()
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
        CPPUNIT_TEST_SUITE_REGISTRATION (NeighbouringDataManagerTests);
      }
    }
  }
}

#endif
