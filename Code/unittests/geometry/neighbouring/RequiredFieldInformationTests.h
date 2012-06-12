#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_REQUIREDFIELDINFORMATIONTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_REQUIREDFIELDINFORMATIONTESTS_H
namespace unittests
{
  namespace hemelb
  {
    namespace geometry
    {
      namespace neighbouring
      {
        using namespace hemelb::geometry::neighbouring;
        class RequiredFieldInformationTests : public CppUnit::TestFixture
        {
            CPPUNIT_TEST_SUITE (RequiredFieldInformationTests);
            CPPUNIT_TEST (TestConstruct);

            CPPUNIT_TEST_SUITE_END();

          public:
            RequiredFieldInformationTests()
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
        // CPPUNIT USES LINENUMBER TO REGISTER MACRO
        // EXTRA LINE
        // EXTRA LINE
        // EXTRA LINE
        CPPUNIT_TEST_SUITE_REGISTRATION (RequiredFieldInformationTests);
      }
    }
  }
}

#endif
