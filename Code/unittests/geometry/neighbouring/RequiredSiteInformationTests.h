#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_REQUIREDSITEINFORMATIONTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_REQUIREDSITEINFORMATIONTESTS_H

#include "geometry/neighbouring/RequiredSiteInformation.h"
namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {
      namespace neighbouring
      {
        using namespace hemelb::geometry::neighbouring;
        class RequiredSiteInformationTests : public CppUnit::TestFixture
        {
            CPPUNIT_TEST_SUITE (RequiredSiteInformationTests);
            CPPUNIT_TEST (TestConstruct);
            CPPUNIT_TEST (TestSpecifyRequirement);
            CPPUNIT_TEST (TestAny);
            CPPUNIT_TEST (TestAnyNonFieldDependent);
            CPPUNIT_TEST (TestAnyFieldDependent);
            CPPUNIT_TEST (TestAnyMacroscopic);

            CPPUNIT_TEST_SUITE_END();

          public:
            RequiredSiteInformationTests()
            {
            }

            void setUp()
            {
              requirements = new RequiredSiteInformation();
            }

            void tearDown()
            {
              delete requirements;
            }

            void TestConstruct()
            {
              CPPUNIT_ASSERT(!requirements->Required(terms::Distribution));
            }

            void TestSpecifyRequirement()
            {
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(requirements->Required(terms::Distribution));
            }

            void TestAny()
            {
              CPPUNIT_ASSERT(!requirements->Any());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(requirements->Any());
            }

            void TestAnyNonFieldDependent()
            {
              CPPUNIT_ASSERT(!requirements->AnyNonFieldDependent());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(!requirements->AnyNonFieldDependent());
              requirements->Require(terms::SiteData);
              CPPUNIT_ASSERT(requirements->AnyNonFieldDependent());
            }

            void TestAnyFieldDependent()
            {
              CPPUNIT_ASSERT(!requirements->AnyFieldDependent());
              requirements->Require(terms::SiteData);
              CPPUNIT_ASSERT(!requirements->AnyFieldDependent());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(requirements->AnyFieldDependent());
            }

            void TestAnyMacroscopic()
            {
              CPPUNIT_ASSERT(!requirements->AnyMacroscopic());
              requirements->Require(terms::SiteData);
              CPPUNIT_ASSERT(!requirements->AnyMacroscopic());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(!requirements->AnyMacroscopic());
              requirements->Require(terms::Velocity);
              CPPUNIT_ASSERT(requirements->AnyMacroscopic());
            }

          private:
            RequiredSiteInformation *requirements;
        };
        // CPPUNIT USES LINENUMBER TO REGISTER MACRO
        // EXTRA LINE
        // EXTRA LINE
        // EXTRA LINE
        CPPUNIT_TEST_SUITE_REGISTRATION (RequiredSiteInformationTests);
      }
    }
  }
}

#endif
