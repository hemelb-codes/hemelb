
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
              CPPUNIT_ASSERT(!requirements->IsRequired(terms::Distribution));
            }

            void TestSpecifyRequirement()
            {
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(requirements->IsRequired(terms::Distribution));
            }

            void TestAny()
            {
              CPPUNIT_ASSERT(!requirements->RequiresAny());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(requirements->RequiresAny());
            }

            void TestAnyNonFieldDependent()
            {
              CPPUNIT_ASSERT(!requirements->RequiresAnyNonFieldDependent());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(!requirements->RequiresAnyNonFieldDependent());
              requirements->Require(terms::SiteData);
              CPPUNIT_ASSERT(requirements->RequiresAnyNonFieldDependent());
            }

            void TestAnyFieldDependent()
            {
              CPPUNIT_ASSERT(!requirements->RequiresAnyFieldDependent());
              requirements->Require(terms::SiteData);
              CPPUNIT_ASSERT(!requirements->RequiresAnyFieldDependent());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(requirements->RequiresAnyFieldDependent());
            }

            void TestAnyMacroscopic()
            {
              CPPUNIT_ASSERT(!requirements->RequiresAnyMacroscopic());
              requirements->Require(terms::SiteData);
              CPPUNIT_ASSERT(!requirements->RequiresAnyMacroscopic());
              requirements->Require(terms::Distribution);
              CPPUNIT_ASSERT(!requirements->RequiresAnyMacroscopic());
              requirements->Require(terms::Velocity);
              CPPUNIT_ASSERT(requirements->RequiresAnyMacroscopic());
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
