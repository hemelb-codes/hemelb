#ifndef HEMELB_FUNCTIONALTESTS_POISEUILLEFLOW_POISEUILLEFLOWTESTS_H
#define HEMELB_FUNCTIONALTESTS_POISEUILLEFLOW_POISEUILLEFLOWTESTS_H

#include <cppunit/TestFixture.h>
//#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace functionaltests
  {
    namespace poiseuilleflow
    {

      class PoiseuilleFlowTests : public CppUnit::TestFixture //FolderTestFixture
      {
          CPPUNIT_TEST_SUITE(PoiseuilleFlowTests);
          CPPUNIT_TEST(TestVelocityProfile);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
          }

          void tearDown()
          {
          }

          void TestVelocityProfile()
          {
            CPPUNIT_ASSERT(true);
          }

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(PoiseuilleFlowTests);
    }
  }
}
#endif  //HEMELB_FUNCTIONALTESTS_POISEUILLEFLOW_POISEUILLEFLOWTESTS_H
