
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_FUNCTIONALTESTS_CPPTESTS_EXAMPLE_EXAMPLETESTS_H
#define HEMELB_FUNCTIONALTESTS_CPPTESTS_EXAMPLE_EXAMPLETESTS_H

#include <cppunit/TestFixture.h>

namespace hemelb
{
  namespace functionaltests
  {
    namespace example
    {

      class ExampleTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(ExampleTests);
          CPPUNIT_TEST(TestTrivialTest);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
          }

          void tearDown()
          {
          }

          void TestTrivialTest()
          {
            CPPUNIT_ASSERT(true);
          }

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(ExampleTests);
    }
  }
}
#endif  //HEMELB_FUNCTIONALTESTS_CPPTESTS_EXAMPLE_EXAMPLETESTS_H
