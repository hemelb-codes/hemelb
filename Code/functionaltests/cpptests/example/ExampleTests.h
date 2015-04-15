// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
          CPPUNIT_TEST_SUITE (ExampleTests);
          CPPUNIT_TEST (TestTrivialTest);CPPUNIT_TEST_SUITE_END();

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

      CPPUNIT_TEST_SUITE_REGISTRATION (ExampleTests);
    }
  }
}
#endif  //HEMELB_FUNCTIONALTESTS_CPPTESTS_EXAMPLE_EXAMPLETESTS_H
