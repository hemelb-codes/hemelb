// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_CONFIGURATION_COMMANDLINETESTS_H
#define HEMELB_UNITTESTS_CONFIGURATION_COMMANDLINETESTS_H

#include <cppunit/TestFixture.h>
#include "configuration/CommandLine.h"
#include "resources/Resource.h"
#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    /**
     * Class to test the simulation master.
     */
    using namespace hemelb::configuration;
    using namespace resources;
    using namespace helpers;
    class CommandLineTests : public FolderTestFixture
    {
        CPPUNIT_TEST_SUITE(CommandLineTests);
        CPPUNIT_TEST(TestConstruct);
        CPPUNIT_TEST_SUITE_END();
      public:
        void setUp()
        {
          configFile = Resource("four_cube.xml").Path();
          argc = 7;
          argv[0] = "hemelb";
          argv[2] = configFile.c_str();
          argv[1] = "-in";
          argv[3] = "-i";
          argv[4] = "1";
          argv[5] = "-ss";
          argv[6] = "1111";
          FolderTestFixture::setUp();
          options = new hemelb::configuration::CommandLine(argc, argv);
        }

        void tearDown()
        {
          FolderTestFixture::tearDown();
          delete options;
        }

        void TestConstruct()
        {
          CPPUNIT_ASSERT(options);
        }

      private:
        int argc;
        std::string configFile;
        hemelb::configuration::CommandLine *options;
        const char* argv[7];

    };
    CPPUNIT_TEST_SUITE_REGISTRATION(CommandLineTests);
  }
}

#endif // ONCE
