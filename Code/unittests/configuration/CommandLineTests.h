
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
          argc = 5;
          argv[0] = "hemelb";
          argv[1] = "-in";
          argv[2] = configFile.c_str();
          argv[3] = "-ss";
          argv[4] = "1111";
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
