// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H
#define HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H

#include <cppunit/TestFixture.h>
#include "SimulationMaster.h"
#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    /**
     * Class to test the simulation master.
     */
    using namespace helpers;
    class SimulationMasterTests : public FolderTestFixture
    {
        CPPUNIT_TEST_SUITE(SimulationMasterTests);
        CPPUNIT_TEST(TestRun);
        CPPUNIT_TEST_SUITE_END();
      public:
        void setUp()
        {
          argc = 9;
          argv[0] = "hemelb";
          argv[2] = "four_cube.xml";
          argv[1] = "-in";
          argv[3] = "-i";
          argv[4] = "1";
          argv[5] = "-s";
          argv[6] = "1";
          argv[7] = "-ss";
          argv[8] = "1111";
          FolderTestFixture::setUp();
          CopyResourceToTempdir("four_cube.xml");
          CopyResourceToTempdir("four_cube.gmy");
          options = new hemelb::configuration::CommandLine(argc, argv);
          master = new SimulationMaster(*options);
        }

        void tearDown()
        {
          FolderTestFixture::tearDown();
          delete master;
          delete options;
        }

        void TestRun()
        {
          master->RunSimulation();
          AssertPresent("results/report.txt");
          AssertPresent("results/report.xml");
          AssertPresent("results/Extracted/wholegeometryvelocityandstress.dat");
          AssertPresent("results/Extracted/centrelinepressure.dat");
          AssertPresent("results/Extracted/centrelineshearrate.dat");
          AssertPresent("results/Extracted/surfaceshearstress.dat");
        }

      private:
        int argc;
        hemelb::configuration::CommandLine *options;
        SimulationMaster *master;
        const char* argv[9];

    };

    CPPUNIT_TEST_SUITE_REGISTRATION(SimulationMasterTests);
  }
}

#endif /* HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H_ */
