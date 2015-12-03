
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H
#define HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H

#include <cppunit/TestFixture.h>
#include "SimulationMaster.h"
#include "unittests/helpers/FolderTestFixture.h"
#include "unittests/helpers/LaddFail.h"

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
        CPPUNIT_TEST_SUITE( SimulationMasterTests);
        CPPUNIT_TEST( TestRun);CPPUNIT_TEST_SUITE_END();
      public:
        void setUp()
        {
          argc = 5;
          argv[0] = "hemelb";
          argv[1] = "-in";
          argv[2] = "four_cube.xml";
          argv[3] = "-ss";
          argv[4] = "1111";
          FolderTestFixture::setUp();
          CopyResourceToTempdir("four_cube.xml");
          CopyResourceToTempdir("four_cube.gmy");
          try {
            options = new hemelb::configuration::CommandLine(argc, argv);
            master = new SimulationMaster(*options, Comms());
          } catch (hemelb::io::xml::ChildError& e) {
            std::cout << e.what() << std::endl;
            throw;
          }
        }

        void tearDown()
        {
          FolderTestFixture::tearDown();
          delete master;
          delete options;
        }

        void TestRun()
        {
          // TODO: This test is fatal if run with LADDIOLET. See ticket #605.
          LADD_FAIL();
          master->RunSimulation();
          AssertPresent("results/report.txt");
          AssertPresent("results/report.xml");
          AssertPresent("results/Extracted/wholegeometryvelocityandstress.dat");
          AssertPresent("results/Extracted/centrelinepressure.dat");
          AssertPresent("results/Extracted/centrelineshearrate.dat");
          AssertPresent("results/Extracted/surfaceshearstress.dat");
          AssertPresent("results/Extracted/surfacetraction.dat");
        }

      private:
        int argc;
        hemelb::configuration::CommandLine *options;
        SimulationMaster *master;
        const char* argv[7];

    };

    CPPUNIT_TEST_SUITE_REGISTRATION( SimulationMasterTests);
  }
}

#endif /* HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H_ */
