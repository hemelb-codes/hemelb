#ifndef HEMELB_FUNCTIONALTESTS_POISEUILLEFLOW_POISEUILLEFLOWTESTS_H
#define HEMELB_FUNCTIONALTESTS_POISEUILLEFLOW_POISEUILLEFLOWTESTS_H

#include "unittests/helpers/FolderTestFixture.h"
#include "SimulationMaster.h"

namespace hemelb
{
  namespace functionaltests
  {
    namespace poiseuilleflow
    {

      /**
       *  @todo This test suite inherits from a fixture defined in the unit tests directory and IMO
       *  creates a dependency that might hit us back at some point. Once both functionaltests and
       *  unittests have been moved into a tests folder, consider creating a helpers directory common
       *  to both, #249
       */
      class PoiseuilleFlowTests : public hemelb::unittests::helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE(PoiseuilleFlowTests);
          CPPUNIT_TEST(TestVelocityProfile);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            argc = 9;
            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = "poiseuille_flow_test.xml";
            argv[3] = "-i";
            argv[4] = "1";
            argv[5] = "-s";
            argv[6] = "1";
            argv[7] = "-ss";
            argv[8] = "1111";
            FolderTestFixture::setUp();
            CopyResourceToTempdir("poiseuille_flow_test.xml");
            CopyResourceToTempdir("poiseuille_flow_test.gmy");
            options = new hemelb::configuration::CommandLine(argc, argv);
            master = new SimulationMaster(*options);
          }

          void tearDown()
          {
            FolderTestFixture::tearDown();
            delete master;
            delete options;
          }

          void TestSimulationRunsAndGeneratesOutput()
          {
            master->RunSimulation();
            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
            AssertPresent("results/Extracted/40mmlinevelocity.dat");
            AssertPresent("results/Extracted/wholegeometryvelocity.dat");
            AssertPresent("results/Extracted/40mmlineshearstress.dat");
          }

          void TestVelocityProfile()
          {
          }

          void TestShearStressProfile()
          {
          }

        private:
          int argc;
          hemelb::configuration::CommandLine *options;
          SimulationMaster *master;
          const char* argv[9];

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(PoiseuilleFlowTests);
    }
  }
}
#endif  //HEMELB_FUNCTIONALTESTS_POISEUILLEFLOW_POISEUILLEFLOWTESTS_H
