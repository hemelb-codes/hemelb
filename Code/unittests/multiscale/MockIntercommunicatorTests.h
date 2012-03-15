#ifndef HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
#define HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
#include "unittests/multiscale/MockIntercommunicand.h"
#include "resources/Resource.h"
#include "unittests/multiscale/MockIntercommunicator.h"
#include "multiscale/MultiscaleSimulationMaster.h"

#include <iostream>

namespace hemelb
{
  namespace unittests
  {
    namespace multiscale
    {
      template<class IntercommuniatorImplementation> class MockHemeLB
      {
        public:
          MockHemeLB(double aspace_resolution, double atime_resolution, std::map<std::string, double> & buffer) :
              inlet(), outlet(), inoutlettype("inoutlet"), intercomms(buffer), time_resolution(atime_resolution), space_resolution(aspace_resolution), current_time(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inoutlettype.template RegisterSharedValue<double>("pressure");
            inoutlettype.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inoutlettype, inlet, "boundary1");
            intercomms.RegisterIntercommunicand(inoutlettype, outlet, "boundary2");
          }
          MockIntercommunicand inlet;
          MockIntercommunicand outlet;
          typename IntercommuniatorImplementation::IntercommunicandTypeT inoutlettype;
          IntercommuniatorImplementation intercomms;
          double time_resolution;
          double space_resolution;
          double current_time;

          void DoLB()
          {
            double resistance = 10.0;
            double velocity = (inlet.GetDensity() - outlet.GetDensity()) / resistance;
            inlet.SetVelocity(velocity);
            outlet.SetVelocity(velocity);
          }

          void Simulate()
          {
            if (!intercomms.ShouldAdvance())
              return;
            intercomms.GetFromMultiscale();
            DoLB();
            current_time += time_resolution;
            intercomms.AdvanceTime(current_time);
            intercomms.SendToMultiscale();
          }

      };
      template<class IntercommuniatorImplementation> class Mock0DModel
      {
        public:
          Mock0DModel(double aspace_resolution, double atime_resolution, std::map<std::string, double> & buffer) :
              inlet(), outlet(), inoutlettype("inoutlet"), intercomms(buffer), time_resolution(atime_resolution), space_resolution(aspace_resolution), current_time(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inoutlettype.template RegisterSharedValue<double>("pressure");
            inoutlettype.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inoutlettype, inlet, "boundary2");
            intercomms.RegisterIntercommunicand(inoutlettype, outlet, "boundary1");
          }
          MockIntercommunicand inlet;
          MockIntercommunicand outlet;
          typename IntercommuniatorImplementation::IntercommunicandTypeT inoutlettype;
          IntercommuniatorImplementation intercomms;
          double time_resolution;
          double space_resolution;
          double current_time;

          void Do1D()
          {
            double capacitance = 10.0;
            double deltap = -1.0 * inlet.GetVelocity() / capacitance;

            outlet.SetDensity(outlet.GetDensity() + deltap);
          }

          void Simulate()
          {
            if (!intercomms.ShouldAdvance())
              return;
            intercomms.GetFromMultiscale();
            Do1D();
            current_time += time_resolution;
            intercomms.AdvanceTime(current_time);
            intercomms.SendToMultiscale();
          }

      };

      class MockIntercommunicatorTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE(MockIntercommunicatorTests);
          CPPUNIT_TEST(TestCRRun);
          CPPUNIT_TEST(TestCHemeRun);CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            pbuffer = new std::map<std::string, double>();
            std::map<std::string, double> &buffer = *pbuffer;
            buffer["boundary1_pressure"] = 81.0;
            buffer["boundary2_pressure"] = 79.0;
            buffer["boundary1_velocity"] = 0.0;
            buffer["boundary2_velocity"] = 0.0;
            buffer["shared_time"] = 0.0;
            mockheme = new MockHemeLB<MockIntercommunicator>(25.0, 0.2, buffer);
            zerod = new Mock0DModel<MockIntercommunicator>(10.0, 0.5, buffer);

          }
          void tearDown()
          {
            delete zerod;
            delete mockheme;
          }
        private:
          void TestCRRun()
          {
            // In this example, we imagine we are coupling inlets in two different instances of HemeLB.
            // Each has it's own instance of an intercommunicator - these would be running in different processes.
            // Our example intercommunicator has been hard coded with how to share the values.
            // In principle, though, this would be defined by whatever configuration file format the implementing intercommunicator wants to use.
            for (int steps = 0; steps < 100; steps++)
            {
              mockheme->Simulate();
              zerod->Simulate();
              //std::cout << "Heme: I: p:" << heme->inlet.GetDensity() << " v:" << heme->inlet.GetVelocity() << " O: p: "
              //    << heme->outlet.GetDensity() << " v:" << heme->outlet.GetVelocity() << std::endl;
              //std::cout << "Zerod: I: p:" << zerod->inlet.GetDensity() << " v:" << zerod->inlet.GetVelocity()
              //    << " O: p: " << zerod->outlet.GetDensity() << " v:" << zerod->outlet.GetVelocity() << std::endl;
            }

            CPPUNIT_ASSERT_DOUBLES_EQUAL(mockheme->current_time, zerod->current_time,1e-6);
          }
          void TestCHemeRun()
          {
            int argc;
            const char* argv[9];
            argc = 9;
            argv[0] = "hemelb";
            argv[2] = "four_cube_multiscale.xml";
            argv[1] = "-in";
            argv[3] = "-i";
            argv[4] = "1";
            argv[5] = "-s";
            argv[6] = "1";
            argv[7] = "-ss";
            argv[8] = "1111";
            FolderTestFixture::setUp();
            CopyResourceToTempdir("four_cube_multiscale.xml");
            CopyResourceToTempdir("four_cube.gmy");
            hemelb::configuration::CommandLine options(argc, argv);
            MockIntercommunicator intercomms(*pbuffer);
            heme = new MultiscaleSimulationMaster<MockIntercommunicator>(options,intercomms);
            // Mock out the behaviour of the simulation master iteration, but with the other model linked in.
            while (heme->GetState()->GetTimeStepsPassed() <= heme->GetState()->GetTotalTimeSteps())
            {
              zerod->Simulate();
              heme->DoTimeStep();
              //std::cout << "Heme: I: p:" << heme->inlet.GetDensity() << " v:" << heme->inlet.GetVelocity() << " O: p: "
              //    << heme->outlet.GetDensity() << " v:" << heme->outlet.GetVelocity() << std::endl;
              std::cout << "Zerod: @" << zerod->current_time << " I: p:"  << zerod->inlet.GetDensity() << " v:" << zerod->inlet.GetVelocity()
                  << " O: p: " << zerod->outlet.GetDensity() << " v:" << zerod->outlet.GetVelocity() << std::endl;
            }
            heme->Finalise();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(heme->GetState()->GetTime(), zerod->current_time,1e-6);
            delete heme;
          }
          MockHemeLB<MockIntercommunicator> *mockheme;
          MultiscaleSimulationMaster<MockIntercommunicator> *heme;
          Mock0DModel<MockIntercommunicator> *zerod;
          std::map<std::string, double> *pbuffer;
      };
      //class
      CPPUNIT_TEST_SUITE_REGISTRATION(MockIntercommunicatorTests);
    } //multiscale
  } //unittests
} //hemelb
#endif  HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
