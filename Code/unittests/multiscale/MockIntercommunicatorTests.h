#ifndef HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
#define HEMELB_UNITTESTS_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
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
      /***
       * Mock intercommunicating entity which looks a bit like a HemeLB conceptually
       * It has an input, and an output, and the flow rate depends on the difference in pressures.
       */
      template<class IntercommuniatorImplementation> class MockHemeLB
      {
        public:
          MockHemeLB(double spaceResolution, double timeResolution, std::map<std::string, double> & buffer) :
              inlet(), outlet(), inOutLetType("inoutlet"), intercomms(buffer), timeResolution(timeResolution), spaceResolution(spaceResolution), currentTime(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inOutLetType.template RegisterSharedValue<double>("pressure");
            inOutLetType.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inOutLetType, inlet, "boundary1");
            intercomms.RegisterIntercommunicand(inOutLetType, outlet, "boundary2");
          }
          MockIntercommunicand inlet;
          MockIntercommunicand outlet;
          typename IntercommuniatorImplementation::IntercommunicandTypeT inOutLetType;
          IntercommuniatorImplementation intercomms;
          double timeResolution;
          double spaceResolution;
          double currentTime;
          void DoLB()
          {
            double resistance = 10.0;
            double velocity = (inlet.GetPressure() - outlet.GetPressure()) / resistance;
            inlet.SetVelocity(velocity);
            outlet.SetVelocity(velocity);
          }

          void Simulate()
          {
            if (!intercomms.ShouldAdvance())
              return;
            intercomms.GetFromMultiscale();
            DoLB();
            currentTime += timeResolution;
            intercomms.AdvanceTime(currentTime);
            intercomms.SendToMultiscale();
          }

      };
      /***
       * Mock intercommunicating entity which is a tank of water.
       * It has a pressure difference which drops depending on the flow rate.
       */
      template<class IntercommuniatorImplementation> class Mock0DModel
      {
        public:
          Mock0DModel(double spaceResolution, double timeResolution, std::map<std::string, double> & buffer) :
              inlet(), outlet(), inOutLetType("inoutlet"), intercomms(buffer), timeResolution(timeResolution), spaceResolution(spaceResolution), currentTime(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inOutLetType.template RegisterSharedValue<double>("pressure");
            inOutLetType.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inOutLetType, inlet, "boundary2");
            intercomms.RegisterIntercommunicand(inOutLetType, outlet, "boundary1");
          }
          MockIntercommunicand inlet;
          MockIntercommunicand outlet;
          typename IntercommuniatorImplementation::IntercommunicandTypeT inOutLetType;
          IntercommuniatorImplementation intercomms;
          double timeResolution;
          double spaceResolution;
          double currentTime;
          double GetOutletPressure()
          {
            return outlet.GetPressure();
          }
          void Do1D()
          {
            double capacitance = 10.0;
            double deltap = -1.0 * inlet.GetVelocity() / capacitance;

            outlet.SetPressure(outlet.GetPressure() + deltap);
          }

          void Simulate()
          {
            if (!intercomms.ShouldAdvance())
              return;
            intercomms.GetFromMultiscale();
            Do1D();
            currentTime += timeResolution;
            intercomms.AdvanceTime(currentTime);
            intercomms.SendToMultiscale();
          }

      };

      // Useful for debugging to have this.
      std::ostream & operator <<(std::ostream & stream, std::map<std::string, double> buffer)
      {
        stream << " { ";
        for (std::map<std::string, double>::iterator entry = buffer.begin(); entry != buffer.end(); entry++)
        {
          stream << entry->first << " : " << entry->second << " , ";
        }
        stream << " } ";
        stream << std::flush;
        return stream;
      }

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
            buffer["boundary1_velocity"] = 0.1;
            buffer["boundary2_velocity"] = 0.1;
            buffer["shared_time"] = 0.0;
            mockheme = new MockHemeLB<MockIntercommunicator>(25.0, 0.2, buffer);
            zerod = new Mock0DModel<MockIntercommunicator>(10.0, 0.5, buffer);

          }
          void tearDown()
          {
            delete zerod;
            delete mockheme;
            delete pbuffer;
          }
        private:
          void TestCRRun()
          {
            // In this example, we imagine we are coupling a mock hemelb with a zero-d capacitative model.
            // Our example intercommunicator has been hard coded with how to share the values.
            // In principle, though, this would be defined by whatever configuration file format the implementing intercommunicator wants to use.
            for (int steps = 0; steps < 100; steps++)
            {
              mockheme->Simulate();
              zerod->Simulate();
            }
            // The asserted value here is a result of discretised exponential decay.
            // I have not thought it worth determining this by a formula.
            CPPUNIT_ASSERT_DOUBLES_EQUAL(80.3379435171394, zerod->GetOutletPressure(), 1e-6);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(mockheme->currentTime, zerod->currentTime, 1e-6);
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
            heme = new MultiscaleSimulationMaster<MockIntercommunicator>(options, intercomms);
            // Mock out the behaviour of the simulation master iteration, but with the other model linked in.
            while (heme->GetState()->GetTimeStep() <= heme->GetState()->GetTotalTimeSteps())
            {
              zerod->Simulate();
              heme->DoTimeStep();
            }
            // In advancing 100 time steps, at 0.2 s per time step, with a 0d model at 0.5s per time step
            // the 0d model with execute 100*2/5=40 times.
            // Each time, the pressure difference will drop by 0.1*0.1=0.01 mmHg.
            // So the final pressure will be 81.0-40*0.01=80.6 mmHg
            CPPUNIT_ASSERT_DOUBLES_EQUAL(80.6, zerod->GetOutletPressure(), 1e-6);
            heme->Finalise();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(heme->GetState()->GetTime(), zerod->currentTime, 1e-6);
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
#endif  //HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
