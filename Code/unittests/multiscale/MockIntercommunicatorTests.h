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
          MockHemeLB(double spaceResolution, double timeResolution, std::map<std::string, double> & buffer,std::map<std::string,bool> &orchestration) :
              inlet(81.0,0.1), outlet(79.0,0.1), inOutLetType("inoutlet"), intercomms(buffer,orchestration), timeResolution(timeResolution), spaceResolution(spaceResolution), currentTime(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inOutLetType.template RegisterSharedValue<double>("pressure");
            inOutLetType.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inOutLetType, inlet, "boundary1");
            intercomms.RegisterIntercommunicand(inOutLetType, outlet, "boundary2");
            intercomms.ShareInitialConditions();
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
            if (intercomms.DoMultiscale(currentTime))
            {
              std::cerr << "MH: " << currentTime << std::endl;
              DoLB();
              currentTime += timeResolution;
            }
          }
      };
      /***
       * Mock intercommunicating entity which is a tank of water.
       * It has a pressure difference which drops depending on the flow rate.
       */
      template<class IntercommuniatorImplementation> class Mock0DModel
      {
        public:
          Mock0DModel(double spaceResolution, double timeResolution, std::map<std::string, double> & buffer, std::map<std::string,bool> &orchestration) :
              inlet(79.0,0.1), outlet(81.0,0.1), inOutLetType("inoutlet"), intercomms(buffer,orchestration), timeResolution(timeResolution), spaceResolution(spaceResolution), currentTime(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inOutLetType.template RegisterSharedValue<double>("pressure");
            inOutLetType.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inOutLetType, inlet, "boundary2");
            intercomms.RegisterIntercommunicand(inOutLetType, outlet, "boundary1");
            intercomms.ShareInitialConditions();
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
          double GetTime()
          {
            return currentTime;
          }
          void Do1D()
          {
            double capacitance = 10.0;
            double deltap = -1.0 * inlet.GetVelocity() / capacitance;

            outlet.SetPressure(outlet.GetPressure() + deltap);
          }

          void Simulate()
          {
            if (intercomms.DoMultiscale(currentTime))
            {
              Do1D();
              currentTime += timeResolution;
              std::cerr << "0D: " << currentTime << " " << GetOutletPressure() << std::endl;
            }
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

            orchestrationOD=new std::map<std::string,bool>();
            orchestrationLB=new std::map<std::string,bool>();
            std::map<std::string,bool> &rorchestrationLB=*orchestrationLB;
            std::map<std::string,bool> &rorchestration0D=*orchestrationOD;
            rorchestrationLB["boundary1_pressure"] = false;
            rorchestrationLB["boundary2_pressure"] = false;
            rorchestrationLB["boundary1_velocity"] = true;
            rorchestrationLB["boundary2_velocity"] = true;
            rorchestration0D["boundary1_pressure"] = true;
            rorchestration0D["boundary2_pressure"] = true;
            rorchestration0D["boundary1_velocity"] = false;
            rorchestration0D["boundary2_velocity"] = false;


            mockheme = new MockHemeLB<MockIntercommunicator>(25.0, 0.2, buffer,rorchestrationLB);
            zerod = new Mock0DModel<MockIntercommunicator>(10.0, 0.5, buffer,rorchestration0D);

          }
          void tearDown()
          {
            delete zerod;
            delete mockheme;
            delete pbuffer;
            //delete orchestrationLB;
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
            CPPUNIT_ASSERT_DOUBLES_EQUAL(80.3258457641491, zerod->GetOutletPressure(), 1e-6);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(mockheme->currentTime, 20.0, 1e-6);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(zerod->currentTime, 20.5, 1e-6);
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
            MockIntercommunicator intercomms(*pbuffer,*orchestrationLB);
            heme = new MultiscaleSimulationMaster<MockIntercommunicator>(options, intercomms);
            // Mock out the behaviour of the simulation master iteration, but with the other model linked in.
            while (heme->GetState()->GetTime() < 20.0 || zerod->GetTime()< 20.0)
            {
              heme->DoTimeStep();
              zerod->Simulate();
            }
            // In advancing 100 time steps, at 0.2 s per time step, with a 0d model at 0.5s per time step
            // the 0d model with execute 100*2/5=40 times, plus one more step, the one where it communicates the previous step.
            // Each time, the pressure difference will drop by 0.1*0.1=0.01 mmHg.
            // So the final pressure will be 81.0-41*0.01=80.59 mmHg
            CPPUNIT_ASSERT_DOUBLES_EQUAL(80.59, zerod->GetOutletPressure(), 1e-6);
            heme->Finalise();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(heme->GetState()->GetTime(), 20.0, 1e-6);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(zerod->currentTime, 20.5, 1e-6); // does one more step, where it sets the shared time.
            delete heme;
          }
          MockHemeLB<MockIntercommunicator> *mockheme;
          MultiscaleSimulationMaster<MockIntercommunicator> *heme;
          Mock0DModel<MockIntercommunicator> *zerod;
          std::map<std::string, double> *pbuffer;
          std::map<std::string,bool> *orchestrationLB;
          std::map<std::string,bool> *orchestrationOD;
      };
      //class
      CPPUNIT_TEST_SUITE_REGISTRATION(MockIntercommunicatorTests);
    } //multiscale
  } //unittests
} //hemelb
#endif  //HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
