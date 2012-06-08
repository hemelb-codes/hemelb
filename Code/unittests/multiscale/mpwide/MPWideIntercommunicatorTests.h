#ifndef HEMELB_UNITTESTS_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATORTESTS_H
#define HEMELB_UNITTESTS_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATORTESTS_H
#include "unittests/multiscale/MockMPWide.h"
#include "unittests/multiscale/MockIntercommunicand.h"
#include <resources/Resource.h>
#include "multiscale/mpwide/MPWideIntercommunicator.h"
#include "multiscale/MultiscaleSimulationMaster.h"


namespace hemelb
{
  namespace unittests
  {
    namespace multiscale
    {
      namespace mpwide
      {
        /***
        * Mock intercommunicating entity which looks a bit like a HemeLB conceptually
        * It has an input, and an output, and the flow rate depends on the difference in pressures.
        */
        template<class IntercommuniatorImplementation> class IIHemeLB
        {
          public:
            IIHemeLB(double spaceResolution, double timeResolution, std::map<std::string, double> & buffer,std::map<std::string,bool> &orchestration) :
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
       
      class MPWideIntercommunicatorTests : public FolderTestFixture
      {
        CPPUNIT_TEST_SUITE(MPWideIntercommunicatorTests);
        //CPPUNIT_TEST(TestMPWidePresent);
        //CPPUNIT_TEST(TestMPWideInit);
        CPPUNIT_TEST(TestMPWideApplication);
        CPPUNIT_TEST_SUITE_END();     
        
        public:
          void setUp()
          {
            pbuffer = new std::map<std::string, double>();
            std::map<std::string, double> &buffer = *pbuffer;
            
            orchestrationLB=new std::map<std::string,bool>();
            std::map<std::string,bool> &rorchestrationLB=*orchestrationLB;
            rorchestrationLB["boundary1_pressure"] = false;
            rorchestrationLB["boundary2_pressure"] = false;
            rorchestrationLB["boundary1_velocity"] = true;
            rorchestrationLB["boundary2_velocity"] = true;
            
            mockheme = new IIHemeLB<MPWideIntercommunicator>(25.0, 0.2, buffer,rorchestrationLB);
          }
          void tearDown()
          {
            delete mockheme;
            delete pbuffer;
            //delete orchestrationLB;
          }
        private:
          IIHemeLB<MPWideIntercommunicator> *mockheme;
          MultiscaleSimulationMaster<MPWideIntercommunicator> *heme;
          std::map<std::string,double> *pbuffer;
          std::map<std::string,bool> *orchestrationLB;
          
          void TestMPWidePresent()
          {
             std::cout << "IP address for localhost is: " << MPW_DNSResolve("localhost") << std::endl;
             std::cout << "MPWide is present." << std::endl;
          }
          void TestMPWideInit()
          {
            
          }
          void TestMPWideApplication()
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
            MPWideIntercommunicator intercomms(*pbuffer,*orchestrationLB);
            std::cout << "Spawning HemeLB..." << std::endl;
            heme = new MultiscaleSimulationMaster<MPWideIntercommunicator>(options, intercomms);
            // Mock out the behaviour of the simulation master iteration, but with the other model linked in.
            std::cout << "HemeLB about to be run..." << std::endl;
            while (heme->GetState()->GetTime() < 20.0)
            {
              heme->DoTimeStep();
              std::cout << "Step taken, going to incrementSharedTime." << std::endl;
              intercomms.UnitTestIncrementSharedTime(); //simple hack func that mocks a 1.0 increase in the 'other' simulation.
            }
            heme->Finalise();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(heme->GetState()->GetTime(), 20.0, 1e-6);
            //CPPUNIT_ASSERT_DOUBLES_EQUAL(zerod->currentTime, 20.5, 1e-6); // does one more step, where it sets the shared time.
            delete heme;
          }       
         
      };
      //class
      CPPUNIT_TEST_SUITE_REGISTRATION(MPWideIntercommunicatorTests);
      
      } //MPWide
    } //multiscale
  } //unittests
} //hemelb
#endif  //HEMELB_UNITTEST_MULTISCALE_MPWIDEINTERCOMMUNICATORTESTS_H
