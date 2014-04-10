// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATORTESTS_H
#define HEMELB_UNITTESTS_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATORTESTS_H
#include "unittests/multiscale/MockMPWide.h"
#include "unittests/multiscale/MockIntercommunicand.h"
#include "unittests/multiscale/mpwide/IntercommunicatingHemeLB.h"
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

        class MPWideIntercommunicatorTests : public helpers::FolderTestFixture
        {
            CPPUNIT_TEST_SUITE (MPWideIntercommunicatorTests);
            // TODO The below test is v important and is not currently being run.
            // CPPUNIT_TEST (testMPWideApplication);
            // CPPUNIT_TEST(testMPWidePresent);
            CPPUNIT_TEST_SUITE_END();

          public:
            void setUp()
            {
              helpers::FolderTestFixture::setUp();

              pbuffer = new std::map<std::string, double>();

              LBorchestration = new std::map<std::string, bool>();
              (*LBorchestration)["boundary1_pressure"] = false;
              (*LBorchestration)["boundary2_pressure"] = false;
              (*LBorchestration)["boundary1_velocity"] = true;
              (*LBorchestration)["boundary2_velocity"] = true;

              std::string configPath = "../../../config_files/MPWSettings.cfg";
              mockheme = new InterCommunicatingHemeLB<MPWideIntercommunicator>(25.0,
                                                                               0.2,
                                                                               *pbuffer,
                                                                               *LBorchestration,
                                                                               configPath);
            }

            void tearDown()
            {
              delete mockheme;
              delete pbuffer;
              delete LBorchestration;
              helpers::FolderTestFixture::tearDown();
            }

          private:
            InterCommunicatingHemeLB<MPWideIntercommunicator> *mockheme;
            std::map<std::string, double> *pbuffer;
            std::map<std::string, bool> *LBorchestration;

            void testMPWidePresent()
            {
              std::string host = "localhost";
              std::cout << "IP address for localhost is: " << MPW_DNSResolve(const_cast<char*>(host.c_str()))
                  << std::endl;
              std::cout << "MPWide is present." << std::endl;
            }
            void testMPWideInit()
            {
              // TODO This test needs writing.
            }
            void testMPWideApplication()
            {
              int argc;
              const char* argv[7];
              argc = 7;
              argv[0] = "hemelb";
              argv[2] = "four_cube_multiscale.xml";
              argv[1] = "-in";
              argv[3] = "-i";
              argv[4] = "1";
              argv[5] = "-ss";
              argv[6] = "1111";

              CopyResourceToTempdir("four_cube_multiscale.xml");
              CopyResourceToTempdir("four_cube.gmy");
              hemelb::configuration::CommandLine options(argc, argv);
              std::string configPath = "../../../config_files/MPWSettings.cfg";
              MPWideIntercommunicator intercomms(Comms().OnIORank(), *pbuffer, *LBorchestration, configPath);

              MultiscaleSimulationMaster<MPWideIntercommunicator> heme(options, Comms(), intercomms);
              // Mock out the behaviour of the simulation master iteration, but with the other model linked in.
              //std::cout << "HemeLB about to be run..." << std::endl;
              while (heme.GetState()->GetTime() < 20.0)
              {
                heme.DoTimeStep();
                //std::cout << "Step taken, going to incrementSharedTime." << std::endl;
                intercomms.UnitTestIncrementSharedTime(); //simple hack func that mocks a 1.0 increase in the 'other' simulation.
              }
              heme.Finalise();
              CPPUNIT_ASSERT_DOUBLES_EQUAL(heme.GetState()->GetTime(), 20.0, 1e-6);
              //CPPUNIT_ASSERT_DOUBLES_EQUAL(zerod->currentTime, 20.5, 1e-6); // does one more step, where it sets the shared time.
              FolderTestFixture::tearDown();
            }

        };
        //class
        CPPUNIT_TEST_SUITE_REGISTRATION (MPWideIntercommunicatorTests);

      } //MPWide
    } //multiscale
  } //unittests
} //hemelb
#endif  //HEMELB_UNITTEST_MULTISCALE_MPWIDEINTERCOMMUNICATORTESTS_H
