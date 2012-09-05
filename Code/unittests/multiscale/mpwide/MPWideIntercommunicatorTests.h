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

        class MPWideIntercommunicatorTests : public FolderTestFixture
        {
            CPPUNIT_TEST_SUITE (MPWideIntercommunicatorTests);
            CPPUNIT_TEST (testMPWideApplication);CPPUNIT_TEST_SUITE_END();

          public:
            void setUp()
            {
              pbuffer = new std::map<std::string, double>();
              std::map<std::string, double> &buffer = *pbuffer;

              LBorchestration = new std::map<std::string, bool>();
              std::map<std::string, bool> &LBorchestrationInstance = *LBorchestration;
              LBorchestrationInstance["boundary1_pressure"] = false;
              LBorchestrationInstance["boundary2_pressure"] = false;
              LBorchestrationInstance["boundary1_velocity"] = true;
              LBorchestrationInstance["boundary2_velocity"] = true;

              mockheme = new InterCommunicatingHemeLB<MPWideIntercommunicator>(25.0, 0.2, buffer, LBorchestrationInstance);
            }
            void tearDown()
            {
              delete mockheme;
              delete pbuffer;
              delete LBorchestration;
            }
          private:
            InterCommunicatingHemeLB<MPWideIntercommunicator> *mockheme;
            std::map<std::string, double> *pbuffer;
            std::map<std::string, bool> *LBorchestration;

            void testMPWidePresent()
            {
              std::cout << "IP address for localhost is: " << MPW_DNSResolve("localhost") << std::endl;
              std::cout << "MPWide is present." << std::endl;
            }
            void testMPWideInit()
            {

            }
            void testMPWideApplication()
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
              MPWideIntercommunicator intercomms(*pbuffer, *LBorchestration);
              //std::cout << "Spawning HemeLB..." << std::endl;

              MultiscaleSimulationMaster<MPWideIntercommunicator> *heme;
              heme = new MultiscaleSimulationMaster<MPWideIntercommunicator>(options, intercomms);
              // Mock out the behaviour of the simulation master iteration, but with the other model linked in.
              //std::cout << "HemeLB about to be run..." << std::endl;
              while (heme->GetState()->GetTime() < 20.0)
              {
                heme->DoTimeStep();
                //std::cout << "Step taken, going to incrementSharedTime." << std::endl;
                intercomms.UnitTestIncrementSharedTime(); //simple hack func that mocks a 1.0 increase in the 'other' simulation.
              }
              heme->Finalise();
              CPPUNIT_ASSERT_DOUBLES_EQUAL(heme->GetState()->GetTime(), 20.0, 1e-6);
              //CPPUNIT_ASSERT_DOUBLES_EQUAL(zerod->currentTime, 20.5, 1e-6); // does one more step, where it sets the shared time.
              delete heme;
            }

        };
        //class
        CPPUNIT_TEST_SUITE_REGISTRATION (MPWideIntercommunicatorTests);

      } //MPWide
    } //multiscale
  } //unittests
} //hemelb
#endif  //HEMELB_UNITTEST_MULTISCALE_MPWIDEINTERCOMMUNICATORTESTS_H
