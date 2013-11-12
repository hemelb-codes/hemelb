// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_CONFIGURATION_SIMCONFIGTESTS_H
#define HEMELB_UNITTESTS_CONFIGURATION_SIMCONFIGTESTS_H
#include "configuration/SimConfig.h"
#include "resources/Resource.h"
#include "unittests/helpers/FolderTestFixture.h"
namespace hemelb
{
  namespace unittests
  {
    namespace configuration
    {
      using namespace hemelb::configuration;
      using namespace resources;
      using namespace helpers;
      class SimConfigTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (SimConfigTests);
          CPPUNIT_TEST (Test_0_2_0_Read);
          CPPUNIT_TEST (Test_0_2_1_Read);
//          CPPUNIT_TEST (Test_0_2_0_Write);
//          CPPUNIT_TEST (Test_0_2_1_Write);
//          CPPUNIT_TEST (TestVelocityInletsWrite);
          CPPUNIT_TEST (TestXMLFileContent);CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {

          }
          void tearDown()
          {

          }

          void Test_0_2_0_Read()
          {
            // smoke test the configuration as having loaded OK
            SimConfig* config = SimConfig::New(Resource("config0_2_0.xml").Path());
            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
            CPPUNIT_ASSERT_EQUAL(0.0001, config->GetTimeStepLength());
            lb::iolets::InOutLetCosine* inlet = dynamic_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0]);
            CPPUNIT_ASSERT(inlet != NULL);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(6000.0, inlet->GetPeriod(), 1e-6);
          }

          void Test_0_2_1_Read()
          {
            // smoke test the configuration as having loaded OK
            SimConfig* config = SimConfig::New(Resource("config.xml").Path());
            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
            CPPUNIT_ASSERT_EQUAL(0.0001, config->GetTimeStepLength());
            lb::iolets::InOutLetCosine* inlet = dynamic_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0]);
            CPPUNIT_ASSERT(inlet != NULL);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(6000.0, inlet->GetPeriod(), 1e-6);
          }
//          void Test_0_2_0_Write()
//          {
//            FolderTestFixture::setUp();
//            //Round trip the config twice.
//            CopyResourceToTempdir("config0_2_0.xml");
//            SimConfig *config = SimConfig::Load("config0_2_0.xml");
//            config->Save("config0_2_0b.xml");
//            delete config;
//            config = SimConfig::Load("config0_2_0b.xml");
//            config->Save("config0_2_0c.xml");
//            delete config;
//            config = SimConfig::Load("config0_2_0c.xml");
//
//            // Assert the values are correct.
//            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
//            CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0 / (70.0 * 1000), config->GetTimeStepLength(), 1e-6);
//
//            CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0 / 70.0,
//                                         static_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0])->GetPeriod(),
//                                         1e-6);
//            FolderTestFixture::tearDown();
//            delete config;
//          }
//          void Test_0_2_1_Write()
//          {
//            FolderTestFixture::setUp();
//            //Round trip the config twice.
//            CopyResourceToTempdir("config.xml");
//            SimConfig *config = SimConfig::Load("config.xml");
//            config->Save("config0_2_1b.xml");
//            delete config;
//            config = SimConfig::Load("config0_2_1b.xml");
//            config->Save("config0_2_1c.xml");
//            delete config;
//            config = SimConfig::Load("config0_2_1c.xml");
//
//            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
//            CPPUNIT_ASSERT_EQUAL(0.0001, config->GetTimeStepLength());
//
//            CPPUNIT_ASSERT_EQUAL(0.6,
//                                 static_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0])->GetPeriod());
//            FolderTestFixture::tearDown();
//            delete config;
//          }

//          void TestVelocityInletsWrite()
//          {
//            FolderTestFixture::setUp();
//            //Round trip the config twice.
//            CopyResourceToTempdir("config_new_velocity_inlets.xml");
//            SimConfig *config = SimConfig::Load("config_new_velocity_inlets.xml");
//            config->Save("config_new_velocity_inlets_b.xml");
//            delete config;
//            config = SimConfig::Load("config_new_velocity_inlets_b.xml");
//            config->Save("config_new_velocity_inlets_c.xml");
//            delete config;
//            config = SimConfig::Load("config_new_velocity_inlets_c.xml");
//
//            lb::iolets::InOutLetWomersleyVelocity* inlet =
//                dynamic_cast<lb::iolets::InOutLetWomersleyVelocity*>(config->GetInlets()[0]);
//            assert(inlet);
//
//            CPPUNIT_ASSERT_EQUAL(10.0, inlet->GetRadius());
//            CPPUNIT_ASSERT_EQUAL(2.5, inlet->GetPressureGradientAmplitude());
//            CPPUNIT_ASSERT_EQUAL(LatticeTime(5), inlet->GetPeriod());
//            CPPUNIT_ASSERT_EQUAL(2.0, inlet->GetWomersleyNumber());
//            FolderTestFixture::tearDown();
//            delete config;
//          }

          void TestXMLFileContent()
          {
            FolderTestFixture::setUp();
            //Round trip the config twice.
            CopyResourceToTempdir("config.xml");
            SimConfig* config = SimConfig::New("config.xml");

            CPPUNIT_ASSERT_EQUAL(80.0, config->GetInitialPressure());
          }

        private:
          std::string exemplar;
      };
      CPPUNIT_TEST_SUITE_REGISTRATION (SimConfigTests);
    }
  }
}
#endif // ONCE
