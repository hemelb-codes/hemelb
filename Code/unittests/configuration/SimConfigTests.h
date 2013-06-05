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
          CPPUNIT_TEST (Test_0_2_0_Write);
          CPPUNIT_TEST (Test_0_2_1_Write);CPPUNIT_TEST_SUITE_END();
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
            SimConfig *config = SimConfig::Load(Resource("config0_2_0.xml").Path().c_str());
            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
            CPPUNIT_ASSERT_EQUAL(60.0 / (70.0 * 1000), config->GetTimeStepLength());

            CPPUNIT_ASSERT_EQUAL(60.0 / 70.0,
                                 static_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0])->GetPeriod());
            delete config;
          }
          void Test_0_2_1_Read()
          {
            // smoke test the configuration as having loaded OK
            SimConfig *config = SimConfig::Load(Resource("config.xml").Path().c_str());
            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
            CPPUNIT_ASSERT_EQUAL(0.0001, config->GetTimeStepLength());

            CPPUNIT_ASSERT_EQUAL(0.6,
                                 static_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0])->GetPeriod());
            delete config;
          }
          void Test_0_2_0_Write()
          {
            FolderTestFixture::setUp();
            //Round trip the config twice.
            CopyResourceToTempdir("config0_2_0.xml");
            SimConfig *config = SimConfig::Load("config0_2_0.xml");
            config->Save("config0_2_0b.xml");
            delete config;
            config = SimConfig::Load("config0_2_0b.xml");
            config->Save("config0_2_0c.xml");
            delete config;
            config = SimConfig::Load("config0_2_0c.xml");

            // Assert the values are correct.
            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0 / (70.0 * 1000), config->GetTimeStepLength(),1e-6);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(60.0 / 70.0,
                                 static_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0])->GetPeriod(),1e-6);
            FolderTestFixture::tearDown();
            delete config;
          }
          void Test_0_2_1_Write()
          {
            FolderTestFixture::setUp();
            //Round trip the config twice.
            CopyResourceToTempdir("config.xml");
            SimConfig *config = SimConfig::Load("config.xml");
            config->Save("config0_2_1b.xml");
            delete config;
            config = SimConfig::Load("config0_2_1b.xml");
            config->Save("config0_2_1c.xml");
            delete config;
            config = SimConfig::Load("config0_2_1c.xml");

            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
            CPPUNIT_ASSERT_EQUAL(0.0001, config->GetTimeStepLength());

            CPPUNIT_ASSERT_EQUAL(0.6,
                                 static_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0])->GetPeriod());
            FolderTestFixture::tearDown();
            delete config;
          }
        private:
          std::string exemplar;
      };
      CPPUNIT_TEST_SUITE_REGISTRATION (SimConfigTests);
    }
  }
}
#endif // ONCE
