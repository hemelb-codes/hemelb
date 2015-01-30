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
            CPPUNIT_ASSERT(inlet != nullptr);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(6000.0, inlet->GetPeriod(), 1e-6);

            // Check that in the absence of the <monitoring> XML element things get initiliased properly
            const hemelb::configuration::SimConfig::MonitoringConfig* monConfig = config->GetMonitoringConfiguration();
            CPPUNIT_ASSERT(!monConfig->doConvergenceCheck);
            CPPUNIT_ASSERT(!monConfig->doIncompressibilityCheck);
            CPPUNIT_ASSERT(!monConfig->convergenceTerminate);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0., monConfig->convergenceRelativeTolerance, 1e-6);
          }

          void Test_0_2_1_Read()
          {
            // smoke test the configuration as having loaded OK
            SimConfig* config = SimConfig::New(Resource("config.xml").Path());
            CPPUNIT_ASSERT_EQUAL(3000lu, config->GetTotalTimeSteps());
            CPPUNIT_ASSERT_EQUAL(0.0001, config->GetTimeStepLength());
            lb::iolets::InOutLetCosine* inlet = dynamic_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0]);
            CPPUNIT_ASSERT(inlet != nullptr);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(6000.0, inlet->GetPeriod(), 1e-6);

            const hemelb::configuration::SimConfig::MonitoringConfig* monConfig = config->GetMonitoringConfiguration();
            CPPUNIT_ASSERT(monConfig->doConvergenceCheck);
            CPPUNIT_ASSERT(monConfig->doIncompressibilityCheck);
            CPPUNIT_ASSERT(monConfig->convergenceTerminate);
            CPPUNIT_ASSERT_EQUAL(1e-9, monConfig->convergenceRelativeTolerance);
            CPPUNIT_ASSERT_EQUAL(monConfig->convergenceVariable, extraction::OutputField::Velocity);
            CPPUNIT_ASSERT_EQUAL(0.01, monConfig->convergenceReferenceValue); // 1 m/s * (delta_t / delta_x) = 0.01
          }

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
