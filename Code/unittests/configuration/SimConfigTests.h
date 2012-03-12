#ifndef HEMELB_UNITTESTS_CONFIGURATION_SIMCONFIGTESTS_H
#define HEMELB_UNITTESTS_CONFIGURATION_SIMCONFIGTESTS_H
#include "configuration/SimConfig.h"
#include "resources/Resource.h"
namespace hemelb
{
  namespace unittests
  {
    namespace configuration
    {
      using namespace hemelb::configuration;
      using namespace resources;
      class SimConfigTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(SimConfigTests);
          CPPUNIT_TEST(Test_0_2_0);
          CPPUNIT_TEST(Test_0_2_1);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {

          }
          void tearDown()
          {

          }

          void Test_0_2_0()
          {
            // smoke test the configuration as having loaded OK
            SimConfig *config=SimConfig::Load(Resource("config0_2_0.xml").Path().c_str());
            CPPUNIT_ASSERT_EQUAL(3000lu, config->TotalTimeSteps);
            CPPUNIT_ASSERT_EQUAL(60.0/(70.0*1000), config->TimeStepLength);
            delete config;
          }
          void Test_0_2_1()
          {
            // smoke test the configuration as having loaded OK
            SimConfig *config=SimConfig::Load(Resource("config.xml").Path().c_str());
            CPPUNIT_ASSERT_EQUAL(3000lu, config->TotalTimeSteps);
            CPPUNIT_ASSERT_EQUAL(0.0001, config->TimeStepLength);
            delete config;
          }
        private:
          std::string exemplar;
      };
      CPPUNIT_TEST_SUITE_REGISTRATION(SimConfigTests);
    }
  }
}
#endif // ONCE
