#ifndef HEMELB_UNITTESTS_CONFIGURATION_SIMCONFIGTESTS_H
#define HEMELB_UNITTESTS_CONFIGURATION_SIMCONFIGTESTS_H
#include "configuration/SimConfig.h"
#include "unittests/resources/Resource.h"
namespace hemelb{
  namespace unittests{
    namespace configuration{
      using namespace hemelb::configuration;
      using namespace resources;
      class SimConfigTests : public CppUnit::TestFixture
      {
        CPPUNIT_TEST_SUITE( SimConfigTests );
        CPPUNIT_TEST( TestConstruct );
        CPPUNIT_TEST_SUITE_END();
        public:
        void setUp(){
          exemplar=Resource("config.xml").Path();
          config = SimConfig::Load( exemplar.c_str());
        }
        void tearDown(){
          delete config;
        }

        void TestConstruct(){
          // smoke test the configuration as having loaded OK
          CPPUNIT_ASSERT_EQUAL(3lu,config->NumCycles);
        }
        private:
        std::string exemplar;
        SimConfig * config;
      };
      CPPUNIT_TEST_SUITE_REGISTRATION( SimConfigTests );
    }
  }
}
#endif // ONCE
