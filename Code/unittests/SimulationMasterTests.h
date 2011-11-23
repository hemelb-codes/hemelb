#ifndef HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H
#define HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H

#include <cppunit/TestFixture.h>



namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      /**
       * Class to test the simulation master.
       */
      class SimulationMasterTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE( SimulationMasterTests );
          CPPUNIT_TEST( TestReport );
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {

          }

          void tearDown()
          {

          }


          void TestReport()
          {

          }

        private:

      };
      CPPUNIT_TEST_SUITE_REGISTRATION( SimulationMasterTests );
    }
  }
}

#endif /* HEMELB_UNITTESTS_SIMULATIONMASTERTESTS_H_ */
