#ifndef HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H
#define HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H

#include <cppunit/TestFixture.h>
#include "reporting/Reporter.h"
#include "reporting/Reporter.hpp"
#include "reporting/Timers.h"

namespace hemelb
{
  namespace unittests
  {
    namespace reporting
    {
      using namespace hemelb::reporting;

      typedef TimersBase<ClockMock, MPICommsMock> TimersMock;

      class ReporterTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(ReporterTests);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            timers = new Timers();
            reporter = new ReporterBase<TimersMock>("examplepath","exampleinputfile",1000,*timers);
          }

          void tearDown()
          {
            delete reporter;
            delete timers;
          }



        private:
          ReporterBase<TimersMock> *reporter;
          TimersMock *timers;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(ReporterTests);


    }
  }
}

#endif /* HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H_ */
