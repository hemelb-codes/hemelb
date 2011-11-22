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
          CPPUNIT_TEST(TestInit);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            timers = new TimersMock();
            reporter = new ReporterBase<TimersMock,WriterMock>("mock_path","exampleinputfile",1000,*timers);
          }

          void tearDown()
          {
            delete reporter;
            delete timers;
          }

          void TestInit()
          {
            CPPUNIT_ASSERT_EQUAL((size_t) 8,reporter->Results().back().find("config file:\n exampleinputfile\n"));
          }

        private:
          ReporterBase<TimersMock,WriterMock> *reporter;
          TimersMock *timers;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(ReporterTests);


    }
  }
}

#endif /* HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H_ */
