#ifndef HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H
#define HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H

#include <cppunit/TestFixture.h>
#include "reporting/Reporter.h"
#include "reporting/Reporter.hpp"
#include "reporting/Timers.h"
#include <functional>
#include <iomanip>
namespace hemelb
{
  namespace unittests
  {
    namespace reporting
    {
      using namespace hemelb::reporting;

      typedef TimersBase<ClockMock, MPICommsMock> TimersMock;
      typedef ReporterBase<TimersMock, WriterMock, MPICommsMock> ReporterMock;

      class ReporterTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(ReporterTests);
          CPPUNIT_TEST(TestInit);
          CPPUNIT_TEST(TestCycle);
          CPPUNIT_TEST(TestImage);
          CPPUNIT_TEST(TestSnapshot);
          CPPUNIT_TEST(TestMainReport);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            timers = new TimersMock();
            reporter = new ReporterMock("mock_path", "exampleinputfile", 1234, *timers);
          }

          void tearDown()
          {
            delete reporter;
            delete timers;
          }

          void TestInit()
          {
            CPPUNIT_ASSERT_EQUAL((size_t) 8,
                                 reporter->Results().back().find("config file:\n exampleinputfile\n"));
          }

          void TestCycle()
          {
            can_count(std::mem_fun(&ReporterMock::Cycle), std::string("cycle id"));
          }
          void TestImage()
          {
            can_count(std::mem_fun(&ReporterMock::Image), std::string("Image written"));
          }
          void TestSnapshot()
          {
            can_count(std::mem_fun(&ReporterMock::Snapshot), std::string("Snapshot written"));
          }

          void can_count(std::mem_fun_t<void, ReporterMock> method,const std::string &message)
          {
            for (int times = 0; times < 5; times++)
            {
              method(reporter);
              std::stringstream expectation;
              expectation << message << ": " << times + 1 << std::endl;
              back_has_substring(expectation.str());
            }
          }

          void back_has_substring(const std::string &expectation)
          {
            std::string &actual = reporter->Results().back();
            CPPUNIT_ASSERT_MESSAGE(expectation + " : " + actual,
                                   std::string::npos != actual.find(expectation));
            reporter->Results().pop_back();
          }
          void TestMainReport()
          {
            // Mock up some timings
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++)
            {
              for (unsigned int j = 0; j < i; j++)
              {
                (*timers)[i].Start();
                (*timers)[i].Stop();
              }
            }
            timers->Reduce(); // invoke the Timers MPI mock
            reporter->Cycle();
            reporter->Cycle();
            reporter->Snapshot();
            reporter->Snapshot();
            reporter->Snapshot();
            reporter->Snapshot();
            reporter->Image();
            reporter->Image();
            reporter->Image();
            for (unsigned int step = 0; step < 1000; step++)
            {
              reporter->TimeStep();
            }
            // Ok, we have our fixture -- now execute it.
            reporter->Write();
            // now we validate the lines from the report
            for (int j = Timers::numberOfTimers - 1; j >= 0; j--)
            {
              table_expectation(j);
            }
            back_has_substring("Local \tMin \tMean \tMax");
            back_has_substring("Per-proc timing data");
            for (int j = 4; j >= 0; j--)
            {
              std::stringstream expectation;
              expectation << "rank: " << j << ", fluid sites: " << j * 1000 << std::flush;
              back_has_substring(expectation.str());
            }
            back_has_substring("Sub-domains info:");
            back_has_substring("total time (s):                            0");
            reporter->Results().pop_back();
            reporter->Results().pop_back();
            reporter->Results().pop_back();
            back_has_substring("time steps per cycle: 500");
            back_has_substring("time steps per second: 11.111");
            back_has_substring("cycles and total time steps: 2, 1000");
            back_has_substring("fluid sites: 1234");
            back_has_substring("topology depths checked: 3");
            back_has_substring("threads: 5, machines checked: 4");
          }

          void table_expectation(unsigned int row)
          {
            std::stringstream expectation;
            expectation << std::setprecision(3);
            double normn = 1.0;
            if (row == Timers::snapshot)
            {
              normn = 4.0;
            }
            if (row == Timers::visualisation)
            {
              normn = 3.0;
            }
            if (row == Timers::lb || row == Timers::mpiSend || row == Timers::mpiWait)
            {
              normn = 2.0;
            }
            expectation << timerNames[row] << "\t\t" << row * 10.0 << "\t" << row * 15.0 / normn
                << "\t" << row * 10.0 / normn << "\t" << row * 5.0 / normn << std::flush;
            back_has_substring(expectation.str());
          }

        private:
          ReporterMock *reporter;
          TimersMock *timers;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(ReporterTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H_ */
