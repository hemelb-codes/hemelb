#ifndef HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H
#define HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H

#include <cppunit/TestFixture.h>
#include "reporting/Reporter.h"
#include "reporting/Reporter.hpp"
#include "reporting/Timers.h"
#include "unittests/FourCubeLatticeData.h"
#include "unittests/lbtests/LbTestsHelper.h"
#include <iomanip>
namespace hemelb
{
  namespace unittests
  {
    namespace reporting
    {
      using namespace hemelb::reporting;

      typedef TimersBase<ClockMock, MPICommsMock> TimersMock;
      typedef ReporterBase<ClockMock,  MPICommsMock, net::BroadcastMock> ReporterMock;
      typedef lb::IncompressibilityChecker<net::BroadcastMock> IncompressibilityCheckerMock;

      class ReporterTests : public CppUnit::TestFixture
      {
        CPPUNIT_TEST_SUITE( ReporterTests);
        CPPUNIT_TEST( TestInit);
        CPPUNIT_TEST( TestMainReport);CPPUNIT_TEST_SUITE_END();
        public:
        void setUp()
        {
          mockTimers = new TimersMock();
          realTimers = new reporting::Timers();
          state = new hemelb::lb::SimulationState(500, 2);
          net = new net::Net();
          latticeData = FourCubeLatticeData::Create(4, 5); // The 5 here is to match the topology size in the MPICommsMock
          lbtests::LbTestsHelper::InitialiseAnisotropicTestData(latticeData);
          latticeData->SwapOldAndNew(); //Needed since InitialiseAnisotropicTestData only initialises FOld
          incompChecker = new IncompressibilityCheckerMock(latticeData,
                                                           net,
                                                           state,
                                                           *realTimers,
                                                           10.0);
          reporter = new ReporterMock("mock_path",
                                      "exampleinputfile",
                                      latticeData->GetFluidSiteCountsOnEachProc(),
                                      1234,
                                      *mockTimers,
                                      *state,
                                      *incompChecker);
        }

        void tearDown()
        {
          delete reporter;
          delete mockTimers;
          delete realTimers;
          delete incompChecker;
          delete net;
        }

        void TestInit()
        {
          AssertValue("exampleinputfile","CONFIG");
        }

        void TestMainReport()
        {
            // Mock up some timings
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++)
            {
              for (unsigned int j = 0; j < i; j++)
              {
                (*mockTimers)[i].Start();
                (*mockTimers)[i].Stop();
              }
            }
            mockTimers->Reduce(); // invoke the Timers MPI mock
            reporter->Snapshot();
            reporter->Snapshot();
            reporter->Snapshot();
            reporter->Snapshot();
            reporter->Image();
            reporter->Image();
            reporter->Image();
            for (unsigned int step = 0; step < 1000; step++)
            {
              state->Increment();
            }
            CPPUNIT_ASSERT_EQUAL(3lu, state->GetCycleId());
            CPPUNIT_ASSERT_EQUAL(1001lu, state->GetTimeStepsPassed());
            reporter->FillDictionary();

            CheckTimingsTable();
            AssertTemplate("","{{#UNSTABLE}} unstable{{/UNSTABLE}}");
            AssertTemplate("R0S64 R1S1000 R2S2000 R3S3000 R4S4000 ","{{#PROCESSOR}}R{{RANK}}S{{SITES}} {{/PROCESSOR}}");
            AssertValue("3","IMAGES");
            AssertValue("4","SNAPSHOTS");
            AssertValue("2","CYCLES");
            AssertValue("1000","STEPS");
            AssertValue("10.000","STEPS_PER_SECOND");
            AssertValue("500","STEPS_PER_CYCLE");
            AssertValue("1234","SITES");
            AssertValue("3","DEPTHS");
            AssertValue("4","MACHINES");
        }

        private:

        void AssertValue(const std::string & expectation, const std::string &symbol){
          AssertTemplate(expectation,"{{"+symbol+"}}");
        }

        void AssertTemplate(const std::string &expectation, const std::string &ttemplate){
          ctemplate::StringToTemplateCache("TestFor"+ttemplate,ttemplate,ctemplate::DO_NOT_STRIP);
          std::string result;
          CPPUNIT_ASSERT( ctemplate::ExpandTemplate("TestFor"+ttemplate, ctemplate::DO_NOT_STRIP, &reporter->GetDictionary(), &result) );
          CPPUNIT_ASSERT_EQUAL(expectation,result);
        }


        void CheckTimingsTable()
        {

          std::stringstream expectation;
          expectation << std::setprecision(3);
          for (unsigned int row = 0; row < Timers::numberOfTimers; row++)
          {
            double normalisation = 1.0;
            if (row == Timers::snapshot)
            {
              normalisation = 4.0;
            }
            if (row == Timers::visualisation)
            {
              normalisation = 3.0;
            }
            if (row == Timers::lb || row == Timers::mpiSend || row == Timers::mpiWait || row
                == Timers::monitoring)
            {
              normalisation = 2.0;
            }
            expectation << "N" << timerNames[row] << "L" << row * 10.0 / normalisation << "MI" << row * 15.0
                / normalisation << "ME" << row * 10.0 / normalisation << "MA" << row * 5.0
                / normalisation << " " << std::flush;
          }
          AssertTemplate(expectation.str(),"{{#TIMER}}N{{NAME}}L{{LOCAL}}MI{{MIN}}ME{{MEAN}}MA{{MAX}} {{/TIMER}}");
        }

        private:
        ReporterMock *reporter;

        // We need two sets of timers because the incompressibility checker is not templated over timing policy.
        TimersMock *mockTimers;
        reporting::Timers* realTimers;

        lb::SimulationState *state;
        IncompressibilityCheckerMock *incompChecker;
        net::Net *net;
        hemelb::geometry::LatticeData *latticeData;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION( ReporterTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_REPORTING_REPORTERTESTS_H_ */
