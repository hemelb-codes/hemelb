#ifndef HEMELB_UNITTESTS_REPORTING_MOCKS_H
#define HEMELB_UNITTESTS_REPORTING_MOCKS_H
namespace hemelb
{
  namespace unittests
  {
    namespace reporting
    {
      class ClockMock
      {
        public:
          ClockMock() :
              fakeTime(0)
          {
          }
          ;
        protected:
          double CurrentTime()
          {
            fakeTime += 10.0;
            return fakeTime;
          }
        private:
          double fakeTime;
      };

      class MPICommsMock
      {
        public:
          MPICommsMock() :
              calls(1)
          {
          }
          ;
        protected:
          int Reduce(double *sendbuf,
                     double *recvbuf,
                     int count,
                     MPI_Datatype datatype,
                     MPI_Op op,
                     int root,
                     MPI_Comm comm)
          {
            CPPUNIT_ASSERT_EQUAL(33, count);
            for (int i = 0; i < count; i++)
            {
              CPPUNIT_ASSERT_EQUAL(10.0 * i, sendbuf[i]);
              recvbuf[i] = 5.0 * i * calls;
            }
            calls++;
            return 0;
          }
          proc_t GetProcessorCount()
          {
            return 5;
          }
          unsigned int GetMachineCount()
          {
            return 4;
          }
          int GetDepths()
          {
            return 3;
          }
        private:
          unsigned int calls;
      };
    }
  }
}
#endif // HEMELB_UNITTESTS_REPORTING_MOCKS_H
