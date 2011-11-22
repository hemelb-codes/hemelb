#ifndef HEMELB_UNITTESTS_REPORTING_MOCKS
#define HEMELB_UNITTESTS_REPORTING_MOCKS
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
              fake_time(0)
          {
          }
          ;
        protected:
          double CurrentTime()
          {
            fake_time += 10.0;
            return fake_time;
          }
        private:
          double fake_time;
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
            CPPUNIT_ASSERT_EQUAL(10, count);
            for (int i = 0; i < count; i++)
            {
              CPPUNIT_ASSERT_EQUAL(10.0 * i, sendbuf[i]);
              recvbuf[i] = 5.0 * i * calls;
            }
            calls++;
            return 0;
          }
        private:
          unsigned int calls;
      };
    }
  }
}
#endif // HEMELB_UNITTESTS_REPORTING_MOCKS
