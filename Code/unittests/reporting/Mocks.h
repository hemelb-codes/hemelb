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
          size_t FluidSitesOnProcessor(int n){
            return n*1000;
          }
          proc_t GetProcessorCount(){
            return 5;
          }
          unsigned int GetMachineCount(){
            return 4;
          }
          int GetDepths(){
            return 3;
          }
        private:
          unsigned int calls;
      };

      class WriterMock
      {
        public:
          WriterMock(const std::string &path):results(0)
          {
            CPPUNIT_ASSERT_EQUAL(std::string("mock_path"),path);
          }
          ~WriterMock()
          {
          }
          std::vector<std::string> & Results(){return results;}
        protected:
          void Print(const char * format, ...)
          {
            char buffer[1000];
            std::va_list arg;
            va_start(arg, format);
            vsprintf(buffer, format, arg);
            va_end(arg);
            results.push_back(std::string(buffer));
          }
        private:
           std::vector<std::string> results;
      };
    }
  }
}
#endif // HEMELB_UNITTESTS_REPORTING_MOCKS_H
