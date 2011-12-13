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
            CPPUNIT_ASSERT_EQUAL(11, count);
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

      class WriterMock
      {
        public:
          WriterMock(const std::string &path) :
              results(0)
          {
            CPPUNIT_ASSERT_EQUAL(std::string("mock_path"), path);
          }
          ~WriterMock()
          {
          }
          std::vector<std::string> & Results()
          {
           Stream();
            return results;
          }

        protected:
          void Print(const char * format, ...)
          {
            char buffer[1000];
            std::va_list arg;
            va_start(arg, format);
            vsprintf(buffer, format, arg);
            va_end(arg);
            Stream() << std::string(buffer) << std::endl;
          }

          std::ostream & Stream()
          {
            if (buffer.str() != "")
            {
              results.push_back(buffer.str());
              buffer.str("");
            }
            return buffer;
          }
        private:
          std::vector<std::string> results;
          std::stringstream buffer;
      };
    }
  }
}
#endif // HEMELB_UNITTESTS_REPORTING_MOCKS_H
