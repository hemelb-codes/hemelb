#ifndef HEMELB_REPORTING_POLICIES_H
#define HEMELB_REPORTING_POLICIES_H
#include <stdarg.h>
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace reporting
  {
    class FileWriterPolicy
    {
      public:
        FileWriterPolicy(const std::string &path)
        {
          file = fopen(path.c_str(), "w");
        }
        ~FileWriterPolicy()
        {
          fclose(file);
        }
      protected:
        void Print(const char * format, ...)
        {
          std::va_list arg;
          va_start(arg, format);
          vfprintf(file, format, arg);
          va_end(arg);
        }
      private:
        FILE* file;
    };

    class MPICommsPolicy
    {
      public:
        MPICommsPolicy() :
            instance(*hemelb::topology::NetworkTopology::Instance())
        {
        }
      protected:
        int Reduce(void *sendbuf,
                   void *recvbuf,
                   int count,
                   MPI_Datatype datatype,
                   MPI_Op op,
                   int root,
                   MPI_Comm comm)
        {
          return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
        }
        size_t FluidSitesOnProcessor(int n)
        {
          return instance.FluidSitesOnEachProcessor[n];
        }
        proc_t GetProcessorCount()
        {
          return instance.GetProcessorCount();
        }
        unsigned int GetMachineCount()
        {
          return instance.GetMachineCount();
        }
        int GetDepths()
        {
          return instance.GetDepths();
        }

      private:
        hemelb::topology::NetworkTopology& instance;
    };

    class HemeLBClockPolicy
    {
      protected:
        static double CurrentTime()
        {
          return hemelb::util::myClock();
        }
    };
  }
}
#endif // ONCE
