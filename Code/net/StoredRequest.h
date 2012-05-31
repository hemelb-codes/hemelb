#ifndef HEMELB_NET_STOREDREQUEST_H
#define HEMELB_NET_STOREDREQUEST_H
#include <cstring>
#include "constants.h"
#include "mpiInclude.h"
namespace hemelb
{
  namespace net
  {

    class BaseRequest
    {
      public:
        void * Pointer;
        int Count;
        MPI_Datatype Type;
        proc_t Rank;
        BaseRequest(void *pointer, int count, MPI_Datatype type, proc_t rank) :
            Pointer(pointer), Count(count), Type(type), Rank(rank)
        {
        }
        virtual bool EnvelopeIdentical(const BaseRequest & other)
        {
          return ( (Count == other.Count) && (Rank == other.Rank) && (Type == other.Type));
        }
        virtual bool PayloadIdentical(const BaseRequest & other)
        {
          // reduction
          bool ok = true;
          for (unsigned int element = 0; element < Count; element++)
          {
            int size;
            MPI_Type_size(other.Type, &size);
            ok = ok && std::memcmp(other.Pointer + size * element, Pointer + size * element, size);
          }
          return ok;
        }
        virtual void Unpack(BaseRequest & other)
        {
          for (unsigned int element = 0; element < Count; element++)
          {
            int size;
            MPI_Type_size(other.Type, &size);
            std::memcpy(other.Pointer + size * element, Pointer + size * element, size);
          }
        }
    };

    class ScalarRequest : public BaseRequest
    {
      public:
        ScalarRequest(void *pointer, MPI_Datatype type, proc_t rank) :
            BaseRequest(pointer, 1, type, rank)
        {
        }
    };

    class GatherVReceiveRequest : public BaseRequest
    {
      public:
        int * Displacements;
        int * Counts;
        GatherVReceiveRequest(void *pointer, int *displacements, int *counts, MPI_Datatype type) :
            BaseRequest(pointer, 0, type, 0), Counts(counts)
        {
        }
    };
  }
}

#endif
