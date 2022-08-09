// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_STOREDREQUEST_H
#define HEMELB_NET_STOREDREQUEST_H
#include <cstring>
#include "constants.h"
#include "net/mpi.h"
namespace hemelb
{
  namespace net
  {

    template <bool ptr_const>
    class SimpleRequest
    {
      public:
        static constexpr bool is_ptr_const = ptr_const;
        using ptr = std::conditional_t<ptr_const, void const*, void*>;
        ptr Pointer;
        int Count;
        MPI_Datatype Type;
        proc_t Rank;
        SimpleRequest(ptr pointer, int count, MPI_Datatype type, proc_t rank) :
            Pointer(pointer), Count(count), Type(type), Rank(rank)
        {
        }
    };
    SimpleRequest(void*, int, MPI_Datatype, proc_t) -> SimpleRequest<false>;
    SimpleRequest(void const*, int, MPI_Datatype, proc_t) -> SimpleRequest<true>;

    template <bool ptr_const>
    class ScalarRequest : public SimpleRequest<ptr_const>
    {
      public:
        using base = SimpleRequest<ptr_const>;
        //using SimpleRequest<ptr_const>::ptr;
        ScalarRequest(typename base::ptr pointer, MPI_Datatype type, proc_t rank) :
            SimpleRequest<ptr_const>(pointer, 1, type, rank)
        {
        }
    };
    ScalarRequest(void*, MPI_Datatype, proc_t) -> ScalarRequest<false>;
    ScalarRequest(void const*, MPI_Datatype, proc_t) -> ScalarRequest<true>;

    class GatherVReceiveRequest : public SimpleRequest<false>
    {
      public:
        int * Counts;
        int * Displacements;
        GatherVReceiveRequest(void *pointer, int *displacements, int *counts, MPI_Datatype type) :
            SimpleRequest(pointer, 0, type, 0), Counts(counts), Displacements(displacements)
        {
        }
    };
  }
}

#endif
