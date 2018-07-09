
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

    class SimpleRequest
    {
      public:
        void * Pointer;
        int Count;
        MPI_Datatype Type;
        proc_t Rank;
        SimpleRequest(void *pointer, int count, MPI_Datatype type, proc_t rank) :
            Pointer(pointer), Count(count), Type(type), Rank(rank)
        {
        }
    };

    class ScalarRequest : public SimpleRequest
    {
      public:
        ScalarRequest(void *pointer, MPI_Datatype type, proc_t rank) :
          SimpleRequest(pointer, 1, type, rank)
        {
        }
    };

    class GatherVReceiveRequest : public SimpleRequest
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
