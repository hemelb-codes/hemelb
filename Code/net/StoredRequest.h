// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_STOREDREQUEST_H
#define HEMELB_NET_STOREDREQUEST_H
#include <cstring>
#include "constants.h"
#include "mpiInclude.h"
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
