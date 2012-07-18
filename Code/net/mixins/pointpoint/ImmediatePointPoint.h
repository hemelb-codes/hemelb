// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_MIXINS_POINTPOINT_IMMEDIATEPOINTPOINT_H
#define HEMELB_NET_MIXINS_POINTPOINT_IMMEDIATEPOINTPOINT_H
#include "net/BaseNet.h"
#include "net/mixins/StoringNet.h"
namespace hemelb
{
  namespace net
  {
    // although ImmediatePointPoint does not use the StoringNet capabilities at all
    // it needs to inherit it
    // so that it becomes the unique final overrider
    class ImmediatePointPoint : public virtual StoringNet
    {

      public:
        ImmediatePointPoint()
        {
        }
        ~ImmediatePointPoint();

        void WaitPointToPoint();
        // we will *NOT* store the requests, so we must provide RequestSendImpl ourselves.
        virtual void RequestSendImpl(void* pointer, int count, proc_t rank, MPI_Datatype type);
        virtual void RequestReceiveImpl(void* pointer, int count, proc_t rank, MPI_Datatype type);
      protected:
        void ReceivePointToPoint(); //PASS
        void SendPointToPoint(); //PASS

    };
  }
}

#endif
