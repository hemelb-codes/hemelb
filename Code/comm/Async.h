// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_ASYNC_H
#define HEMELB_COMM_ASYNC_H

#include "comm/Communicator.h"
#include "comm/Request.h"

namespace hemelb
{
  namespace comm
  {
    // Simple class that manages a bunch of asynchronous
    // communications. Construct it and add the comms. Their
    // completion is waited on when the object is destructed.
    class Async
    {
    public:
      typedef std::shared_ptr<Async> Ptr;
      
      inline static Ptr New(Communicator::ConstPtr c)
      {
	return std::make_shared<Async>(c);
      }
      
      inline Async(Communicator::ConstPtr c) : comms(c), q(c->MakeRequestList())
      {
      }
      inline ~Async()
      {
	q->WaitAll();
      }

      inline void Wait() {
	q->WaitAll();
	q->clear();
      }
      
      inline Communicator::ConstPtr GetComm() const {
	return comms;
      }
      
      template <typename... Ts>
      void Isend(Ts... args)
      {
	q->push_back(comms->Isend(args...));
      }
      template <typename... Ts>
      void Issend(Ts... args)
      {
	q->push_back(comms->Issend(args...));
      }
      template <typename... Ts>
      void Irecv(Ts... args)
      {
	q->push_back(comms->Irecv(args...));
      }
      
    private:
      Communicator::ConstPtr comms;
      RequestList::Ptr q;
    };
  
  }
}

#endif
