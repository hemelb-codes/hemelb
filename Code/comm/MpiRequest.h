//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MPIREQUEST_H
#define HEMELB_NET_MPIREQUEST_H

#include <vector>
#include "comm/Request.h"
#include "comm/MpiError.h"

namespace hemelb
{
  namespace comm
  {
    /**
     *
     */
    class MpiRequest : public Request
    {
      public:
        MpiRequest();
        MpiRequest(MPI_Request req);
        MpiRequest(MpiRequest&& req);
        MpiRequest& operator=(MpiRequest&& req);
        MpiRequest(const MpiRequest&) = delete;
        MpiRequest& operator=(const MpiRequest&) = delete;
      
        /**
         * Allow implicit casts to MPI_Request
         * @return The underlying MPI_Request
         */
        operator MPI_Request() const
        {
          return req;
        }
        operator bool() const;

        virtual void Wait();

        virtual bool Test();

      private:
        friend class MpiRequestList;
	MPI_Request req;
    };

    class MpiRequestList : public RequestList
    {
    public:
      virtual size_t size() const;
      virtual void resize(size_t i);
      virtual void push_back(Request::Ptr);
      virtual void set(size_t i, Request::Ptr);
      
      virtual void WaitAll();
      virtual bool TestAll();
    private:
      std::vector<MPI_Request> reqs;
    };
    
  }
}
#endif // HEMELB_NET_MPIREQUEST_H
