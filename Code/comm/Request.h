//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_COMM_REQUEST_H
#define HEMELB_COMM_REQUEST_H

#include <vector>

namespace hemelb
{
  namespace comm
  {
    class RequestList;
    /**
     *
     */
    class Request
    {
      public:
        typedef std::shared_ptr<Request> Ptr;

        Request() = default;
        Request(Request&& req) = default;
        Request& operator=(Request&& req) = default;
        Request(const Request&) = delete;
        Request& operator=(const Request&) = delete;
        operator bool() const;

        virtual void Wait() = 0;
        virtual bool Test() = 0;
    };

    class RequestList
    {
    public:
      typedef std::shared_ptr<RequestList> Ptr;
      
      virtual size_t size() const = 0;
      virtual void resize(size_t i) = 0;
      virtual void push_back(Request&&) = 0;
      virtual void set(size_t i, Request&&) = 0;
      
      virtual void WaitAll() = 0;
      virtual bool TestAll() = 0;
    };
  }
}
#endif // HEMELB_COMM_REQUEST_H
