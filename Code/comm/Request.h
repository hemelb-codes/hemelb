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
    // class Status;
    /**
     *
     */
    class Request
    {
      public:
        typedef std::shared_ptr<Request> Ptr;
        typedef std::vector<Ptr> ReqVec;
      //typedef std::vector<Status> StatVec;

        Request();
      //Request(MPI_Request req);

        /**
         * Allow implicit casts to MPI_Request
         * @return The underlying MPI_Request
         */
        // operator MPI_Request() const
        // {
        //   return req;
        // }
        operator bool() const;

        virtual void Wait() = 0;
        // void Wait(Status& stat);

        static void WaitAll(ReqVec& reqs);
        // static void WaitAll(ReqVec& reqs, StatVec& stats);

        virtual bool Test() = 0;
        static bool TestAll(ReqVec& reqs);

    };

  }
}
#endif // HEMELB_COMM_REQUEST_H
