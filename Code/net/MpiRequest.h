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
#include "net/MpiError.h"
#include <boost/shared_ptr.hpp>

namespace hemelb
{
  namespace net
  {
    class MpiStatus;
    /**
     *
     */
    class MpiRequest
    {
      private:
        typedef std::vector<MpiRequest> ReqVec;
        typedef std::vector<MpiStatus> StatVec;

      public:
        MpiRequest();
        MpiRequest(MPI_Request req);

        /**
         * Allow implicit casts to MPI_Request
         * @return The underlying MPI_Request
         */
        operator MPI_Request() const
        {
          return *reqPtr;
        }

        operator bool() const;

        void Wait();
        void Wait(MpiStatus& stat);

        static void WaitAll(ReqVec& reqs);
        static void WaitAll(ReqVec& reqs, StatVec& stats);

        bool Test();
        static bool TestAll(ReqVec& reqs);

      private:

        boost::shared_ptr<MPI_Request> reqPtr;
    };

  }
}
#endif // HEMELB_NET_MPIREQUEST_H
