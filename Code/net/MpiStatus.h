//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MPISTATUS_H
#define HEMELB_NET_MPISTATUS_H

#include "net/MpiError.h"
#include <memory>

namespace hemelb
{
  namespace net
  {
    class MpiRequest;

    class MpiStatus
    {
      public:
        MpiStatus();
        MpiStatus(MPI_Status stat);

        /**
         * Allow implicit casts to MPI_Status
         * @return The underlying MPI_Status
         */
        operator MPI_Status() const
        {
          return *statPtr;
        }

      private:
        std::shared_ptr<MPI_Status> statPtr;
        // Request needs to be able to access statPtr.
        friend class MpiRequest;
    };
  }
}
#endif // HEMELB_NET_MPISTATUS_H
