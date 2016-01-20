
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_IOCOMMUNICATOR_H
#define HEMELB_NET_IOCOMMUNICATOR_H

//#include <vector>
//#include <cstdio>
//
//#include "constants.h"
#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace net
  {
    /**
     * An MPI communicator with a special I/O rank.
     */
    class IOCommunicator : public MpiCommunicator
    {
      public:
        IOCommunicator(const MpiCommunicator& comm);
        bool OnIORank() const;
        int GetIORank() const;
    };
  }
}

#endif /* HEMELB_NET_IOCOMMUNICATOR_H */
