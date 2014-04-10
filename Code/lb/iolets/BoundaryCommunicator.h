//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_IOLETS_BOUNDARYCOMMUNICATOR_H
#define HEMELB_LB_IOLETS_BOUNDARYCOMMUNICATOR_H

#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class BoundaryCommunicator : public net::MpiCommunicator
      {
        public:
          BoundaryCommunicator(const net::MpiCommunicator& parent);
          bool IsCurrentProcTheBCProc() const;
          int GetBCProcRank() const;
      };
    }
  }
}

#endif // HEMELB_LB_IOLETS_BOUNDARYCOMMUNICATOR_H
