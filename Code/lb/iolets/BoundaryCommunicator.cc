//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "lb/iolets/BoundaryCommunicator.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      BoundaryCommunicator::BoundaryCommunicator(const net::MpiCommunicator& parent) :
          MpiCommunicator(parent.Duplicate())
      {

      }
      bool BoundaryCommunicator::IsCurrentProcTheBCProc() const
      {
        return Rank() == GetBCProcRank();
      }
      int BoundaryCommunicator::GetBCProcRank() const
      {
        return 0;
      }
    }
  }
}
