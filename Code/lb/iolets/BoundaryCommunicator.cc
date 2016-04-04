
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
