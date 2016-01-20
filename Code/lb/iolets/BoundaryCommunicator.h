
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
