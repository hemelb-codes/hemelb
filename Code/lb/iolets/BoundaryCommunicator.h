
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_BOUNDARYCOMMUNICATOR_H
#define HEMELB_LB_IOLETS_BOUNDARYCOMMUNICATOR_H

#include "comm/Communicator.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class BoundaryCommunicator
      {
      public:
	BoundaryCommunicator(comm::Communicator::ConstPtr parent);
	bool IsCurrentProcTheBCProc() const;
	int GetBCProcRank() const;
	comm::Communicator::ConstPtr GetComm() const;
      private:
	comm::Communicator::Ptr comm;
	
      };
    }
  }
}

#endif // HEMELB_LB_IOLETS_BOUNDARYCOMMUNICATOR_H
