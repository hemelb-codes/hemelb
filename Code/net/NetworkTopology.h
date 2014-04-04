// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_NETWORKTOPOLOGY_H
#define HEMELB_NET_NETWORKTOPOLOGY_H

#include <vector>
#include <cstdio>

#include "constants.h"
#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace net
  {
    class NetworkTopology
    {
      public:
        static NetworkTopology* Instance();
        ~NetworkTopology();

        void Init(MpiCommunicator& commun);
        bool IsCurrentProcTheIOProc() const;
        int GetIOProcRank() const;

        // Functions for getting the rank of this processor and the total size
        // of the topology.
        proc_t GetLocalRank() const;
        proc_t GetProcessorCount() const;

        const MpiCommunicator & GetComms() const {
          return comms;
        }

      private:
        NetworkTopology();

        MpiCommunicator comms;

        /**
         * This variable is necessary, because the destructor for this static object will always
         * be called, regardless of whether the init method (that actually initialises the MPI
         * environment) is called.
         */
        static bool initialised;
        static NetworkTopology instance;
    };
  }
}

#endif /* HEMELB_NET_NETWORKTOPOLOGY_H */
