// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
        /**
         * Get the singleton instance.
         * @return
         */
        static IOCommunicator* Instance();
        /**
         * Initalise the singleton instance.
         * @param commun
         */
        static void Init(MpiCommunicator& commun);

        bool OnIORank() const;
        int GetIORank() const;

      private:
        IOCommunicator(const MpiCommunicator& comm);
        IOCommunicator();

        /**
         * This variable is necessary, because the destructor for this static object will always
         * be called, regardless of whether the init method (that actually initialises the MPI
         * environment) is called.
         */
        static bool initialised;
        static IOCommunicator instance;
    };
  }
}

#endif /* HEMELB_NET_IOCOMMUNICATOR_H */
