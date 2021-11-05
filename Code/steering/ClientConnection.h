// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_STEERING_CLIENTCONNECTION_H
#define HEMELB_STEERING_CLIENTCONNECTION_H

#include <netinet/in.h>
#include <mutex>
#include "reporting/Timers.h"

namespace hemelb
{
  namespace steering
  {
    class ClientConnection
    {
      public:
        ClientConnection(int iSteeringSessionId, reporting::Timers & timings);
        ~ClientConnection();

        int GetWorkingSocket();

        void ReportBroken(int iSocketNum);

      private:
        static const in_port_t MYPORT = 65250;
        static const unsigned int CONNECTION_BACKLOG = 10;

        int mCurrentSocket;
        int mListeningSocket;
        bool mIsBroken;
        // Use a semaphore to make sure that we don't create two new connections
        // when a broken one is reported simultaneously by two separate threads
        // (for example).
        std::mutex mIsBusy;
        reporting::Timers & timers;

    };
  }
}

#endif /* HEMELB_STEERING_CLIENTCONNECTION_H */
