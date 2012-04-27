#ifndef HEMELB_STEERING_CLIENTCONNECTION_H
#define HEMELB_STEERING_CLIENTCONNECTION_H

#include <netinet/in.h>
#include <semaphore.h>
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
        reporting::Timers & timers;
        // Use a semaphore to make sure that we don't create two new connections
        // when a broken one is reported simultaneously by two separate threads
        // (for example).
        sem_t mIsBusy;
    };
  }
}

#endif /* HEMELB_STEERING_CLIENTCONNECTION_H */
