#ifndef HEMELB_STEERING_BASIC_CLIENTCONNECTION_H
#define HEMELB_STEERING_BASIC_CLIENTCONNECTION_H

#include <semaphore.h>

namespace hemelb
{
  namespace steering
  {
    class ClientConnection
    {
      public:
        ClientConnection(int iSteeringSessionId);
        ~ClientConnection();

        int GetWorkingSocket();

        void ReportBroken(int iSocketNum);

      private:
        static const unsigned int MYPORT = 65250;
        static const unsigned int CONNECTION_BACKLOG = 10;

        int mCurrentSocket;
        bool mIsBroken;
        // Use a semaphore to make sure that we don't create two new connections
        // when a broken one is reported simultaneously by two separate threads
        // (for example).
        sem_t mIsBusy;
    };
  }
}

#endif /* HEMELB_STEERING_BASIC_CLIENTCONNECTION_H */
