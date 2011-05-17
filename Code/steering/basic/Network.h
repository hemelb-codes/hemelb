#ifndef HEMELB_STEERING_NETWORK_H
#define HEMELB_STEERING_NETWORK_H

#include "steering/ClientConnection.h"
#include "net/IteratedAction.h"

namespace hemelb
{
  namespace steering
  {

    class Network : public net::IteratedAction
    {
      public:
        Network(int iSteeringSessionId);

        // Receive a bytestream of known length from a socket into a buffer.
        bool recv_all(char *buf, const int length);

        // Send all bytes from a buffer of known length over a socket.
        bool send_all(const char *buf, const int length);

        /**
         * Use the time between MPI send and MPI receives being completed to perform the
         * Network sends and receives.
         */
        void PreReceive();

        bool IsConnected();

      private:
        ClientConnection clientConnection;

        char* sendBuf;
        char* recvBuf;

        int sendBufLen;
        int sendBufUsed;

        int recvBufLen;
        int recvBufUsed;
    };

  }
}

#endif // HEMELB_STEERING_NETWORK_H
