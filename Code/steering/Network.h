#ifndef HEMELB_STEERING_NETWORK_H
#define HEMELB_STEERING_NETWORK_H

#include <string>

#include "steering/ClientConnection.h"
#include "net/IteratedAction.h"

namespace hemelb
{
  namespace steering
  {

    class Network : public net::IteratedAction
    {
      public:
        Network(int iSteeringSessionId, reporting::Timers & timings);

        // Receive a bytestream of known length from a socket into a buffer.
        bool recv_all(char *buf, const int length);

        // Send all bytes from a buffer of known length over a socket.
        bool send_all(const char *buf, const int length);

        /**
         * Use the time between MPI send and MPI receives being completed to perform the
         * Network sends.
         */
        void PreReceive();

        bool IsConnected();

      private:
        void Break(int socket);

        long sendInternal(const char* data, long length, int socket);

        ClientConnection clientConnection;

        // Buffers to keep the data from partial sends and receives.
        std::string sendBuf;
        std::string recvBuf;
    };

  }
}

#endif // HEMELB_STEERING_NETWORK_H
