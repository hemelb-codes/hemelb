#ifndef HEMELB_STEERING_NETWORK_H
#define HEMELB_STEERING_NETWORK_H

namespace hemelb
{
  namespace steering
  {

    class Network
    {
      public:
        // Receive a bytestream of known length from a socket into a buffer.
        static int recv_all(int sockid, char *buf, const int length);

        // Send all bytes from a buffer of known length over a socket.
        static int send_all(int sockid, const char *buf, const int length);
    };

  }
}

#endif // HEMELB_STEERING_NETWORK_H
