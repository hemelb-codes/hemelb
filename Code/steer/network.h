#ifndef __steer_network_h_
#define __steer_network_h_

#ifndef NO_STEER

class Network
{
  public:
    // Receive a bytestream of known length from a socket into a buffer.
    static int recv_all (int sockid, char *buf, int *length);
    // Send all bytes from a buffer of known length over a socket.
    static int send_all (int sockid, char *buf, int *length);
};

#endif// NO_STEER

#endif//__steer_network_h_
