#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <sys/stat.h>
#include <stdio.h>

#include "steering/basic/Network.h"

namespace hemelb
{
  namespace steering
  {
    /**
     * Receive a bytestream of known length from a socket into a buffer.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns the number of bytes recieved or -1 on failure.
     */
    int Network::recv_all(int sockid, char *buf, const int length)
    {
      int received_bytes = 0;
      int bytes_left_to_receive = length;
      int n = 0;

      // While some data left to be received...
      while (received_bytes < length)
      {
        // Receive some data (up to the remaining length)
        n = recv(sockid, buf + received_bytes, bytes_left_to_receive, NULL);

        if (n == -1)
          break;

        received_bytes += n;
        bytes_left_to_receive -= n;
      }

      return n == -1
        ? -1
        : received_bytes;
    }

    /**
     * Send all bytes from a buffer of known length over a socket.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns the number of bytes sent or -1 on failure.
     */
    int Network::send_all(int sockid, const char *buf, const int length)
    {
      int sent_bytes = 0;
      int bytes_left_to_send = length;
      int n = 0;

      while (sent_bytes < length)
      {
        n = send(sockid, buf + sent_bytes, bytes_left_to_send, 0);

        if (n == -1)
          break;

        sent_bytes += n;
        bytes_left_to_send -= n;
      }

      return n == -1
        ? -1
        : sent_bytes;
    }

  }
}
