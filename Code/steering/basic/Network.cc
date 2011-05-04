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
    ssize_t Network::recv_all(int sockid, char *buf, const int length)
    {
      ssize_t received_bytes = 0;
      ssize_t bytes_left_to_receive = length;
      ssize_t n = 0;

      // TODO: Make this better.
      // While some data left to be received...
      while (received_bytes < length)
      {
        // Receive some data (up to the remaining length)
        n = recv(sockid, buf + received_bytes, bytes_left_to_receive, 0);

        if (n <= 0)
        {
          // Distinguish between cases where the pipe fails because it'd block
          // (No problem, we'll try again later) or because the pipe is broken.
          return n;
        }
        else
        {
          received_bytes += n;
          bytes_left_to_receive -= n;
        }
      }

      return received_bytes;
    }

    /**
     * Send all bytes from a buffer of known length over a socket.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns the number of bytes sent or -1 on failure.
     */
    ssize_t Network::send_all(int sockid, const char *buf, const int length)
    {
      ssize_t sent_bytes = 0;
      ssize_t bytes_left_to_send = length;
      ssize_t n = 0;

      // TODO: Make this better.
      while (sent_bytes < length)
      {
        n = send(sockid, buf + sent_bytes, bytes_left_to_send, 0);

        if (n <= 0)
        {
          // Distinguish between cases where the pipe fails because it'd block
          // (No problem, we'll try again later) or because the pipe is broken.
          if (errno == EWOULDBLOCK)
          {
            // Wait a bit before trying again.
            usleep(10);
          }
          else
          {
            return n;
          }
        }
        else
        {
          sent_bytes += n;
          bytes_left_to_send -= n;
        }
      }
      return sent_bytes;
    }

  }
}
