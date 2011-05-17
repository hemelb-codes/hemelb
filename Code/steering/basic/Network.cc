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

#include "log/Logger.h"
#include "steering/basic/Network.h"

namespace hemelb
{
  namespace steering
  {
    Network::Network(int steeringSessionId) :
      clientConnection(steeringSessionId)
    {

    }

    /**
     * Receive a bytestream of known length from a socket into a buffer.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns the number of bytes recieved or -1 on failure.
     */
    bool Network::recv_all(char *buf, const int length)
    {
      int socketToClient = clientConnection.GetWorkingSocket();

      if (socketToClient < 0)
      {
        return false;
      }

      ssize_t received_bytes = 0;
      ssize_t bytes_left_to_receive = length;
      ssize_t n = 0;

      // TODO: Make this better.
      // While some data left to be received...
      while (received_bytes < length)
      {
        // Receive some data (up to the remaining length)
        n = recv(socketToClient, buf + received_bytes, bytes_left_to_receive, 0);

        // If there was an error, report it and return.
        if (n <= 0)
        {
          // If there was no data and it wasn't simply that the socket would block,
          // raise an error.
          if (errno != EAGAIN)
          {
            log::Logger::Log<log::Warning, log::Singleton>("Steering component: broken network pipe... (%s)",
                                                           strerror(errno));
            clientConnection.ReportBroken(socketToClient);
          }

          return false;
        }
        else
        {
          received_bytes += n;
          bytes_left_to_receive -= n;
        }
      }

      return true;
    }

    void Network::PreReceive()
    {

    }

    bool Network::IsConnected()
    {
      return clientConnection.GetWorkingSocket() > 0;
    }

    /**
     * Send all bytes from a buffer of known length over a socket.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns the number of bytes sent or -1 on failure.
     */
    bool Network::send_all(const char *buf, const int length)
    {
      // Get a socket.
      int socketToClient = clientConnection.GetWorkingSocket();

      // If there's no such socket, we don't have a connection. Nothing to do.
      if (socketToClient < 0)
      {
        return false;
      }

      ssize_t sent_bytes = 0;
      ssize_t bytes_left_to_send = length;
      ssize_t n = 0;

      // TODO: Make this better.
      while (sent_bytes < length)
      {
        n = send(socketToClient, buf + sent_bytes, bytes_left_to_send, 0);

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
            log::Logger::Log<log::Warning, log::Singleton>("Network send had broken pipe... (%s)",
                                                           strerror(errno));
            clientConnection.ReportBroken(socketToClient);

            return false;
          }
        }
        else
        {
          sent_bytes += n;
          bytes_left_to_send -= n;
        }
      }

      return true;
    }

  }
}
