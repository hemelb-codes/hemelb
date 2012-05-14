#include <unistd.h>
#include <cerrno>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <csignal>
#include <sys/stat.h>
#include <cstdio>

#include "debug/Debugger.h"
#include "log/Logger.h"
#include "steering/Network.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace steering
  {
    Network::Network(int steeringSessionId, reporting::Timers & timings) :
      clientConnection(steeringSessionId, timings)
    {

    }

    /**
     * Do nothing.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns true if we have successfully provided that much data.
     */
    bool Network::recv_all(char *buf, const int length)
    {
      return false;
    }

    void Network::PreReceive()
    {
    }

    bool Network::IsConnected()
    {
      return false;
    }

    /**
     * Do nothing.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns the number of bytes sent or -1 on failure.
     */
    bool Network::send_all(const char *buf, const int length)
    {
      return false;
    }

    void Network::Break(int socket)
    {
    }

    /**
     * Do nothing.
     *
     * @param data
     * @param length
     * @param socket
     * @return
     */
    long Network::sendInternal(const char* data, long length, int socket)
    {
      return -1;
    }

  }
}
