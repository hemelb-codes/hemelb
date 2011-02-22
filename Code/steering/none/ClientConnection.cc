#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <signal.h>
#include <netinet/in.h>

#include "steering/ClientConnection.h"
#include "HttpPost.h"

namespace hemelb
{
  namespace steering
  {

    /**
     * In the 'no-steering' code, this should do nothing.
     *
     * @param iSteeringSessionId
     * @return
     */
    ClientConnection::ClientConnection(int iSteeringSessionId)
    {
    }

    ClientConnection::~ClientConnection()
    {
    }

    int ClientConnection::GetWorkingSocket()
    {
      return -1;
    }

    void ClientConnection::ReportBroken(int iSocketNum)
    {
    }

  }
}

