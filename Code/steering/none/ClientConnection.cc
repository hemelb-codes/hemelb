#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <csignal>
#include <netinet/in.h>

#include "steering/ClientConnection.h"

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
    ClientConnection::ClientConnection(int iSteeringSessionId, reporting::Timers & timings):
        mIsBusy(), timers(timings)
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

