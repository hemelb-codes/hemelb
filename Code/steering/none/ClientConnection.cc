// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
    ClientConnection::ClientConnection(int iSteeringSessionId, reporting::Timers & timings) :
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

