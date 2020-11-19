// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <csignal>
#include <fcntl.h>
#include <netinet/in.h>

#include "log/Logger.h"
#include "steering/ClientConnection.h"
#include "steering/basic/HttpPost.h"

namespace hemelb
{
  namespace steering
  {
    ClientConnection::ClientConnection(int iSteeringSessionId, reporting::Timers & timings) :
        mIsBusy(), timers(timings)
    {
      // Write the name of this machine to a file.

      {
        char thisMachineName[255];
        gethostname(thisMachineName, 255);
        std::FILE *f = std::fopen("env_details.asc", "w");
        std::fprintf(f, "%s\n", thisMachineName);
        std::fclose(f);
      }

      mCurrentSocket = -1;
      mIsBroken = false;

      // Create the socket.
      mListeningSocket = socket(AF_INET, SOCK_STREAM, 0);
      if (mListeningSocket == -1)
      {
        perror("socket");
        exit(1);
      }

      // Make the socket reusable.
      int yes = 1;
      if (setsockopt(mListeningSocket, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1)
      {
        perror("setsockopt");
        exit(1);
      }

      // Bind to the socket.
      {
        struct sockaddr_in my_address;

        my_address.sin_family = AF_INET;
        my_address.sin_port = htons((in_port_t ) MYPORT);
        my_address.sin_addr.s_addr = INADDR_ANY;
        memset(my_address.sin_zero, '\0', sizeof my_address.sin_zero);

        if (bind(mListeningSocket, (struct sockaddr *) &my_address, sizeof my_address) == -1)
        {
          perror("bind");
          exit(1);
        }
      }

      // Mark the socket as accepting incoming connections.
      if (listen(mListeningSocket, CONNECTION_BACKLOG) == -1)
      {
        perror("listen");
        exit(1);
      }
    }

    ClientConnection::~ClientConnection()
    {
      close(mCurrentSocket);
      close(mListeningSocket);
    }

    int ClientConnection::GetWorkingSocket()
    {
      // Lock the mutex and release it when this goes out of scope
      std::lock_guard<std::mutex> lock(mIsBusy);

      // If we haven't yet had a socket, or the current one is broken, open a new one.
      if (mCurrentSocket < 0 || mIsBroken)
      {
        // Accept an incoming connection from the client.
        struct sockaddr_in clientAddress;
        socklen_t socketSize = sizeof (clientAddress);

        int lOldSocket = mCurrentSocket;

        // Make the socket non-blocking, just while we try to accept on it, then set
        // it back to what it was before.
        {
          int flags = fcntl(mListeningSocket, F_GETFL, 0);
          if (flags == -1)
          {
            flags = 0;
          }
#ifndef HEMELB_WAIT_ON_CONNECT
          if (fcntl(mListeningSocket, F_SETFL, flags | O_NONBLOCK) < 0)
          {
            perror("flags");
          }
#else
          log::Logger::Log<log::Info, log::Singleton>("Waiting for steering client connection");
          timers[reporting::Timers::steeringWait].Start();

#endif
          // Try to accept a socket (from the non-blocking socket)
          mCurrentSocket = accept(mListeningSocket,
                                  (struct sockaddr *) &clientAddress,
                                  &socketSize);
#ifdef HEMELB_WAIT_ON_CONNECT
          timers[reporting::Timers::steeringWait].Stop();
          log::Logger::Log<log::Debug, log::Singleton>("Continuing after receiving steering connection.");
#endif
          // We've got a socket - make that socket non-blocking too.
          if (mCurrentSocket > 0)
          {
            log::Logger::Log<log::Info, log::Singleton>("Steering client connected");
            flags = fcntl(mCurrentSocket, F_GETFL, 0);
            if (flags == -1)
            {
              flags = 0;
            }
            if (fcntl(mCurrentSocket, F_SETFL, flags | O_NONBLOCK) < 0)
            {
              perror("flags");
            }
          }
        }

        // If we had a socket before, close it.
        if (lOldSocket > 0)
        {
          close(lOldSocket);
        }

        // We've only just created the socket so it shouldn't be broken.
        mIsBroken = false;
      }

      int lRet = mCurrentSocket;

      return lRet;
    }

    void ClientConnection::ReportBroken(int iSocketNum)
    {
      // Lock the mutex and release it when this goes out of scope
      std::lock_guard<std::mutex> lock(mIsBusy);
      if (mCurrentSocket == iSocketNum)
      {
        mIsBroken = true;
      }
    }

  }
}
