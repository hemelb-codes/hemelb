#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <signal.h>
#include <netinet/in.h>

#include "steering/basic/ClientConnection.h"
#include "HttpPost.h"

namespace hemelb
{
  namespace steering
  {
    ClientConnection::ClientConnection(int iSteeringSessionId)
    {
      // Write the name of this machine to a file.
      {
        char thisMachineName[255];
        gethostname(thisMachineName, 255);
        FILE *f = fopen("env_details.asc", "w");
        fprintf(f, "%s\n", thisMachineName);
        fclose(f);
      }

      // Send the steering session id we're using to the rendezvous resource.
      {
        char steering_session_id_char[255];
        std::sprintf(steering_session_id_char, "%i", iSteeringSessionId);

        HttpPost::request("bunsen.chem.ucl.ac.uk", 28080, "/ahe/test/rendezvous/",
                          steering_session_id_char);
      }

      mCurrentSocket = -1;
      mIsBroken = false;
    }

    int ClientConnection::GetWorkingSocket()
    {
      // If the socket's broken, close it.
      if (mIsBroken)
      {
        close(mCurrentSocket);
      }

      // If we haven't yet had a socket, or the current one is broken, open a new one.
      if (mCurrentSocket < 0 || mIsBroken)
      {
        // Create the socket.
        int openSocket = socket(AF_INET, SOCK_STREAM, 0);
        if (openSocket == -1)
        {
          perror("socket");
          exit(1);
        }

        // Make the socket reusable.
        int yes = 1;
        if (setsockopt(openSocket, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1)
        {
          perror("setsockopt");
          exit(1);
        }

        // Bind to the socket.
        {
          struct sockaddr_in my_address;

          my_address.sin_family = AF_INET;
          my_address.sin_port = htons(MYPORT);
          my_address.sin_addr.s_addr = INADDR_ANY;
          memset(my_address.sin_zero, '\0', sizeof my_address.sin_zero);

          if (bind(openSocket, (struct sockaddr *) &my_address, sizeof my_address) == -1)
          {
            perror("bind");
            exit(1);
          }
        }

        // Mark the socket as accepting incoming connections.
        if (listen(openSocket, CONNECTION_BACKLOG) == -1)
        {
          perror("listen");
          exit(1);
        }

        // Accept an incoming connection from the client.
        struct sockaddr_in clientAddress;
        socklen_t socketSize = sizeof (clientAddress);

        mCurrentSocket = accept(openSocket, (struct sockaddr *) &clientAddress, &socketSize);

        if (mCurrentSocket == -1)
        {
          perror("accept");
        }

        // Close the open socket (we only want the client-specific one).
        close(openSocket);

        // We've only just created the socket so it shouldn't be broken.
        mIsBroken = false;
      }

      return mCurrentSocket;
    }

    void ClientConnection::ReportBroken()
    {
      mIsBroken = true;
    }

  }
}

