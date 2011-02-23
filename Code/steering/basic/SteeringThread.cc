#include <pthread.h>
#include <unistd.h>

#include "io/XdrMemReader.h"

#include "steering/basic/SteeringThread.h"
#include "steering/basic/Network.h"
#include "steering/common/common.h"

namespace hemelb
{
  namespace steering
  {

    SteeringThread::SteeringThread(ClientConnection* iClientConn, sem_t* bSteeringVariableLock) :
      mClientConn(iClientConn), mSteeringVariableLock(bSteeringVariableLock)
    {
    }

    // The work routine of the thread - reads data from the socket to the client.
    // This is original code with minimal tweaks to make it work with
    // the new (Oct 2010) structure.
    void SteeringThread::DoWork(void)
    {
      // Create a buffer for the data received.
      const int num_chars = STEERABLE_PARAMETERS * sizeof(float) / sizeof(char);
      const int bytes = sizeof(char) * num_chars;

      char steeringRecvBuffer[bytes];

      // This outer loop is infinite - we keep receiving frame after frame.
      while (true)
      {
        // Initialise the stream reader.
        io::XdrMemReader steeringStream(steeringRecvBuffer, bytes);

        bool readOK = false;

        // Keep looping until we receive a frame.
        while (true)
        {
          int mSocketFileDescriptor = mClientConn->GetWorkingSocket();

          // Add the socket to a list to be checked for read-readiness.
          fd_set readableDescriptors;
          {
            struct timeval tv;
            tv.tv_sec = 0;
            tv.tv_usec = 0;

            FD_ZERO(&readableDescriptors);
            FD_SET(mSocketFileDescriptor, &readableDescriptors);

            select(mSocketFileDescriptor + 1, &readableDescriptors, NULL, NULL, &tv);
          }

          // If the socket is readable, read from it.
          if (FD_ISSET(mSocketFileDescriptor, &readableDescriptors))
          {
            // Try to receive all the data.
            int steerDataRecvB = Network::recv_all(mSocketFileDescriptor, steeringRecvBuffer,
                                                   num_chars);

            // If there was an error, report it and return.
            if (steerDataRecvB < 0)
            {
              printf("Steering thread: broken network pipe...\n");
              mClientConn->ReportBroken(mSocketFileDescriptor);
              break;
            }

            readOK = true;

            // Otherwise yield the thread, and break from the loop as we've received the data we need.
            sched_yield();
            break;
          }
          // If the socket isn't readable, sleep for a bit before trying again.
          else
          {
            usleep(5000);
            sched_yield();
          }
        }

        // We've got a frame. Lock the steering variables, and
        // write to them.
        if (readOK)
        {
          sem_wait(mSteeringVariableLock);

          for (int i = 0; i < STEERABLE_PARAMETERS; i++)
            steeringStream.readFloat(steering::steer_par[i]);

          sem_post(mSteeringVariableLock);
        }
      }
    }

  }
}
