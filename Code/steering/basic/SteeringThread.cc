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

    SteeringThread::SteeringThread(int fd, sem_t* bSteeringVariableLock) :
      mSteeringVariableLock(bSteeringVariableLock), mFdInt(fd)
    {
    }

    // The work routine of the thread.
    // This is original code with minimal tweaks to make it work with
    // the new (Oct 2010) structure.
    void SteeringThread::DoWork(void)
    {
      // printf("Kicking off steering thread with FD %i\n", read_fd);

      int num_chars = STEERABLE_PARAMETERS * sizeof(float) / sizeof(char);
      int bytes = sizeof(char) * num_chars;

      char* steeringRecvBuffer = new char[bytes];

      while (1)
      {
        //XDR xdr_steering_stream;

        //xdrmem_create(&xdr_steering_stream, xdr_steering_data, bytes,
        //              XDR_DECODE);
        io::XdrMemReader* steeringStream = new io::XdrMemReader(steeringRecvBuffer, bytes);

        while (1)
        {
          struct timeval tv;
          fd_set readfds;
          int steerDataRecvB = 0;

          tv.tv_sec = 0;
          tv.tv_usec = 0;

          FD_ZERO(&readfds);
          FD_SET(mFdInt, &readfds);

          select(mFdInt + 1, &readfds, NULL, NULL, &tv);
          // printf("STEERING: Polling..\n"); fflush(0x0);

          if (FD_ISSET(mFdInt, &readfds))
          {
            /* If there's something to read, read it... */
            //        printf("STEERING: Got data\n"); fflush(0x0);
            steerDataRecvB = Network::recv_all(mFdInt, steeringRecvBuffer, num_chars);
            sched_yield();
            break;
          }
          else
          {
            usleep(5000);
          }

          if (steerDataRecvB < 0)
          {
            printf("Steering thread: broken network pipe...\n");
            delete steeringStream;
            delete[] steeringRecvBuffer;
            return;
          }
          sched_yield();
        }

        sem_wait(mSteeringVariableLock);

        for (int i = 0; i < STEERABLE_PARAMETERS; i++)
          steeringStream->readFloat(steering::steer_par[i]);

        sem_post(mSteeringVariableLock);

        delete steeringStream;
      }
      delete[] steeringRecvBuffer;

      return;

    }

  }
}
