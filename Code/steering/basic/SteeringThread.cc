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

    SteeringThread::SteeringThread(int fd, Control* controller) :
      mSteeringController(controller), mFdInt(fd)
    {
    }

    // Make it joinable
    pthread_attr_t* SteeringThread::GetPthreadAttributes()
    {
      pthread_attr_t* attr = new pthread_attr_t;

      pthread_attr_init(attr);
      pthread_attr_setdetachstate(attr, PTHREAD_CREATE_JOINABLE);
      return attr;
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
        io::XdrMemReader* steeringStream =
            new io::XdrMemReader(steeringRecvBuffer, bytes);

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
            steerDataRecvB = Network::recv_all(mFdInt, steeringRecvBuffer,
                                               num_chars);
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
        // pthread_mutex_lock(&steer_param_lock);

        sem_wait(&mSteeringController->steering_var_lock);

        for (int i = 0; i < STEERABLE_PARAMETERS; i++)
          //xdr_float(&xdr_steering_stream, &steering::steer_par[i]);
          steeringStream->readFloat(steering::steer_par[i]);

        sem_post(&mSteeringController->steering_var_lock);

        if (steering::steer_par[14] > -1.0 && steering::steer_par[15] > -1.0)
          mSteeringController->updated_mouse_coords = 1;

        delete steeringStream;
      }
      delete[] steeringRecvBuffer;

      return;

    }

  }
}
