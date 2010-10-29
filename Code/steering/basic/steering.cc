#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#include <netdb.h>

#include <sys/time.h>
#include <sched.h>

#include <sys/stat.h>
#include <string.h>

#ifdef _AIX
#include <fcntl.h>
#endif

#include "constants.h"

#include "io/XdrWriter.h"
#include "io/XdrMemWriter.h"

#include "vis/ColourPalette.h"
#include "vis/visthread.h"

#include "steering/common/common.h"

#include "steering/basic/steering.h"
#include "steering/basic/Network.h"
#include "steering/basic/SimulationParameters.h"
#include "steering/basic/HttpPost.h"

#define MYPORT 65250
#define CONNECTION_BACKLOG 10

//pthread_mutex_t steer_param_lock = PTHREAD_MUTEX_INITIALIZER;

namespace hemelb
{
  namespace steering
  {

    double frameTiming()
    {
      struct timeval time_data;
      gettimeofday(&time_data, NULL);
      return (double) time_data.tv_sec + (double) time_data.tv_usec / 1.0e6;
    }

    void* hemeLB_steer(void* ptr)
    {
      int read_fd = *(int*) ptr;

      // printf("Kicking off steering thread with FD %i\n", read_fd);

      int num_chars = STEERABLE_PARAMETERS * sizeof(float) / sizeof(char);
      int bytes = sizeof(char) * num_chars;

      char* xdr_steering_data = new char[bytes];

      while (1)
      {
        XDR xdr_steering_stream;

        xdrmem_create(&xdr_steering_stream, xdr_steering_data, bytes,
                      XDR_DECODE);

        while (1)
        {
          struct timeval tv;
          fd_set readfds;
          int steerDataRecvB = 0;

          tv.tv_sec = 0;
          tv.tv_usec = 0;

          FD_ZERO(&readfds);
          FD_SET(read_fd, &readfds);

          select(read_fd + 1, &readfds, NULL, NULL, &tv);
          // printf("STEERING: Polling..\n"); fflush(0x0);

          if (FD_ISSET(read_fd, &readfds))
          {
            /* If there's something to read, read it... */
            //				printf("STEERING: Got data\n"); fflush(0x0);
            steerDataRecvB = Network::recv_all(read_fd, xdr_steering_data,
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
            delete[] xdr_steering_data;
            xdr_destroy(&xdr_steering_stream);
            return NULL;
          }
          sched_yield();
        }
        // pthread_mutex_lock(&steer_param_lock);

        sem_wait(&Control::Get()->steering_var_lock);

        for (int i = 0; i < STEERABLE_PARAMETERS; i++)
          xdr_float(&xdr_steering_stream, &steering::steer_par[i]);

        sem_post(&Control::Get()->steering_var_lock);

        /* printf("Got steering params ");
         for (int i = 0; i < STEERABLE_PARAMETERS; i++)
         printf("%0.4f ", steer_par[i]);
         printf("\n"); */

        if (steering::steer_par[14] > -1.0 && steering::steer_par[15] > -1.0)
          Control::Get()->updated_mouse_coords = 1;

        // pthread_mutex_unlock(&steer_param_lock);

        xdr_destroy(&xdr_steering_stream);
      }
      delete[] xdr_steering_data;

      return NULL;
    }
  }
}

