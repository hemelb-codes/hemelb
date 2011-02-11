#include <pthread.h>
#include <cstdio>
#include <cstring>
#include <netinet/in.h>
#include <sys/time.h>
#include <unistd.h>
#include <signal.h>

#include <sys/types.h>

#include "io/XdrMemWriter.h"

#include "steering/basic/NetworkThread.h"
#include "steering/basic/SteeringThread.h"
#include "steering/basic/Control.h"
#include "steering/basic/HttpPost.h"
#include "steering/basic/Network.h"
#include "steering/basic/SimulationParameters.h"

namespace hemelb
{
  namespace steering
  {
    pthread_mutex_t NetworkThread::var_lock = PTHREAD_MUTEX_INITIALIZER;

    // Use initialisation list to do the work.
    NetworkThread::NetworkThread(lb::LBM* lbm,
                                 Control* steeringController,
                                 lb::SimulationState* iSimState,
                                 const lb::LbmParameters* iLbmParams) :
      mLbm(lbm), mSteeringController(steeringController), mSimState(iSimState),
          mLbmParams(iLbmParams)
    {
      /* Storing references to lbm and steeringController; we
       * won't want to destroy them.
       */

      xdrSendBuffer_pixel_data = new char[pixel_data_bytes];
      xdrSendBuffer_frame_details = new char[frame_details_bytes];
    }

    NetworkThread::~NetworkThread()
    {
      delete[] xdrSendBuffer_frame_details;
      delete[] xdrSendBuffer_pixel_data;
    }

    // Override the base class to make sure this thread's joinable.
    pthread_attr_t* NetworkThread::GetPthreadAttributes(void)
    {
      pthread_attr_t* attr = new pthread_attr_t;
      pthread_attr_init(attr);
      pthread_attr_setdetachstate(attr, PTHREAD_CREATE_JOINABLE);
      return attr;
    }

    void NetworkThread::setRenderState(int val)
    {
      pthread_mutex_lock(&var_lock);
      mSimState->DoRendering = val;
      pthread_mutex_unlock(&var_lock);
    }

    // Return seconds since epoch to microsec precision.
    double NetworkThread::frameTiming()
    {
      struct timeval time_data;
      gettimeofday(&time_data, NULL);
      return (double) time_data.tv_sec + (double) time_data.tv_usec / 1.0e6;
    }

    // The work routine of the thread.
    // This is original code with minimal tweaks to make it work with
    // the new (Oct 2010) structure.
    void NetworkThread::DoWork(void)
    {
      char steering_session_id_char[255];

      // PG compiler sticks to the letter of the standard: snprintf should not
      // be exported in std:: namespace as GCC does.
      // std::snprintf(steering_session_id_char, 255, "%i",
      // mLbm->steering_session_id);

      std::sprintf(steering_session_id_char, "%i", mLbm->steering_session_id);

      setRenderState(0);

      gethostname(host_name, 255);

      FILE *f = fopen("env_details.asc", "w");

      fprintf(f, "%s\n", host_name);
      fclose(f);

      // fprintf (timings_ptr, "MPI 0 Hostname -> %s\n\n", host_name);

      //printf("kicking off network thread.....\n"); fflush(0x0);

      int sock_fd;
      int new_fd;
      int yes = 1;

      int is_broken_pipe = 0;
      int frame_number = 0;

      SteeringThread* steering_thread = NULL;

      static char ip_addr[16];
      static char rank_0_host_details[1024];

      signal(SIGPIPE, SIG_IGN); // Ignore a broken pipe

      HttpPost::get_host_details(rank_0_host_details, ip_addr);

      HttpPost::request("bunsen.chem.ucl.ac.uk", 28080, "/ahe/test/rendezvous/",
                        steering_session_id_char, rank_0_host_details);

      while (1)
      {
        setRenderState(0);

        // sem_wait( &nrl );

        struct sockaddr_in my_address;
        struct sockaddr_in their_addr; // client address

        socklen_t sin_size;

        if ( (sock_fd = socket(AF_INET, SOCK_STREAM, 0)) == -1)
        {
          perror("socket");
          exit(1);
        }
        if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1)
        {
          perror("setsockopt");
          exit(1);
        }
        my_address.sin_family = AF_INET;
        my_address.sin_port = htons(MYPORT);
        my_address.sin_addr.s_addr = INADDR_ANY;
        memset(my_address.sin_zero, '\0', sizeof my_address.sin_zero);

        if (bind(sock_fd, (struct sockaddr *) &my_address, sizeof my_address) == -1)
        {
          perror("bind");
          exit(1);
        }
        if (listen(sock_fd, CONNECTION_BACKLOG) == -1)
        {
          perror("listen");
          exit(1);
        }
        sin_size = sizeof (their_addr);

        // socklen_t _length = sizeof(my_address);
        //
        // getsockname (sock_fd, (struct sockaddr *) &my_address,&_length);
        // printf("Server Port is: %d\n", ntohs(my_address.sin_port));

        if ( (new_fd = accept(sock_fd, (struct sockaddr *) &their_addr, &sin_size)) == -1)
        {
          perror("accept");
          // continue;
        }

        // fprintf (timings_ptr, "server: got connection from %s (FD %i)\n", inet_ntoa (their_addr.sin_addr), new_fd);
        // printf ("RG thread: server: got connection from %s (FD %i)\n", inet_ntoa (their_addr.sin_addr), new_fd);
        steering_thread = new SteeringThread(new_fd, mSteeringController);
        steering_thread->Run();

        close(sock_fd);

        is_broken_pipe = 0;

        mSteeringController->isConnected.SetValue(true);

        // setRenderState(1);
        // At this point we're ready to send a frame...
        // setRendering=1;
        // sem_wait( &nrl );

        while (!is_broken_pipe)
        {
          // printf("THREAD: waiting for signal that frame is ready to send..\n"); fflush(0x0);

          bool is_frame_ready_local = 0;

          while (!is_frame_ready_local)
          {
            usleep(5000);
            sem_wait(&mSteeringController->nrl);
            is_frame_ready_local = mSteeringController->is_frame_ready;
            sem_post(&mSteeringController->nrl);
            // printf("THREAD is_frame_ready_local %i\n", is_frame_ready_local);
          }
          sem_wait(&mSteeringController->nrl);
          mSteeringController->sending_frame = 1;
          // printf("THREAD sending frame = 1\n");
          // setRenderState(0);
          // printf("THREAD: received signal that frame is ready to send..\n"); fflush(0x0);

          double frameTimeStart = frameTiming();

          int bytesSent = 0;

          {
            const int pixeldatabytes = 8;
            char *xdr_pixel = new char[pixeldatabytes];
            io::XdrMemWriter pixelWriter = io::XdrMemWriter(xdr_pixel, pixeldatabytes);

            pixelWriter << vis::controller->mScreen.PixelsX << vis::controller->mScreen.PixelsY;

            int pixelDataBytesSent = Network::send_all(new_fd, xdr_pixel, pixeldatabytes);

            if (pixelDataBytesSent < 0)
            {
              HandleBrokenPipe();
              is_broken_pipe = 1;
              break;
            }

            bytesSent += pixelDataBytesSent;

            delete xdr_pixel;
          }

          io::XdrMemWriter pixelDataWriter = io::XdrMemWriter(xdrSendBuffer_pixel_data,
                                                              pixel_data_bytes);
          io::XdrMemWriter frameDetailsWriter = io::XdrMemWriter(xdrSendBuffer_frame_details,
                                                                 frame_details_bytes);

          for (int i = 0; i < vis::controller->col_pixels_recv[RECV_BUFFER_A]; i++)
          {
            pixelDataWriter.writePixel(&vis::controller->col_pixel_recv[RECV_BUFFER_A][i],
                                       vis::ColourPalette::pickColour, mLbmParams->StressType);
          }

          int frameBytes = pixelDataWriter.getCurrentStreamPosition();

          frameDetailsWriter << frameBytes;

          int detailsBytes = frameDetailsWriter.getCurrentStreamPosition();

          int detailsBytesSent = Network::send_all(new_fd, xdrSendBuffer_frame_details,
                                                   detailsBytes);

          if (detailsBytesSent < 0)
          {
            HandleBrokenPipe();
            is_broken_pipe = 1;
            break;
          }
          else
          {
            bytesSent += detailsBytesSent;
          }

          int frameBytesSent = Network::send_all(new_fd, xdrSendBuffer_pixel_data, frameBytes);

          if (frameBytesSent < 0)
          {
            HandleBrokenPipe();
            is_broken_pipe = 1;
            break;
          }
          else
          {
            bytesSent += frameBytesSent;
          }

          SimulationParameters* sim = new SimulationParameters();
          sim->collectGlobalVals(mLbm, mSimState);
          int sizeToSend = sim->paramsSizeB;
          int simParamsBytesSent = Network::send_all(new_fd, sim->pack(), sizeToSend);

          if (simParamsBytesSent < 0)
          {
            HandleBrokenPipe();
            is_broken_pipe = 1;
            break;
          }

          bytesSent += simParamsBytesSent;

          // printf ("Sim bytes sent %i\n", sizeToSend);
          delete sim;

          // fprintf (timings_ptr, "bytes sent %i\n", bytesSent);
          // printf ("RG thread: bytes sent %i\n", bytesSent);

          setRenderState(1);

          double frameTimeSend = frameTiming() - frameTimeStart;

          // printf("Time to send frame = %0.6f s\n", frameTimeSend);

          double timeDiff = (1.0 / 25.0) - frameTimeSend;

          if (timeDiff > 0.0)
          {
            // printf("Sleeping for %0.6f s\n", timeDiff);

            usleep(timeDiff * 1.0e6);
          }

          // sem_post(&nrl);

          mSteeringController->sending_frame = 0;
          mSteeringController->is_frame_ready = 0;
          sem_post(&mSteeringController->nrl);

          frame_number++;

        } // while (is_broken_pipe == 0)

        close(new_fd);

        mSteeringController->isConnected.SetValue(false);

        // pthread_join(steering_thread, NULL);

      } // while(1)
    }

    void NetworkThread::HandleBrokenPipe()
    {
      printf("RG thread: broken network pipe...\n");
      mSteeringController->sending_frame = false;
      sem_post(&mSteeringController->nrl);
      setRenderState(0);
    }
  }

}
