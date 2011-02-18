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
#include "steering/Control.h"
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
                                 const lb::LbmParameters* iLbmParams,
                                 ClientConnection* iClientConnection) :
      mClientConnection(iClientConnection), mLbm(lbm), mSteeringController(steeringController),
          mSimState(iSimState), mLbmParams(iLbmParams)
    {
      xdrSendBuffer_pixel_data = new char[pixel_data_bytes];
    }

    NetworkThread::~NetworkThread()
    {
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
      if(val > 0)
      {
        sem_wait(&mSimState->Rendering);
      }
      else
      {
        sem_post(&mSimState->Rendering);
      }
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
      // Suppress signals from a broken pipe.
      signal(SIGPIPE, SIG_IGN);

      // Loop forever (no breaks), creating new sockets whenever the current pipe breaks.
      while (true)
      {
        int socketToClient = mClientConnection->GetWorkingSocket();

        // Start the steering thread on the connection to the client.
        {
          SteeringThread* steering_thread =
              new SteeringThread(socketToClient, &mSteeringController->steering_var_lock);
          steering_thread->Run();
        }

        // Tell the steering controller that we have a connection.
        mSteeringController->isConnected.SetValue(true);

        bool is_broken_pipe = false;

        // While the connection is still live, keep sending frames.
        while (!is_broken_pipe)
        {
          // Wait until a frame is ready.
          {
            bool is_frame_ready_local = false;

            while (!is_frame_ready_local)
            {
              usleep(5000);
              sem_wait(&mSteeringController->nrl);
              is_frame_ready_local = mSteeringController->is_frame_ready;
              sem_post(&mSteeringController->nrl);
            }
          }

          // Take control of the steering controller.
          sem_wait(&mSteeringController->nrl);
          mSteeringController->sending_frame = 1;
          // Turn off rendering.
          setRenderState(0);

          double frameTimeStart = frameTiming();

          int bytesSent = 0;

          // Send the dimensions of the image, in terms of pixel count.
          {
            // 2 ints for the X and Y dimensions = 2*4 bytes
            const int pixeldatabytes = 8;
            char xdr_pixel[pixeldatabytes];
            io::XdrMemWriter pixelWriter = io::XdrMemWriter(xdr_pixel, pixeldatabytes);

            pixelWriter << vis::controller->mScreen.PixelsX << vis::controller->mScreen.PixelsY;

            int pixelDataBytesSent = Network::send_all(socketToClient, xdr_pixel, pixeldatabytes);

            if (pixelDataBytesSent < 0)
            {
              HandleBrokenPipe();
              is_broken_pipe = true;
              break;
            }

            bytesSent += pixelDataBytesSent;
          }

          io::XdrMemWriter pixelDataWriter(xdrSendBuffer_pixel_data, pixel_data_bytes);

          for (int i = 0; i < vis::controller->col_pixels_recv[RECV_BUFFER_A]; i++)
          {
            pixelDataWriter.writePixel(&vis::controller->col_pixel_recv[RECV_BUFFER_A][i],
                                       vis::ColourPalette::pickColour, mLbmParams->StressType);
          }

          // Send the number of bytes being used on pixel data.
          int frameBytes = pixelDataWriter.getCurrentStreamPosition();
          {
            char xdrSendBuffer_frame_details[frame_details_bytes];

            io::XdrMemWriter frameDetailsWriter = io::XdrMemWriter(xdrSendBuffer_frame_details,
                                                                   frame_details_bytes);
            frameDetailsWriter << frameBytes;

            int frameDetailsBytes = frameDetailsWriter.getCurrentStreamPosition();

            int frameDetailsBytesSent = Network::send_all(socketToClient,
                                                          xdrSendBuffer_frame_details,
                                                          frameDetailsBytes);

            if (frameDetailsBytesSent < 0)
            {
              HandleBrokenPipe();
              is_broken_pipe = true;
              break;
            }
            else
            {
              bytesSent += frameDetailsBytesSent;
            }
          }

          int frameBytesSent = Network::send_all(socketToClient, xdrSendBuffer_pixel_data,
                                                 frameBytes);

          if (frameBytesSent < 0)
          {
            HandleBrokenPipe();
            is_broken_pipe = true;
            break;
          }
          else
          {
            bytesSent += frameBytesSent;
          }

          // Send the numerical data from the simulation, wanted by the client.
          {
            SimulationParameters sim;
            sim.collectGlobalVals(mLbm, mSimState);
            int sizeToSend = sim.paramsSizeB;
            int simParamsBytesSent = Network::send_all(socketToClient, sim.pack(), sizeToSend);

            if (simParamsBytesSent < 0)
            {
              HandleBrokenPipe();
              is_broken_pipe = true;
              break;
            }

            bytesSent += simParamsBytesSent;
          }

          // All data is sent so we can start rendering again.
          setRenderState(1);

          // Send this thread to sleep to aim for an approximate 25Hz of sends.
          {
            double sendingTime = frameTiming() - frameTimeStart;

            if ( (1.0 / 25.0) > sendingTime)
            {
              usleep( ( (1.0 / 25.0) - sendingTime) * 1.0e6);
            }
          }

          mSteeringController->sending_frame = 0;
          mSteeringController->is_frame_ready = 0;
          sem_post(&mSteeringController->nrl);
        } // while (is_broken_pipe == 0)

        mSteeringController->isConnected.SetValue(false);
      } // while(1)
    }

    void NetworkThread::HandleBrokenPipe()
    {
      mClientConnection->ReportBroken();
      printf("RG thread: broken network pipe...\n");
      mSteeringController->sending_frame = false;
      sem_post(&mSteeringController->nrl);
      setRenderState(0);
    }
  }

}
