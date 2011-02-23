#include "steering/ImageSendComponent.h"
#include "steering/basic/SimulationParameters.h"
#include "steering/basic/Network.h"
#include "io/XdrMemWriter.h"
#include <signal.h>

namespace hemelb
{
  namespace steering
  {
    // Use initialisation list to do the work.
    ImageSendComponent::ImageSendComponent(lb::LBM* lbm,
                                           lb::SimulationState* iSimState,
                                           const lb::LbmParameters* iLbmParams,
                                           ClientConnection* iClientConnection) :
      mClientConnection(iClientConnection), mLbm(lbm), mSimState(iSimState), mLbmParams(iLbmParams)
    {
      xdrSendBuffer_pixel_data = new char[pixel_data_bytes];

      // Suppress signals from a broken pipe.
      signal(SIGPIPE, SIG_IGN);
      isFrameReady = false;
      lastRender = 0.0;
      double lastSend = 0.0;
    }

    ImageSendComponent::~ImageSendComponent()
    {
      delete[] xdrSendBuffer_pixel_data;
    }

    // TODO The timing of network sends (i.e. the frame rate) needs to be checked.
    // TODO Need to check whether the current method (when the image is there, send the whole
    // thing in one go) is OK performance-wise, or whether we should instead try on another iteration
    // if we can't send it immediately.

    // This is original code with minimal tweaks to make it work with
    // the new (Feb 2011) structure.
    void ImageSendComponent::DoWork(void)
    {
      debug::Debugger::Get()->BreakHere();

      // If no frame is ready for sending, return.
      if (!isFrameReady)
      {
        return;
      }

      // Get a socket.
      int socketToClient = mClientConnection->GetWorkingSocket();

      // If it's non-existent, we don't have a connection. Nothing to do.
      if (socketToClient < 0)
      {
        isConnected = false;
        return;
      }

      // Tell the steering controller that we have a connection.
      isConnected = true;

      int bytesSent = 0;

      // Send the dimensions of the image, in terms of pixel count.
      {
        // 2 ints for the X and Y dimensions = 2*4 bytes
        const int pixeldatabytes = 8;
        char xdr_pixel[pixeldatabytes];
        io::XdrMemWriter pixelWriter = io::XdrMemWriter(xdr_pixel, pixeldatabytes);

        pixelWriter << vis::controller->mScreen.PixelsX << vis::controller->mScreen.PixelsY;

        int pixelDataBytesSent = SendSuccess(socketToClient, xdr_pixel, pixeldatabytes);

        if (pixelDataBytesSent < 0)
        {
          return;
        }
        else
        {
          bytesSent += pixelDataBytesSent;
        }
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

        int frameDetailsBytesSent = SendSuccess(socketToClient, xdrSendBuffer_frame_details,
                                                frameDetailsBytes);

        if (frameDetailsBytesSent < 0)
        {
          return;
        }
        else
        {
          bytesSent += frameDetailsBytesSent;
        }
      }

      int frameBytesSent = SendSuccess(socketToClient, xdrSendBuffer_pixel_data, frameBytes);

      if (frameBytesSent < 0)
      {
        return;
      }
      else
      {
        bytesSent += frameBytesSent;
      }

      // Send the numerical data from the simulation, wanted by the client.
      {
        SimulationParameters sim;

        sim.pressureMin = mLbm->GetMinPhysicalPressure();
        sim.pressureMax = mLbm->GetMaxPhysicalPressure();
        sim.velocityMin = mLbm->GetMinPhysicalVelocity();
        sim.velocityMax = mLbm->GetMaxPhysicalVelocity();
        sim.stressMax = mLbm->GetMaxPhysicalStress();
        sim.timeStep = mSimState->TimeStep;
        sim.time = mSimState->IntraCycleTime;
        sim.cycle = mSimState->CycleId;
        sim.nInlets = mLbm->inlets;

        sim.mousePressure = vis::controller->mouse_pressure;
        sim.mouseStress = vis::controller->mouse_stress;

        int sizeToSend = sim.paramsSizeB;
        int simParamsBytesSent = SendSuccess(socketToClient, sim.pack(), sizeToSend);

        if (simParamsBytesSent < 0)
        {
          return;
        }
        else
        {
          bytesSent += simParamsBytesSent;
        }
      }

      isFrameReady = false;
    }

    int ImageSendComponent::SendSuccess(int iSocket, char * data, int length)
    {
      // Try to send all the data.
      int pixelDataBytesSent = Network::send_all(iSocket, data, length);

      // We couldn't send. The pipe is broken.
      if (pixelDataBytesSent < 0)
      {
        mClientConnection->ReportBroken(iSocket);
        isConnected = false;
        isFrameReady = false;
        return -1;
      }
      else
      {
        return pixelDataBytesSent;
      }
    }

    bool ImageSendComponent::ShouldRenderNewNetworkImage()
    {
      if (!isConnected)
      {
        return false;
      }

      // If we're going to exceed 25Hz by rendering now, wait until next iteration.
      {
        double frameTimeStart = frameTiming();

        double deltaTime = frameTimeStart - lastRender;

        if ( (1.0 / 25.0) > deltaTime)
        {
          return false;
        }
        else
        {
          // TODO this isn't ideal. Change it.
          lastRender = frameTimeStart;
          return true;
        }
      }
    }

  }
}
