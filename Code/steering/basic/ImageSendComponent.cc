#include <errno.h>
#include <signal.h>

#include "log/Logger.h"
#include "steering/ImageSendComponent.h"
#include "steering/basic/Network.h"
#include "io/XdrMemWriter.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace steering
  {
    // Use initialisation list to do the work.
    ImageSendComponent::ImageSendComponent(lb::LBM* lbm,
                                           lb::SimulationState* iSimState,
                                           vis::Control* iControl,
                                           const lb::LbmParameters* iLbmParams,
                                           ClientConnection* iClientConnection) :
      mClientConnection(iClientConnection), mLbm(lbm), mSimState(iSimState), mVisControl(iControl),
          mLbmParams(iLbmParams)
    {
      xdrSendBuffer = new char[maxSendSize];

      // Suppress signals from a broken pipe.
      signal(SIGPIPE, SIG_IGN);
      isFrameReady = false;
      isConnected = false;
      lastRender = 0.0;
    }

    ImageSendComponent::~ImageSendComponent()
    {
      delete[] xdrSendBuffer;
    }

    // TODO Need to check whether the current method (when the image is there, send the whole
    // thing in one go) is OK performance-wise, or whether we should instead try on another iteration
    // if we can't send it immediately.

    // This is original code with minimal tweaks to make it work with
    // the new (Feb 2011) structure.
    void ImageSendComponent::DoWork(void)
    {
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
      else
      {
        // Tell the steering controller that we have a connection.
        isConnected = true;
      }

      io::XdrMemWriter imageWriter = io::XdrMemWriter(xdrSendBuffer, maxSendSize);

      unsigned int initialPosition = imageWriter.getCurrentStreamPosition();

      // Write the dimensions of the image, in terms of pixel count.
      imageWriter << mVisControl->mScreen.GetPixelsX() << mVisControl->mScreen.GetPixelsY();

      // Write the length of the pixel data
      imageWriter << (int) (mVisControl->mScreen.GetPixelCount() * bytes_per_pixel_data);

      // Write the pixels themselves
      mVisControl->mScreen.WritePixels(&mVisControl->mDomainStats,
                                       &mVisControl->mVisSettings,
                                       &imageWriter);

      // Write the numerical data from the simulation, wanted by the client.
      {
        SimulationParameters sim;

        sim.timeStep = (int) mSimState->TimeStep;
        sim.time = mSimState->IntraCycleTime;
        sim.cycle = (int) mSimState->CycleId;
        sim.nInlets = mLbm->inlets;

        sim.mousePressure = mVisControl->mVisSettings.mouse_pressure;
        sim.mouseStress = mVisControl->mVisSettings.mouse_stress;

        mVisControl->mVisSettings.mouse_pressure = -1.0;
        mVisControl->mVisSettings.mouse_stress = -1.0;

        sim.pack(&imageWriter);
      }

      // Send to the client.
      if (SendSuccess(socketToClient, xdrSendBuffer, imageWriter.getCurrentStreamPosition()
          - initialPosition) < 0)
      {
        return;
      }

      isFrameReady = false;
    }

    ssize_t ImageSendComponent::SendSuccess(int iSocket, char * data, int length)
    {
      // Try to send all the data.
      ssize_t pixelDataBytesSent = Network::send_all(iSocket, data, length);

      // We couldn't send. The pipe is broken.
      if (pixelDataBytesSent <= 0)
      {
        if (errno != EAGAIN)
        {
          log::Logger::Log<log::Warning, log::Singleton>("Image send component: broken network pipe... (%s)",
                                                         strerror(errno));
          mClientConnection->ReportBroken(iSocket);
          isConnected = false;
        }
        return -1;
      }
      else
      {
        return pixelDataBytesSent;
      }
    }

    bool ImageSendComponent::ShouldRenderNewNetworkImage()
    {
      isConnected = mClientConnection->GetWorkingSocket() > 0;

      if (!isConnected)
      {
        return false;
      }

      // If we're going to exceed 25Hz by rendering now, wait until next iteration.
      {
        double frameTimeStart = util::myClock();

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
