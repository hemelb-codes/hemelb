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
                                           Network* iNetwork) :
      mNetwork(iNetwork), mLbm(lbm), mSimState(iSimState), mVisControl(iControl),
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

      isConnected = mNetwork->IsConnected();

      if (!isConnected)
      {
        return;
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
      bool success = mNetwork->send_all(xdrSendBuffer, imageWriter.getCurrentStreamPosition()
          - initialPosition);

      if (!success)
      {
        isConnected = false;
      }

      isFrameReady = false;
    }

    bool ImageSendComponent::ShouldRenderNewNetworkImage()
    {
      isConnected = mNetwork->IsConnected();

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
