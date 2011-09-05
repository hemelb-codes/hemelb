#include <errno.h>
#include <signal.h>

#include "log/Logger.h"
#include "steering/ImageSendComponent.h"
#include "steering/Network.h"
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

      isConnected = false;
      lastRender = 0.0;
    }

    ImageSendComponent::~ImageSendComponent()
    {
      delete[] xdrSendBuffer;
    }

    // This is original code with minimal tweaks to make it work with
    // the new (Feb 2011) structure.
    void ImageSendComponent::DoWork(const vis::ScreenPixels<vis::RayDataType_t>* pix)
    {
      isConnected = mNetwork->IsConnected();

      if (!isConnected)
      {
        return;
      }

      io::XdrMemWriter imageWriter = io::XdrMemWriter(xdrSendBuffer, maxSendSize);

      unsigned int initialPosition = imageWriter.getCurrentStreamPosition();

      // Write the dimensions of the image, in terms of pixel count.
      imageWriter << pix->GetPixelsX() << pix->GetPixelsY();

      // Write the length of the pixel data
      imageWriter << (int) (pix->GetStoredPixelCount() * bytes_per_pixel_data);

      // Write the pixels themselves
      pix->WritePixels(&imageWriter, &mVisControl->mDomainStats, &mVisControl->mVisSettings);

      // Write the numerical data from the simulation, wanted by the client.
      {
        SimulationParameters sim;

        sim.timeStep = (int) mSimState->GetTimeStep();
        sim.time = mSimState->GetIntraCycleTime();
        sim.cycle = (int) mSimState->GetCycleId();
        sim.nInlets = mLbm->inlets;

        sim.mousePressure = mVisControl->mVisSettings.mouse_pressure;
        sim.mouseStress = mVisControl->mVisSettings.mouse_stress;

        mVisControl->mVisSettings.mouse_pressure = -1.0;
        mVisControl->mVisSettings.mouse_stress = -1.0;

        sim.pack(&imageWriter);
      }

      // Send to the client.
      mNetwork->send_all(xdrSendBuffer, imageWriter.getCurrentStreamPosition() - initialPosition);
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
          log::Logger::Log<log::Debug, log::OnePerCore>("Image-send component requesting new render, %f seconds since last one.",
                                                        deltaTime);
          lastRender = frameTimeStart;
          return true;
        }
      }
    }

  }
}
