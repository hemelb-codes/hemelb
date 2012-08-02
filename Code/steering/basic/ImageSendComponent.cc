#include <cerrno>
#include <csignal>

#include "log/Logger.h"
#include "steering/ImageSendComponent.h"
#include "steering/Network.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace steering
  {
    // Use initialisation list to do the work.
    ImageSendComponent::ImageSendComponent(lb::SimulationState* iSimState,
                                           vis::Control* iControl,
                                           const lb::LbmParameters* iLbmParams,
                                           Network* iNetwork,
                                           unsigned inletCountIn) :
        mNetwork(iNetwork), mSimState(iSimState), mVisControl(iControl), mLbmParams(iLbmParams), inletCount(inletCountIn), MaxFramerate(25.0)
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
    void ImageSendComponent::DoWork(const vis::PixelSet<vis::ResultPixel>* pix)
    {
      isConnected = mNetwork->IsConnected();

      if (!isConnected)
      {
        return;
      }

      io::writers::xdr::XdrMemWriter imageWriter = io::writers::xdr::XdrMemWriter(xdrSendBuffer, maxSendSize);

      unsigned int initialPosition = imageWriter.getCurrentStreamPosition();

      // Write the dimensions of the image, in terms of pixel count.
      imageWriter << mVisControl->GetPixelsX() << mVisControl->GetPixelsY();

      // Write the length of the pixel data
      imageWriter << (int) (pix->GetPixelCount() * bytes_per_pixel_data);

      // Write the pixels themselves
      mVisControl->WritePixels(&imageWriter, *pix, mVisControl->domainStats, mVisControl->visSettings);

      // Write the numerical data from the simulation, wanted by the client.
      {
        SimulationParameters sim;

        sim.timeStep = (int) mSimState->GetTimeStep();
        sim.time = mSimState->GetTime();
        sim.nInlets = inletCount;

        sim.mousePressure = mVisControl->visSettings.mouse_pressure;
        sim.mouseStress = mVisControl->visSettings.mouse_stress;

        mVisControl->visSettings.mouse_pressure = -1.0;
        mVisControl->visSettings.mouse_stress = -1.0;

        sim.pack(&imageWriter);
      }

      // Send to the client.
      log::Logger::Log<log::Debug, log::Singleton>("Sending network image at timestep %d",mSimState->GetTimeStep());
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

        if ( (1.0 / MaxFramerate) > deltaTime)
        {
          return false;
        }
        else
        {
          log::Logger::Log<log::Trace, log::Singleton>("Image-send component requesting new render, %f seconds since last one at step %d max rate is %f.",
                                                        deltaTime, mSimState->GetTimeStep(), MaxFramerate);
          lastRender = frameTimeStart;
          return true;
        }
      }
    }

  }
}
