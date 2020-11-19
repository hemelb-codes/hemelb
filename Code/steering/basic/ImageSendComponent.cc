// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
    ImageSendComponent::ImageSendComponent(lb::SimulationState* iSimState, vis::Control* iControl,
                                           const lb::LbmParameters* iLbmParams, Network* iNetwork,
                                           unsigned inletCountIn) :
      isConnected(false),
      mNetwork(iNetwork), mSimState(iSimState), mVisControl(iControl),
      inletCount(inletCountIn), MaxFramerate(25.0),
      xdrSendBuffer(new char[maxSendSize]),
      lastRender(0.0)
    {
      // Suppress signals from a broken pipe.
      signal(SIGPIPE, SIG_IGN);
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

      auto imageWriter = io::writers::xdr::XdrMemWriter(xdrSendBuffer.get(), maxSendSize);

      unsigned int initialPosition = imageWriter.getCurrentStreamPosition();

      // Write the dimensions of the image, in terms of pixel count.
      imageWriter << mVisControl->GetPixelsX() << mVisControl->GetPixelsY();

      // Write the length of the pixel data
      imageWriter << (int) (pix->GetPixelCount() * bytes_per_pixel_data);

      // Write the pixels themselves
      mVisControl->WritePixels(&imageWriter,
                               *pix,
                               mVisControl->domainStats,
                               mVisControl->visSettings);

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
      mNetwork->send_all(xdrSendBuffer.get(), imageWriter.getCurrentStreamPosition() - initialPosition);
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
                                                       deltaTime,
                                                       mSimState->GetTimeStep(),
                                                       MaxFramerate);
          lastRender = frameTimeStart;
          return true;
        }
      }
    }

  }
}
