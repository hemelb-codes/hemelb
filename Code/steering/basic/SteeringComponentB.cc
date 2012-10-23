// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <cerrno>

#include "steering/SteeringComponent.h"
#include "steering/Network.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace steering
  {
    SteeringComponent::SteeringComponent(Network* iNetwork,
                                         vis::Control* iVisControl,
                                         steering::ImageSendComponent* imageSendComponent,
                                         net::Net * iNet,
                                         lb::SimulationState * iSimState,
                                         configuration::SimConfig* iSimConfig,
                                         util::UnitConverter* iUnits) :
        net::PhasedBroadcastRegular<false, 1, 0, true, false>(iNet, iSimState, SPREADFACTOR),
        mNetwork(iNetwork), mSimState(iSimState), mVisControl(iVisControl), imageSendComponent(imageSendComponent),
        mUnits(iUnits),simConfig(iSimConfig)
    {
      ClearValues();
      AssignValues();
    }

    void SteeringComponent::ProgressFromParent(unsigned long splayNumber)
    {
      ReceiveFromParent<float>(privateSteeringParams, STEERABLE_PARAMETERS + 1);
    }

    void SteeringComponent::ProgressToChildren(unsigned long splayNumber)
    {
      SendToChildren<float>(privateSteeringParams, STEERABLE_PARAMETERS + 1);
    }

    bool SteeringComponent::RequiresSeparateSteeringCore()
    {
      return true;
    }

    void SteeringComponent::TopNodeAction()
    {
      /*
       * The final steering parameter is DoRendering, which is true if we're connected and ready
       * for the next frame.
       */
      {
        privateSteeringParams[STEERABLE_PARAMETERS] = (float) (isConnected && readyForNextImage);
      }

      // Create a buffer for the data received.
      const int num_chars = STEERABLE_PARAMETERS * sizeof(float) / sizeof(char);
      const int bytes = sizeof(char) * num_chars;

      char steeringRecvBuffer[bytes];

      // Get the open socket.
      isConnected = mNetwork->IsConnected();

      if (!isConnected)
      {
        return;
      }

      bool newSteeringDataExists = false;

      // Keep looping through all data that we can read - if multiple steering commands have been
      // received we only want to act on the last set.
      while (true)
      {
        // Try to receive all the data.
        if (mNetwork->recv_all(steeringRecvBuffer, num_chars))
        {
          // We read data!
          newSteeringDataExists = true;
        }
        // If it wasn't readable, we've exhausted all the received data at the socket, so
        // break and actually use the last read data.
        else
        {
          break;
        }
      }

      if (newSteeringDataExists)
      {
        // Initialise the stream reader.
        io::writers::xdr::XdrMemReader steeringStream(steeringRecvBuffer, bytes);

        for (int i = 0; i < STEERABLE_PARAMETERS; i++)
          steeringStream.readFloat(privateSteeringParams[i]);
      }
    }

    void SteeringComponent::Effect()
    {
      // TODO we need to make sure that doing this doesn't overwrite the values in the config.xml file.
      // At the moment, it definitely does.
      AssignValues();
    }
  }
}
