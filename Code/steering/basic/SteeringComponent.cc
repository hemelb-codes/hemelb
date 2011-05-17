#include <errno.h>

#include "steering/SteeringComponent.h"
#include "Network.h"
#include "io/XdrMemReader.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace steering
  {
    SteeringComponent::SteeringComponent(int imagesPeriod,
                                         Network* iNetwork,
                                         vis::Control* iVisControl,
                                         lb::LBM* iLbm,
                                         net::Net * iNet,
                                         lb::SimulationState * iSimState) :
      net::PhasedBroadcast<false, 1, 0, true, false>(iNet, iSimState, SPREADFACTOR),
          imagesPeriod(imagesPeriod), mNetwork(iNetwork), mLbm(iLbm), mSimState(iSimState),
          mVisControl(iVisControl)
    {
      Reset();
    }

    void SteeringComponent::ProgressFromParent(unsigned int splayNumber)
    {
      ReceiveFromParent<float> (privateSteeringParams, STEERABLE_PARAMETERS + 1);
    }

    void SteeringComponent::ProgressToChildren(unsigned int splayNumber)
    {
      SendToChildren<float> (privateSteeringParams, STEERABLE_PARAMETERS + 1);
    }

    bool SteeringComponent::RequiresSeparateSteeringCore()
    {
      return true;
    }

    void SteeringComponent::TopNodeAction()
    {
      /*
       * The final steering parameter is DoRendering, which is true if we're connected and not in
       * the middle of sending another image, OR if we're writing an image file. Get that out of
       * the way early on.
       */
      {
        privateSteeringParams[STEERABLE_PARAMETERS] = (float) ( (isConnected)
            || (mSimState->TimeStep % imagesPeriod < GetRoundTripLength()));

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
        io::XdrMemReader steeringStream(steeringRecvBuffer, bytes);

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
