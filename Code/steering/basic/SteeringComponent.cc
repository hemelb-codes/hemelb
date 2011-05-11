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
                                         ClientConnection* iClientConnection,
                                         vis::Control* iVisControl,
                                         lb::LBM* iLbm,
                                         net::Net * iNet,
                                         lb::SimulationState * iSimState) :
      net::PhasedBroadcast<false, 1, 0, true, false>(iNet, iSimState, SPREADFACTOR),
          imagesPeriod(imagesPeriod), mClientConnection(iClientConnection), mLbm(iLbm),
          mSimState(iSimState), mVisControl(iVisControl)
    {
      Reset();
    }

    /**
     * Does nothing - data only travels from parent to children in steering.
     */
    void SteeringComponent::ProgressFromChildren()
    {

    }

    void SteeringComponent::ProgressFromParent()
    {
      ReceiveFromParent<float> (privateSteeringParams, STEERABLE_PARAMETERS + 1);
    }

    void SteeringComponent::ProgressToChildren()
    {
      SendToChildren<float> (privateSteeringParams, STEERABLE_PARAMETERS + 1);
    }

    /**
     * Does nothing - data only travels from parent to children in steering.
     */
    void SteeringComponent::ProgressToParent()
    {

    }

    void SteeringComponent::PostReceiveFromChildren()
    {

    }

    void SteeringComponent::PostReceiveFromParent()
    {

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
      int socket = mClientConnection->GetWorkingSocket();

      if (socket < 0)
      {
        isConnected = false;
        return;
      }
      else
      {
        isConnected = true;
      }

      bool newSteeringDataExists = false;

      // Keep looping through all data that we can read - if multiple steering commands have been
      // received we only want to act on the last set.
      while (true)
      {
        // Create a set of sockets to check for readability, just containing the single socket
        // of interest.
        fd_set readableDescriptors;
        {
          struct timeval tv;
          tv.tv_sec = 0;
          tv.tv_usec = 0;

          FD_ZERO(&readableDescriptors);
          FD_SET(socket, &readableDescriptors);

          select(socket + 1, &readableDescriptors, NULL, NULL, &tv);
        }

        // If the socket is readable, read from it.
        if (FD_ISSET(socket, &readableDescriptors))
        {
          // Try to receive all the data.
          ssize_t steerDataRecvB = Network::recv_all(socket, steeringRecvBuffer, num_chars);

          // If there was an error, report it and return.
          if (steerDataRecvB < 0)
          {
            // If there was no data and it wasn't simply that the socket would block,
            // raise an error.
            if (errno != EAGAIN)
            {
              log::Logger::Log<log::Warning, log::Singleton>("Steering component: broken network pipe... (%s)",
                                                             strerror(errno));
              mClientConnection->ReportBroken(socket);
              isConnected = false;
            }
            break;
          }
          else if (steerDataRecvB == 0)
          {
            break;
          }

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
