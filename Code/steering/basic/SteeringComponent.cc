#include "steering/SteeringComponent.h"
#include "Network.h"
#include "io/XdrMemReader.h"
#include <errno.h>

namespace hemelb
{
  namespace steering
  {
    SteeringComponent::SteeringComponent(bool* iIsNetworkSending,
                                         int imagesPeriod,
                                         sem_t* varLock,
                                         ClientConnection* iClientConnection,
                                         vis::Control* iVisControl,
                                         lb::LBM* iLbm,
                                         net::Net * iNet,
                                         const topology::NetworkTopology *iNetTop,
                                         lb::SimulationState * iSimState) :
      net::PhasedBroadcast(iNet, iNetTop, iSimState, SPREADFACTOR), imagesPeriod(imagesPeriod),
          isNetworkSending(iIsNetworkSending), mVariableLock(varLock),
          mClientConnection(iClientConnection), mLbm(iLbm), mSimState(iSimState),
          mVisControl(iVisControl)
    {
      Reset();
    }

    /*
     * This function initialises all of the steering parameters, on all nodes.
     */
    void SteeringComponent::Reset()
    {
      isConnected = false;
      updatedMouseCoords = false;

      // scene center (dx,dy,dz)
      privateSteeringParams[0] = 0.0;
      privateSteeringParams[1] = 0.0;
      privateSteeringParams[2] = 0.0;

      // longitude and latitude
      privateSteeringParams[3] = 45.0;
      privateSteeringParams[4] = 45.0;

      // zoom and brightness
      privateSteeringParams[5] = 1.0;
      privateSteeringParams[6] = 0.03;

      // velocity and stress ranges
      privateSteeringParams[7] = 0.1;
      privateSteeringParams[8] = 0.1;

      // Minimum pressure and maximum pressure for Colour mapping
      privateSteeringParams[9] = 80.0;
      privateSteeringParams[10] = 120.0;

      // Glyph length
      privateSteeringParams[11] = 1.0;

      // Rendered frame size, pixel x and pixel y
      privateSteeringParams[12] = 512;
      privateSteeringParams[13] = 512;

      // x-y position of the mouse of the client
      privateSteeringParams[14] = -1.0;
      privateSteeringParams[15] = -1.0;

      // signal useful to terminate the simulation
      privateSteeringParams[16] = 0.0;

      // Vis_mode
      privateSteeringParams[17] = 0.0;

      // vis_streaklines_per_pulsatile_period
      privateSteeringParams[18] = 5.0;

      // vis_streakline_length
      privateSteeringParams[19] = 100.0;

      // Value of DoRendering
      privateSteeringParams[20] = 0.0;
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
        privateSteeringParams[STEERABLE_PARAMETERS] = (float) ( (isConnected && !*isNetworkSending)
            || (mSimState->TimeStep % imagesPeriod >= 0 && mSimState->TimeStep % imagesPeriod
                < (int) GetRoundTripLength()));

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
          int steerDataRecvB = Network::recv_all(socket, steeringRecvBuffer, num_chars);

          // If there was an error, report it and return.
          if (steerDataRecvB < 0)
          {
            // If there was no data and it wasn't simply that the socket would block,
            // raise an error.
            if (steerDataRecvB != EAGAIN)
            {
              printf("Steering component: broken network pipe...\n");
              mClientConnection->ReportBroken(socket);
              isConnected = false;
            }
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
      mVisControl->ctr_x += privateSteeringParams[0];
      mVisControl->ctr_y += privateSteeringParams[1];
      mVisControl->ctr_z += privateSteeringParams[2];

      float longitude = privateSteeringParams[3];
      float latitude = privateSteeringParams[4];

      float zoom = privateSteeringParams[5];

      mVisControl->brightness = privateSteeringParams[6];

      // The minimum value here is by default 0.0 all the time
      mVisControl->physical_velocity_threshold_max = privateSteeringParams[7];

      // The minimum value here is by default 0.0 all the time
      mVisControl->physical_stress_threshold_max = privateSteeringParams[8];

      mVisControl->physical_pressure_threshold_min = privateSteeringParams[9];
      mVisControl->physical_pressure_threshold_max = privateSteeringParams[10];

      vis::GlyphDrawer::glyph_length = privateSteeringParams[11];

      float pixels_x = privateSteeringParams[12];
      float pixels_y = privateSteeringParams[13];

      int newMouseX = int (privateSteeringParams[14]);
      int newMouseY = int (privateSteeringParams[15]);

      if (newMouseX != mVisControl->mouse_x || newMouseY != mVisControl->mouse_y)
      {
        updatedMouseCoords = true;
        mVisControl->mouse_x = newMouseX;
        mVisControl->mouse_y = newMouseY;
      }

      mSimState->IsTerminating = int (privateSteeringParams[16]);

      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      mVisControl->mode = int (privateSteeringParams[17]);

      mVisControl->streaklines_per_pulsatile_period = privateSteeringParams[18];
      mVisControl->streakline_length = privateSteeringParams[19];

      mSimState->DoRendering = int (privateSteeringParams[20]);

      mVisControl->updateImageSize(pixels_x, pixels_y);

      float lattice_density_min =
          mLbm->ConvertPressureToLatticeUnits(mVisControl->physical_pressure_threshold_min) / Cs2;
      float lattice_density_max =
          mLbm->ConvertPressureToLatticeUnits(mVisControl->physical_pressure_threshold_max) / Cs2;
      float lattice_velocity_max =
          mLbm->ConvertVelocityToLatticeUnits(mVisControl->physical_velocity_threshold_max);
      float lattice_stress_max =
          mLbm->ConvertStressToLatticeUnits(mVisControl->physical_stress_threshold_max);

      mVisControl->SetProjection(pixels_x, pixels_y, mVisControl->ctr_x, mVisControl->ctr_y,
                                 mVisControl->ctr_z, longitude, latitude, zoom);

      mVisControl->density_threshold_min = lattice_density_min;
      mVisControl->density_threshold_minmax_inv = 1.0F
          / (lattice_density_max - lattice_density_min);
      mVisControl->velocity_threshold_max_inv = 1.0F / lattice_velocity_max;
      mVisControl->stress_threshold_max_inv = 1.0F / lattice_stress_max;
    }

  }
}
