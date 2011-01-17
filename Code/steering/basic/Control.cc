/*
 * Control.cc
 *
 *  Created on: Oct 27, 2010
 *      Author: rupert
 */
#include <pthread.h>
#include <semaphore.h>

#include "lb/lb.h"
#include "steering/basic/Control.h"
#include "steering/common/common.h"

#include "vis/visthread.h"

namespace hemelb
{

  namespace steering
  {

    Control::Control(bool isCurrentProcTheSteeringProc) :
      is_frame_ready(false), sending_frame(false), isConnected(false),
          updated_mouse_coords(false),
          mIsCurrentProcTheSteeringProc(isCurrentProcTheSteeringProc)
    {
      // From main()
      sem_init(&nrl, 0, 1);
      //sem_init(&connected_sem, 0, 1);
      sem_init(&steering_var_lock, 0, 1);
    }

    Control::~Control()
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        delete[] hemelb::vis::xdrSendBuffer_frame_details;
        delete[] hemelb::vis::xdrSendBuffer_pixel_data;
        delete mNetworkThread;
      }
      sem_destroy(&nrl);
      //sem_destroy(&connected_sem);
      sem_destroy(&steering_var_lock);
    }

    // Kick off the networking thread
    void Control::StartNetworkThread(lb::LBM* lbm,
                                     lb::SimulationState *iSimState,
                                     const lb::LbmParameters *iLbmParams)
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        hemelb::vis::xdrSendBuffer_pixel_data
            = new char[hemelb::vis::pixel_data_bytes];
        hemelb::vis::xdrSendBuffer_frame_details
            = new char[hemelb::vis::frame_details_bytes];

        pthread_mutex_init(&LOCK, NULL);
        pthread_cond_init(&network_send_frame, NULL);

        mNetworkThread = new NetworkThread(lbm, this, iSimState, iLbmParams);
        mNetworkThread->Run();
      }
    }

    // Broadcast the steerable parameters to all tasks.
    void Control::UpdateSteerableParameters(bool shouldRenderForSnapshot,
                                            int* perform_rendering,
                                            hemelb::lb::SimulationState &iSimulationState,
                                            hemelb::vis::Control* visController,
                                            lb::LBM* lbm)
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        hemelb::vis::doRendering = ShouldRenderForNetwork()
            || shouldRenderForSnapshot;
        sem_wait(&steering_var_lock);
      }

      BroadcastSteerableParameters(perform_rendering, iSimulationState,
                                   visController, lbm);

      if (mIsCurrentProcTheSteeringProc)
        sem_post(&steering_var_lock);
    }

    // Broadcast the steerable parameters to all tasks.
    void Control::BroadcastSteerableParameters(int *perform_rendering,
                                               hemelb::lb::SimulationState &lSimulationState,
                                               vis::Control *visControl,
                                               lb::LBM* lbm)
    {
      steer_par[STEERABLE_PARAMETERS] = *perform_rendering;

      MPI_Bcast(steer_par, STEERABLE_PARAMETERS + 1, MPI_FLOAT, 0,
                MPI_COMM_WORLD);

      float longitude, latitude;
      float zoom;
      float lattice_velocity_max, lattice_stress_max;
      float lattice_density_min, lattice_density_max;

      int pixels_x, pixels_y;

      visControl->ctr_x += steer_par[0];
      visControl->ctr_y += steer_par[1];
      visControl->ctr_z += steer_par[2];

      longitude = steer_par[3];
      latitude = steer_par[4];

      zoom = steer_par[5];

      visControl->brightness = steer_par[6];

      // The minimum value here is by default 0.0 all the time
      visControl->physical_velocity_threshold_max = steer_par[7];

      // The minimum value here is by default 0.0 all the time
      visControl->physical_stress_threshold_max = steer_par[8];

      visControl->physical_pressure_threshold_min = steer_par[9];
      visControl->physical_pressure_threshold_max = steer_par[10];

      vis::GlyphDrawer::glyph_length = steer_par[11];

      pixels_x = steer_par[12];
      pixels_y = steer_par[13];

      visControl->mouse_x = int (steer_par[14]);
      visControl->mouse_y = int (steer_par[15]);

      lSimulationState.IsTerminating = int (steer_par[16]);

      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      visControl->mode = int (steer_par[17]);

      visControl->streaklines_per_pulsatile_period = steer_par[18];
      visControl->streakline_length = steer_par[19];

      *perform_rendering = int (steer_par[20]);

      visControl->updateImageSize(pixels_x, pixels_y);

      lattice_density_min
          = lbm->ConvertPressureToLatticeUnits(
                                               visControl->physical_pressure_threshold_min)
              / Cs2;
      lattice_density_max
          = lbm->ConvertPressureToLatticeUnits(
                                               visControl->physical_pressure_threshold_max)
              / Cs2;
      lattice_velocity_max
          = lbm->ConvertVelocityToLatticeUnits(
                                               visControl->physical_velocity_threshold_max);
      lattice_stress_max
          = lbm->ConvertStressToLatticeUnits(
                                             visControl->physical_stress_threshold_max);

      visControl->SetProjection(pixels_x, pixels_y, visControl->ctr_x,
                                visControl->ctr_y, visControl->ctr_z,
                                longitude, latitude, zoom);

      visControl->density_threshold_min = lattice_density_min;
      visControl->density_threshold_minmax_inv = 1.0F / (lattice_density_max
          - lattice_density_min);
      visControl->velocity_threshold_max_inv = 1.0F / lattice_velocity_max;
      visControl->stress_threshold_max_inv = 1.0F / lattice_stress_max;
    }

    // Do we need to render a frame for the client?
    bool Control::ShouldRenderForNetwork()
    {
      // Iff we're connected and not sending a
      // frame, we'd better render a new one!
      return isConnected.GetValue() && !sending_frame;
    }

  }

}
