#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "usage.h"
#include "fileutils.h"
#include "utilityFunctions.h"
#include "lb.h"

#include "SimulationMaster.h"

#include "steering/steering.h"

#include "vis/visthread.h"
#include "vis/ColourPalette.h"

#include "debug/Debugger.h"

int cycle_id;
int time_step;
double intra_cycle_time;

#ifndef NO_STEER

#endif

FILE *timings_ptr;

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

#ifndef NO_STEER
  sem_init(&hemelb::steering::nrl, 0, 1);
  sem_init(&hemelb::steering::connected_sem, 0, 1);
  sem_init(&hemelb::steering::steering_var_lock, 0, 1);

  hemelb::steering::is_frame_ready = 0;
  hemelb::steering::connected = 0;
  hemelb::steering::sending_frame = 0;
  hemelb::steering::updated_mouse_coords = 0;
#endif

  double simulation_time = 0.;
  int total_time_steps, stability = STABLE;
  int depths;

  // #ifdef NO_STEER
  //   int doRendering;
  // #endif

  int snapshots_period;
  int images_period;
  int is_unstable = 0;

#ifndef NO_STEER
  pthread_t network_thread;
  pthread_attr_t pthread_attrib;
#endif

  SimulationMaster lMaster = SimulationMaster(argc, argv);

  double total_time = hemelb::util::myClock();

  char* input_file_path(argv[1]);

  char input_config_name[256];
  char input_parameters_name[256];
  char vis_parameters_name[256];
  char snapshot_directory[256];
  char image_directory[256];
  char complete_image_name[256];

  strcpy(input_config_name, input_file_path);
  strcat(input_config_name, "/config.dat");
  hemelb::util::check_file(input_config_name);

  strcpy(input_parameters_name, input_file_path);
  strcat(input_parameters_name, "/pars.asc");
  hemelb::util::check_file(input_parameters_name);

  strcpy(vis_parameters_name, input_file_path);
  strcat(vis_parameters_name, "/rt_pars.asc");
  hemelb::util::check_file(vis_parameters_name);

  // Create directory path for the output images
  strcpy(image_directory, input_file_path);
  strcat(image_directory, "/Images/");

  //Create directory path for the output snapshots
  strcpy(snapshot_directory, input_file_path);
  strcat(snapshot_directory, "/Snapshots/");

  // Actually create the directories.
  if (lMaster.GetNet()->IsCurrentProcTheIOProc())
  {
    hemelb::util::MakeDirAllRXW(image_directory);
    hemelb::util::MakeDirAllRXW(snapshot_directory);

    char timings_name[256];
    char procs_string[256];

    sprintf(procs_string, "%i", lMaster.GetNet()->procs);
    strcpy(timings_name, input_file_path);
    strcat(timings_name, "/timings");
    strcat(timings_name, procs_string);
    strcat(timings_name, ".asc");

    timings_ptr = fopen(timings_name, "w");

    fprintf(timings_ptr, "***********************************************************\n");
    fprintf(timings_ptr, "Opening parameters file:\n %s\n", input_parameters_name);
    fprintf(timings_ptr, "Opening config file:\n %s\n", input_config_name);
    fprintf(timings_ptr, "Opening vis parameters file:\n %s\n\n", vis_parameters_name);

#ifndef NO_STEER
    hemelb::vis::xdrSendBuffer_pixel_data = new char[hemelb::vis::pixel_data_bytes];
    hemelb::vis::xdrSendBuffer_frame_details = new char[hemelb::vis::frame_details_bytes];

    pthread_mutex_init(&hemelb::steering::LOCK, NULL);
    pthread_cond_init(&hemelb::steering::network_send_frame, NULL);

    //    pthread_mutex_lock (&LOCK);

    pthread_attr_init(&pthread_attrib);
    pthread_attr_setdetachstate(&pthread_attrib, PTHREAD_CREATE_JOINABLE);

    pthread_create(&network_thread, &pthread_attrib, hemelb::steering::hemeLB_network,
                   (void*) lMaster.GetLBM());
#endif // NO_STEER
  }

  lMaster.GetLBM()->lbmInit(input_config_name, input_parameters_name, lMaster.GetNet());

  if (lMaster.GetNet()->netFindTopology(&depths) == 0)
  {
    fprintf(timings_ptr, "MPI_Attr_get failed, aborting\n");
    lMaster.GetNet()->Abort();
  }

  lMaster.GetNet()->netInit(lMaster.GetLBM()->total_fluid_sites);

  lMaster.GetLBM()->lbmSetInitialConditions(lMaster.GetNet());

  hemelb::vis::controller = new hemelb::vis::Control();
  hemelb::vis::controller->initLayers(lMaster.GetNet());

  hemelb::vis::controller->readParameters(vis_parameters_name, lMaster.GetLBM(),
                                          lMaster.GetNet());
#ifndef NO_STEER
  hemelb::steering::UpdateSteerableParameters(&hemelb::vis::doRendering,
                                              hemelb::vis::controller, lMaster.GetLBM());
#endif

  hemelb::util::DeleteDirContents(snapshot_directory);
  hemelb::util::DeleteDirContents(image_directory);

  total_time_steps = 0;

  int is_finished = 0;

  simulation_time = hemelb::util::myClock();

  if (lMaster.mSnapshotsPerCycle == 0)
    snapshots_period = 1e9;
  else
    snapshots_period = hemelb::util::max(1, lMaster.GetLBM()->period
        / lMaster.mSnapshotsPerCycle);

  if (lMaster.mImagesPerCycle == 0)
    images_period = 1e9;
  else
    images_period = hemelb::util::max(1, lMaster.GetLBM()->period / lMaster.mImagesPerCycle);

  for (cycle_id = 1; cycle_id <= lMaster.mMaxCycleCount && !is_finished; cycle_id++)
  {
    int restart = 0;

    lMaster.GetLBM()->lbmInitMinMaxValues();

    for (time_step = 1; time_step <= lMaster.GetLBM()->period; time_step++)
    {
      ++total_time_steps;
      intra_cycle_time = (PULSATILE_PERIOD * time_step) / lMaster.GetLBM()->period;

      int write_snapshot_image = (time_step % images_period == 0)
        ? 1
        : 0;
#ifndef NO_STEER
      int render_for_network_stream = 0;

      /* In the following two if blocks we do the core magic to ensure we only render
       when (1) we are not sending a frame or (2) we need to output to disk */

      if (lMaster.GetNet()->IsCurrentProcTheIOProc())
      {
        sem_wait(&hemelb::steering::connected_sem);
        bool local_connected = hemelb::steering::connected;
        sem_post(&hemelb::steering::connected_sem);
        if (local_connected)
        {
          render_for_network_stream = (hemelb::steering::sending_frame == 0)
            ? 1
            : 0;
        }
        else
        {
          render_for_network_stream = 0;
        }
      }

      if (total_time_steps % BCAST_FREQ == 0)
      {
        if (lMaster.GetNet()->IsCurrentProcTheIOProc())
        {
          hemelb::vis::doRendering = (render_for_network_stream || write_snapshot_image)
            ? 1
            : 0;
          sem_wait(&hemelb::steering::steering_var_lock);
        }
        hemelb::steering::UpdateSteerableParameters(&hemelb::vis::doRendering,
                                                    hemelb::vis::controller,
                                                    lMaster.GetLBM());
        if (lMaster.GetNet()->IsCurrentProcTheIOProc())
          sem_post(&hemelb::steering::steering_var_lock);
      }

      /* for debugging purposes we want to ensure we capture the variables in a single
       instant of time since variables might be altered by the thread half way through?
       This is to be done. */

      if (lMaster.GetNet()->IsCurrentProcTheIOProc() && time_step % 100 == 0)
        printf(
               "time step %i sending_frame %i render_network_stream %i write_snapshot_image %i rendering %i\n",
               time_step, hemelb::steering::sending_frame, render_for_network_stream,
               write_snapshot_image, hemelb::vis::doRendering);

#endif // NO_STEER
      lMaster.GetLBM()->lbmUpdateBoundaryDensities(cycle_id, time_step);

      stability = lMaster.GetLBM()->lbmCycle(hemelb::vis::doRendering, lMaster.GetNet());

      if ( (restart = lMaster.GetLBM()->IsUnstable(lMaster.GetNet())) != 0)
      {
        break;
      }
      lMaster.GetLBM()->lbmUpdateInletVelocities(time_step, lMaster.GetNet());

#ifndef NO_STREAKLINES
      hemelb::vis::controller->streaklines(time_step, lMaster.GetLBM()->period,
                                           lMaster.GetNet());
#endif
#ifndef NO_STEER
      if (total_time_steps % BCAST_FREQ == 0 && hemelb::vis::doRendering
          && !write_snapshot_image)
      {
        hemelb::vis::controller->render(RECV_BUFFER_A,
                                        hemelb::vis::ColourPalette::pickColour,
                                        lMaster.GetNet());

        if (hemelb::vis::controller->mouse_x >= 0 && hemelb::vis::controller->mouse_y
            >= 0 && hemelb::steering::updated_mouse_coords)
        {
          for (int i = 0; i < hemelb::vis::controller->col_pixels_recv[RECV_BUFFER_A]; i++)
          {
            if (hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.isRt
                && int(hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.i)
                    == hemelb::vis::controller->mouse_x
                && int(hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.j)
                    == hemelb::vis::controller->mouse_y)
            {
              hemelb::vis::controller->calculateMouseFlowField(
                                                               &hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i],
                                                               lMaster.GetLBM());
              break;
            }
          }
          hemelb::steering::updated_mouse_coords = 0;
        }
        if (lMaster.GetNet()->IsCurrentProcTheIOProc())
        {
          hemelb::steering::is_frame_ready = 1;
          sem_post(&hemelb::steering::nrl); // let go of the lock
        }
      }
#endif // NO_STEER
      if (write_snapshot_image)
      {
        hemelb::vis::controller->render(RECV_BUFFER_B,
                                        hemelb::vis::ColourPalette::pickColour,
                                        lMaster.GetNet());

        if (lMaster.GetNet()->IsCurrentProcTheIOProc())
        {
          char image_filename[255];

          snprintf(image_filename, 255, "%08i.dat", time_step);
          strcpy(complete_image_name, image_directory);
          strcat(complete_image_name, image_filename);

          hemelb::vis::controller->writeImage(RECV_BUFFER_B, complete_image_name,
                                              hemelb::vis::ColourPalette::pickColour);
        }
      }
      if (time_step % snapshots_period == 0)
      {
        char snapshot_filename[255];
        char complete_snapshot_name[255];

        snprintf(snapshot_filename, 255, "snapshot_%06i.asc", time_step);
        strcpy(complete_snapshot_name, snapshot_directory);
        strcat(complete_snapshot_name, snapshot_filename);

        lMaster.GetLBM()->lbmWriteConfig(stability, complete_snapshot_name,
                                         lMaster.GetNet());
      }
#ifndef NO_STEER
      if (lMaster.GetNet()->IsCurrentProcTheIOProc())
      {
        if (render_for_network_stream == 1)
        {
          // printf("sending signal to thread that frame is ready to go...\n"); fflush(0x0);
          sched_yield();
          sem_post(&hemelb::steering::nrl);
          //pthread_mutex_unlock (&LOCK);
          //pthread_cond_signal (&network_send_frame);
        }
      }
#endif
      if (stability == STABLE_AND_CONVERGED)
      {
        is_finished = 1;
        break;
      }
      if (lbm_terminate_simulation)
      {
        is_finished = 1;
        break;
      }
      if (lMaster.GetLBM()->period > 400000)
      {
        is_unstable = 1;
        break;
      }
    }

    if (restart)
    {
      hemelb::util::DeleteDirContents(snapshot_directory);
      hemelb::util::DeleteDirContents(image_directory);

      lMaster.GetLBM()->lbmRestart(lMaster.GetNet());
#ifndef NO_STREAKLINES
      hemelb::vis::controller->restart();
#endif
      if (lMaster.GetNet()->IsCurrentProcTheIOProc())
      {
        printf("restarting: period: %i\n", lMaster.GetLBM()->period);
        fflush(0x0);
      }
      snapshots_period = (lMaster.mSnapshotsPerCycle == 0)
        ? 1e9
        : hemelb::util::max(1, lMaster.GetLBM()->period / lMaster.mSnapshotsPerCycle);

      images_period = (lMaster.mImagesPerCycle == 0)
        ? 1e9
        : hemelb::util::max(1, lMaster.GetLBM()->period / lMaster.mImagesPerCycle);

      cycle_id = 0;
      continue;
    }
    lMaster.GetLBM()->lbmCalculateFlowFieldValues();

    if (lMaster.GetNet()->IsCurrentProcTheIOProc())
    {
      fprintf(timings_ptr, "cycle id: %i\n", cycle_id);
      printf("cycle id: %i\n", cycle_id);

      fflush(NULL);
    }
  }
  simulation_time = hemelb::util::myClock() - simulation_time;

  time_step = hemelb::util::min(time_step, lMaster.GetLBM()->period);
  cycle_id = hemelb::util::min(cycle_id, lMaster.mMaxCycleCount);
  time_step = time_step * cycle_id;

  if (lMaster.GetNet()->IsCurrentProcTheIOProc())
  {
    fprintf(timings_ptr, "\n");
    fprintf(timings_ptr, "threads: %i, machines checked: %i\n\n",
            lMaster.GetNet()->procs, lMaster.GetNet()->GetMachineCount());
    fprintf(timings_ptr, "topology depths checked: %i\n\n", depths);
    fprintf(timings_ptr, "fluid sites: %i\n\n", lMaster.GetLBM()->total_fluid_sites);
    fprintf(timings_ptr, "cycles and total time steps: %i, %i \n\n", cycle_id,
            total_time_steps);
    fprintf(timings_ptr, "time steps per second: %.3f\n\n", total_time_steps
        / simulation_time);
  }

  if (is_unstable)
  {
    if (lMaster.GetNet()->IsCurrentProcTheIOProc())
    {
      fprintf(timings_ptr, "Attention: simulation unstable with %i timesteps/cycle\n",
              lMaster.GetLBM()->period);
      fprintf(timings_ptr, "Simulation is terminated\n");
      fclose(timings_ptr);
    }
  }
  else
  {
    if (lMaster.GetNet()->IsCurrentProcTheIOProc())
    {

      fprintf(timings_ptr, "time steps per cycle: %i\n", lMaster.GetLBM()->period);
      fprintf(timings_ptr, "pressure min, max (mmHg): %e, %e\n",
              lMaster.GetLBM()->GetMinPhysicalPressure(),
              lMaster.GetLBM()->GetMaxPhysicalPressure());
      fprintf(timings_ptr, "velocity min, max (m/s) : %e, %e\n",
              lMaster.GetLBM()->GetMinPhysicalVelocity(),
              lMaster.GetLBM()->GetMaxPhysicalVelocity());
      fprintf(timings_ptr, "stress   min, max (Pa)  : %e, %e\n",
              lMaster.GetLBM()->GetMinPhysicalStress(),
              lMaster.GetLBM()->GetMaxPhysicalStress());
      fprintf(timings_ptr, "\n");

      for (int n = 0; n < lMaster.GetLBM()->inlets; n++)
      {
        fprintf(timings_ptr, "inlet id: %i, average / peak velocity (m/s): %e / %e\n", n,
                lMaster.GetLBM()->GetAverageInletVelocity(n),
                lMaster.GetLBM()->GetPeakInletVelocity(n));
      }
      fprintf(timings_ptr, "\n");

      fprintf(timings_ptr, "\n");
      fprintf(timings_ptr, "domain decomposition time (s):             %.3f\n",
              lMaster.GetNet()->dd_time);
      fprintf(timings_ptr, "pre-processing buffer management time (s): %.3f\n",
              lMaster.GetNet()->bm_time);
      fprintf(timings_ptr, "input configuration reading time (s):      %.3f\n",
              lMaster.GetNet()->fr_time);

      total_time = hemelb::util::myClock() - total_time;
      fprintf(timings_ptr, "total time (s):                            %.3f\n\n",
              total_time);

      fprintf(timings_ptr, "Sub-domains info:\n\n");

      for (int n = 0; n < lMaster.GetNet()->procs; n++)
      {
        fprintf(timings_ptr, "rank: %i, fluid sites: %i\n", n,
                lMaster.GetNet()->fluid_sites[n]);
      }

      fclose(timings_ptr);
    }
  }
  delete hemelb::vis::controller;

#ifndef NO_STEER
  if (lMaster.GetNet()->IsCurrentProcTheIOProc())
  {
    // there are some problems if the following function is called

    //pthread_join (network_thread, NULL);
    delete[] hemelb::vis::xdrSendBuffer_frame_details;
    delete[] hemelb::vis::xdrSendBuffer_pixel_data;
  }
#endif

  return (0);
}

