#include "SimulationMaster.h"
#include "SimConfig.h"

#include "vis/visthread.h"
#include "vis/ColourPalette.h"
#include "utilityFunctions.h"
#include "usage.h"
#include "debug/Debugger.h"
#include "fileutils.h"

#include <stdlib.h>

SimulationMaster::SimulationMaster(int iArgCount, char *iArgList[])
{
  mLbm = new LBM();
  mNet = new Net(iArgCount, iArgList);

  hemelb::debug::Debugger::Init(iArgList[0]);
}

SimulationMaster::~SimulationMaster()
{
  delete mNet;
  delete mLbm;
}

void SimulationMaster::Initialise(hemelb::SimConfig *iSimConfig,
                                  int iSteeringSessionid,
                                  FILE *bTimingsFile)
{
  steeringController
      = hemelb::steering::Control::Init(GetNet()->IsCurrentProcTheIOProc());

  if (GetNet()->IsCurrentProcTheIOProc())
  {
    steeringController->StartNetworkThread(GetLBM(), &mSimulationState);
  }

  GetLBM()->lbmInit(iSimConfig, iSteeringSessionid,
                    (int) (iSimConfig->StepsPerCycle), iSimConfig->VoxelSize,
                    GetNet());
  if (GetNet()->netFindTopology() == 0)
  {
    fprintf(bTimingsFile, "MPI_Attr_get failed, aborting\n");
    GetNet()->Abort();
  }
  GetNet()->Initialise(GetLBM()->total_fluid_sites);
  GetLBM()->lbmSetInitialConditions(GetNet());
  hemelb::vis::controller = new hemelb::vis::Control(lbm_stress_type);
  hemelb::vis::controller->initLayers(GetNet());
  GetLBM()->ReadVisParameters(GetNet(), hemelb::vis::controller);
  steeringController->UpdateSteerableParameters(false,
                                                &hemelb::vis::doRendering,
                                                hemelb::vis::controller,
                                                GetLBM());
}

void SimulationMaster::RunSimulation(FILE *iTimingsFile,
                                     hemelb::SimConfig *& lSimulationConfig,
                                     double iStartTime,
                                     std::string image_directory,
                                     std::string snapshot_directory,
                                     unsigned int lSnapshotsPerCycle,
                                     unsigned int lImagesPerCycle)
{
  double simulation_time = hemelb::util::myClock();
  bool is_unstable = false;
  int total_time_steps = 0;

  int snapshots_period = (lSnapshotsPerCycle == 0)
    ? 1e9
    : hemelb::util::max(1, GetLBM()->period / lSnapshotsPerCycle);

  int images_period = (lImagesPerCycle == 0)
    ? 1e9
    : hemelb::util::max(1, GetLBM()->period / lImagesPerCycle);

  bool is_finished = false;
  int stability = STABLE;

  for (mSimulationState.CycleId = 1; mSimulationState.CycleId
      <= lSimulationConfig->NumCycles && !is_finished; mSimulationState.CycleId++)
  {
    bool restart = false;

    GetLBM()->lbmInitMinMaxValues();

    for (mSimulationState.TimeStep = 1; mSimulationState.TimeStep
        <= GetLBM()->period; mSimulationState.TimeStep++)
    {
      ++total_time_steps;
      mSimulationState.IntraCycleTime = (PULSATILE_PERIOD
          * mSimulationState.TimeStep) / GetLBM()->period;

      bool write_snapshot_image = ( (mSimulationState.TimeStep % images_period)
          == 0)
        ? true
        : false;

      /* In the following two if blocks we do the core magic to ensure we only Render
       when (1) we are not sending a frame or (2) we need to output to disk */

      int render_for_network_stream = 0;
      if (GetNet()->IsCurrentProcTheIOProc())
      {
        render_for_network_stream
            = steeringController->ShouldRenderForNetwork();
      }

      if (mSimulationState.TimeStep % BCAST_FREQ == 0)
      {
        steeringController->UpdateSteerableParameters(
                                                      write_snapshot_image,
                                                      &hemelb::vis::doRendering,
                                                      hemelb::vis::controller,
                                                      GetLBM());

      }

      /* for debugging purposes we want to ensure we capture the variables in a single
       instant of time since variables might be altered by the thread half way through?
       This is to be done. */

      if (GetNet()->IsCurrentProcTheIOProc() && mSimulationState.TimeStep % 100
          == 0)
        printf(
               "time step %i sending_frame %i render_network_stream %i write_snapshot_image %i rendering %i\n",
               mSimulationState.TimeStep, steeringController->sending_frame,
               render_for_network_stream, write_snapshot_image,
               hemelb::vis::doRendering);

      GetLBM()->lbmUpdateBoundaryDensities(mSimulationState.CycleId,
                                           mSimulationState.TimeStep);

      stability = GetLBM()->lbmCycle(hemelb::vis::doRendering, GetNet());

      if ( (restart = GetLBM()->IsUnstable(GetNet())) != false)
      {
        break;
      }
      GetLBM()->lbmUpdateInletVelocities(mSimulationState.TimeStep, GetNet());

#ifndef NO_STREAKLINES
      hemelb::vis::controller->streaklines(mSimulationState.TimeStep,
                                           GetLBM()->period, GetNet());
#endif
#ifndef NO_STEER

      if (total_time_steps % BCAST_FREQ == 0 && hemelb::vis::doRendering
          && !write_snapshot_image)
      {
        hemelb::vis::controller->render(RECV_BUFFER_A, GetNet());

        if (hemelb::vis::controller->mouse_x >= 0
            && hemelb::vis::controller->mouse_y >= 0
            && steeringController->updated_mouse_coords)
        {
          for (int i = 0; i
              < hemelb::vis::controller->col_pixels_recv[RECV_BUFFER_A]; i++)
          {
            if (hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.isRt
                && int(
                       hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.i)
                    == hemelb::vis::controller->mouse_x
                && int(
                       hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.j)
                    == hemelb::vis::controller->mouse_y)
            {
              double mouse_pressure, mouse_stress;
              GetLBM()->CalculateMouseFlowField(
                                                &hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i],
                                                mouse_pressure,
                                                mouse_stress,
                                                hemelb::vis::controller->density_threshold_min,
                                                hemelb::vis::controller->density_threshold_minmax_inv,
                                                hemelb::vis::controller->stress_threshold_max_inv);

              hemelb::vis::controller->setMouseParams(mouse_pressure,
                                                      mouse_stress);

              break;
            }
          }
          steeringController->updated_mouse_coords = 0;
        }
        if (GetNet()->IsCurrentProcTheIOProc())
        {
          steeringController->is_frame_ready = 1;
          sem_post(&steeringController->nrl); // let go of the lock
        }
      }
#endif // NO_STEER
      if (write_snapshot_image)
      {
        hemelb::vis::controller->render(RECV_BUFFER_B, GetNet());

        if (GetNet()->IsCurrentProcTheIOProc())
        {
          char image_filename[255];
          snprintf(image_filename, 255, "%08i.dat", mSimulationState.TimeStep);

          hemelb::vis::controller->writeImage(
                                              RECV_BUFFER_B,
                                              image_directory
                                                  + std::string(image_filename),
                                              hemelb::vis::ColourPalette::pickColour);
        }
      }
      if (mSimulationState.TimeStep % snapshots_period == 0)
      {
        char snapshot_filename[255];
        snprintf(snapshot_filename, 255, "snapshot_%06i.asc",
                 mSimulationState.TimeStep);

        GetLBM()->lbmWriteConfig(stability, snapshot_directory
            + std::string(snapshot_filename), GetNet());
      }
#ifndef NO_STEER
      if (GetNet()->IsCurrentProcTheIOProc())
      {
        if (render_for_network_stream == 1)
        {
          // printf("sending signal to thread that frame is ready to go...\n"); fflush(0x0);
          sched_yield();
          sem_post(&steeringController->nrl);
          //pthread_mutex_unlock (&LOCK);
          //pthread_cond_signal (&network_send_frame);
        }
      }
#endif
      if (stability == STABLE_AND_CONVERGED)
      {
        is_finished = true;
        break;
      }
      if (lbm_terminate_simulation)
      {
        is_finished = true;
        break;
      }
      if (GetLBM()->period > 400000)
      {
        is_unstable = true;
        break;
      }
    }

    if (restart)
    {
      hemelb::util::DeleteDirContents(snapshot_directory);
      hemelb::util::DeleteDirContents(image_directory);

      GetLBM()->lbmRestart(GetNet());
#ifndef NO_STREAKLINES
      hemelb::vis::controller->restart();
#endif
      if (GetNet()->IsCurrentProcTheIOProc())
      {
        printf("restarting: period: %i\n", GetLBM()->period);
        fflush(0x0);
      }
      snapshots_period = (lSnapshotsPerCycle == 0)
        ? 1e9
        : hemelb::util::max(1, GetLBM()->period / lSnapshotsPerCycle);

      images_period = (lImagesPerCycle == 0)
        ? 1e9
        : hemelb::util::max(1, GetLBM()->period / lImagesPerCycle);

      mSimulationState.CycleId = 0;
      continue;
    }
    GetLBM()->lbmCalculateFlowFieldValues();

    if (GetNet()->IsCurrentProcTheIOProc())
    {
      fprintf(iTimingsFile, "cycle id: %i\n", mSimulationState.CycleId);
      printf("cycle id: %i\n", mSimulationState.CycleId);

      fflush(NULL);
    }
  }

  mSimulationState.CycleId = hemelb::util::min(mSimulationState.CycleId,
                                               lSimulationConfig->NumCycles);

  PostSimulation(total_time_steps, hemelb::util::myClock() - simulation_time,
                 iTimingsFile, is_unstable, iStartTime);
}

void SimulationMaster::PostSimulation(int iTotalTimeSteps,
                                      double iSimulationTime,
                                      FILE *timings_ptr,
                                      bool iIsUnstable,
                                      double iStartTime)
{
  if (GetNet()->IsCurrentProcTheIOProc())
  {
    fprintf(timings_ptr, "\n");
    fprintf(timings_ptr, "threads: %i, machines checked: %i\n\n",
            GetNet()->mProcessorCount, GetNet()->GetMachineCount());
    fprintf(timings_ptr, "topology depths checked: %i\n\n", mNet->GetDepths());
    fprintf(timings_ptr, "fluid sites: %i\n\n", GetLBM()->total_fluid_sites);
    fprintf(timings_ptr, "cycles and total time steps: %i, %i \n\n",
            mSimulationState.CycleId, iTotalTimeSteps);
    fprintf(timings_ptr, "time steps per second: %.3f\n\n", iTotalTimeSteps
        / iSimulationTime);
  }

  if (iIsUnstable)
  {
    if (GetNet()->IsCurrentProcTheIOProc())
    {
      fprintf(timings_ptr,
              "Attention: simulation unstable with %i timesteps/cycle\n",
              GetLBM()->period);
      fprintf(timings_ptr, "Simulation is terminated\n");
    }
  }
  else
  {
    if (GetNet()->IsCurrentProcTheIOProc())
    {

      fprintf(timings_ptr, "time steps per cycle: %i\n", GetLBM()->period);
      fprintf(timings_ptr, "pressure min, max (mmHg): %e, %e\n",
              GetLBM()->GetMinPhysicalPressure(),
              GetLBM()->GetMaxPhysicalPressure());
      fprintf(timings_ptr, "velocity min, max (m/s) : %e, %e\n",
              GetLBM()->GetMinPhysicalVelocity(),
              GetLBM()->GetMaxPhysicalVelocity());
      fprintf(timings_ptr, "stress   min, max (Pa)  : %e, %e\n",
              GetLBM()->GetMinPhysicalStress(),
              GetLBM()->GetMaxPhysicalStress());
      fprintf(timings_ptr, "\n");

      for (int n = 0; n < GetLBM()->inlets; n++)
      {
        fprintf(timings_ptr,
                "inlet id: %i, average / peak velocity (m/s): %e / %e\n", n,
                GetLBM()->GetAverageInletVelocity(n),
                GetLBM()->GetPeakInletVelocity(n));
      }
      fprintf(timings_ptr, "\n");

      fprintf(timings_ptr, "\n");
      fprintf(timings_ptr, "domain decomposition time (s):             %.3f\n",
              GetNet()->dd_time);
      fprintf(timings_ptr, "pre-processing buffer management time (s): %.3f\n",
              GetNet()->bm_time);
      fprintf(timings_ptr, "input configuration reading time (s):      %.3f\n",
              GetNet()->fr_time);

      fprintf(timings_ptr,
              "total time (s):                            %.3f\n\n",
               (hemelb::util::myClock() - iStartTime));

      fprintf(timings_ptr, "Sub-domains info:\n\n");

      for (int n = 0; n < GetNet()->mProcessorCount; n++)
      {
        fprintf(timings_ptr, "rank: %i, fluid sites: %i\n", n,
                GetNet()->mFluidSitesOnEachProcessor[n]);
      }
    }
  }
}

LBM *SimulationMaster::GetLBM()
{
  return mLbm;
}

Net *SimulationMaster::GetNet()
{
  return mNet;
}
