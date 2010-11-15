#include "SimulationMaster.h"
#include "SimConfig.h"

#include "vis/visthread.h"
#include "vis/ColourPalette.h"
#include "utilityFunctions.h"
#include "usage.h"
#include "debug/Debugger.h"
#include "fileutils.h"

#include <limits>
#include <stdlib.h>

SimulationMaster::SimulationMaster(int iArgCount, char *iArgList[])
{
  mLbm = new LBM();
  mNet = new Net(iArgCount, iArgList);

  hemelb::debug::Debugger::Init(iArgList[0]);

  mLbTime = 0.0;
  mMPISendTime = 0.0;
  mMPIWaitTime = 0.0;
  mImagingTime = 0.0;
  mSnapshotTime = 0.0;

  mImagesWritten = 0;
  mSnapshotsWritten = 0;
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
  mTimingsFile = bTimingsFile;
  mSimulationState.IsTerminating = 0;

  steeringController
      = hemelb::steering::Control::Init(GetNet()->IsCurrentProcTheIOProc());

  if (GetNet()->IsCurrentProcTheIOProc())
  {
    steeringController->StartNetworkThread(GetLBM(), &mSimulationState,
                                           GetLBM()->GetLbmParams());
  }

  GetLBM()->lbmInit(iSimConfig, mGlobLatDat, iSteeringSessionid,
                    (int) (iSimConfig->StepsPerCycle), iSimConfig->VoxelSize,
                    GetNet());
  if (GetNet()->netFindTopology() == 0)
  {
    fprintf(bTimingsFile, "MPI_Attr_get failed, aborting\n");
    GetNet()->Abort();
  }
  GetNet()->Initialise(GetLBM()->total_fluid_sites, mGlobLatDat, mLocalLatDat);
  GetLBM()->lbmSetInitialConditions(GetNet(), mLocalLatDat);
  hemelb::vis::controller
      = new hemelb::vis::Control(mLbm->GetLbmParams()->StressType, mGlobLatDat);
  hemelb::vis::controller->initLayers(mGlobLatDat, mLocalLatDat, GetNet());
  GetLBM()->ReadVisParameters(GetNet());

  hemelb::vis::controller->SetProjection(512, 512, iSimConfig->VisCentre.x,
                                         iSimConfig->VisCentre.y,
                                         iSimConfig->VisCentre.z,
                                         iSimConfig->VisLongitude,
                                         iSimConfig->VisLatitude,
                                         iSimConfig->VisZoom);

  steeringController->UpdateSteerableParameters(false,
                                                &hemelb::vis::doRendering,
                                                mSimulationState,
                                                hemelb::vis::controller,
                                                GetLBM());
}

void SimulationMaster::RunSimulation(hemelb::SimConfig *& lSimulationConfig,
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
                                                      mSimulationState,
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

      stability = GetLBM()->lbmCycle(hemelb::vis::doRendering, GetNet(),
                                     mLocalLatDat, mLbTime, mMPISendTime,
                                     mMPIWaitTime);

      if ( (restart = GetLBM()->IsUnstable(mLocalLatDat, GetNet())) != false)
      {
        break;
      }
      GetLBM()->lbmUpdateInletVelocities(mSimulationState.TimeStep,
                                         mLocalLatDat, GetNet());

#ifndef NOMPI
      double lPreImageTime = MPI_Wtime();
#endif

#ifndef NO_STREAKLINES
      hemelb::vis::controller->streaklines(mSimulationState.TimeStep,
                                           GetLBM()->period, mGlobLatDat,
                                           mLocalLatDat, GetNet());
#endif
#ifndef NO_STEER

      if (total_time_steps % BCAST_FREQ == 0 && hemelb::vis::doRendering
          && !write_snapshot_image)
      {
        hemelb::vis::controller->render(RECV_BUFFER_A, mGlobLatDat, GetNet());

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
        hemelb::vis::controller->render(RECV_BUFFER_B, mGlobLatDat, GetNet());
        mImagesWritten++;

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

#ifndef NOMPI
      double lPreSnapshotTime = MPI_Wtime();
      mImagingTime += (lPreSnapshotTime - lPreImageTime);
#endif

      if (mSimulationState.TimeStep % snapshots_period == 0)
      {
        char snapshot_filename[255];
        snprintf(snapshot_filename, 255, "snapshot_%06i.asc",
                 mSimulationState.TimeStep);

        mSnapshotsWritten++;
        GetLBM()->lbmWriteConfig(stability, snapshot_directory
            + std::string(snapshot_filename), GetNet(), mGlobLatDat,
                                 mLocalLatDat);
      }

#ifndef NOMPI
      mSnapshotTime += (MPI_Wtime() - lPreSnapshotTime);
#endif

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
      if (mSimulationState.IsTerminating)
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

      GetLBM()->lbmRestart(mLocalLatDat, GetNet());
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
      fprintf(mTimingsFile, "cycle id: %i\n", mSimulationState.CycleId);
      printf("cycle id: %i\n", mSimulationState.CycleId);

      fflush(NULL);
    }
  }

  hemelb::debug::Debugger::Get()->BreakHere();

  mSimulationState.CycleId = hemelb::util::min(mSimulationState.CycleId,
                                               lSimulationConfig->NumCycles);

  PostSimulation(total_time_steps, hemelb::util::myClock() - simulation_time,
                 is_unstable, iStartTime);
}

void SimulationMaster::PostSimulation(int iTotalTimeSteps,
                                      double iSimulationTime,
                                      bool iIsUnstable,
                                      double iStartTime)
{
  if (GetNet()->IsCurrentProcTheIOProc())
  {
    fprintf(mTimingsFile, "\n");
    fprintf(mTimingsFile, "threads: %i, machines checked: %i\n\n",
            GetNet()->mProcessorCount, GetNet()->GetMachineCount());
    fprintf(mTimingsFile, "topology depths checked: %i\n\n", mNet->GetDepths());
    fprintf(mTimingsFile, "fluid sites: %i\n\n", GetLBM()->total_fluid_sites);
    fprintf(mTimingsFile, "cycles and total time steps: %i, %i \n\n",
            mSimulationState.CycleId, iTotalTimeSteps);
    fprintf(mTimingsFile, "time steps per second: %.3f\n\n", iTotalTimeSteps
        / iSimulationTime);
  }

  if (iIsUnstable)
  {
    if (GetNet()->IsCurrentProcTheIOProc())
    {
      fprintf(mTimingsFile,
              "Attention: simulation unstable with %i timesteps/cycle\n",
              GetLBM()->period);
      fprintf(mTimingsFile, "Simulation is terminated\n");
    }
  }
  else
  {
    if (GetNet()->IsCurrentProcTheIOProc())
    {

      fprintf(mTimingsFile, "time steps per cycle: %i\n", GetLBM()->period);
      fprintf(mTimingsFile, "pressure min, max (mmHg): %e, %e\n",
              GetLBM()->GetMinPhysicalPressure(),
              GetLBM()->GetMaxPhysicalPressure());
      fprintf(mTimingsFile, "velocity min, max (m/s) : %e, %e\n",
              GetLBM()->GetMinPhysicalVelocity(),
              GetLBM()->GetMaxPhysicalVelocity());
      fprintf(mTimingsFile, "stress   min, max (Pa)  : %e, %e\n",
              GetLBM()->GetMinPhysicalStress(),
              GetLBM()->GetMaxPhysicalStress());
      fprintf(mTimingsFile, "\n");

      for (int n = 0; n < GetLBM()->inlets; n++)
      {
        fprintf(mTimingsFile,
                "inlet id: %i, average / peak velocity (m/s): %e / %e\n", n,
                GetLBM()->GetAverageInletVelocity(n),
                GetLBM()->GetPeakInletVelocity(n));
      }
      fprintf(mTimingsFile, "\n");

      fprintf(mTimingsFile, "\n");
      fprintf(mTimingsFile,
              "domain decomposition time (s):             %.3f\n",
              GetNet()->dd_time);
      fprintf(mTimingsFile,
              "pre-processing buffer management time (s): %.3f\n",
              GetNet()->bm_time);
      fprintf(mTimingsFile,
              "input configuration reading time (s):      %.3f\n",
              GetNet()->fr_time);

      fprintf(mTimingsFile,
              "total time (s):                            %.3f\n\n",
               (hemelb::util::myClock() - iStartTime));

      fprintf(mTimingsFile, "Sub-domains info:\n\n");

      for (int n = 0; n < GetNet()->mProcessorCount; n++)
      {
        fprintf(mTimingsFile, "rank: %i, fluid sites: %i\n", n,
                GetNet()->mFluidSitesOnEachProcessor[n]);
      }
    }
  }

  PrintTimingData(-1);
}

// NB. This function takes into account that the process numbered 0
// does not always carry data.
void SimulationMaster::PrintTimingData(int iSignal)
{
#ifndef NOMPI
  double lTimings[5] = { mLbTime, mMPISendTime, mMPIWaitTime, mImagingTime,
                         mSnapshotTime };
  std::string lNames[5] = { "LBM", "MPISend", "MPIWait", "Images", "Snaps" };

  if (mSimulationState.CycleId > 0)
  {
    mLbTime /= (double) mSimulationState.CycleId;
    mMPISendTime /= (double) mSimulationState.CycleId;
    mMPIWaitTime /= (double) mSimulationState.CycleId;
  }

  if (mImagesWritten > 0)
  {
    mImagingTime /= (double) mImagesWritten;
  }

  if (mSnapshotsWritten > 0)
  {
    mSnapshotTime /= (double) mSnapshotsWritten;
  }

  double lMins[5];
  double lMaxes[5];
  double lMeans[5];

  if (mNet->IsCurrentProcTheIOProc())
  {
    for (int ii = 0; ii < 3; ii++)
      lTimings[ii] = 0.0;
  }

  MPI_Reduce(lTimings, lMaxes, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(lTimings, lMeans, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // Change the values for LBM and MPI on process 0 so they don't interfere with the min
  // operation (previously values were 0.0 so they won't affect max / mean
  // calc).
  if (mNet->IsCurrentProcTheIOProc())
  {
    for (int ii = 0; ii < 3; ii++)
      lTimings[ii] = std::numeric_limits<double>::max();
  }
  MPI_Reduce(lTimings, lMins, 5, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if (mNet->IsCurrentProcTheIOProc())
  {
    if (mNet->mProcessorCount > 1)
    {
      for (int ii = 0; ii < 3; ii++)
        lMeans[ii] /= (double) (mNet->mProcessorCount - 1);
      for (int ii = 3; ii < 5; ii++)
        lMeans[ii] /= (double) mNet->mProcessorCount;
    }

    fprintf(mTimingsFile,
            "\n\nPer-proc timing data (secs per [cycle,image,snapshot]): \n\n");
    fprintf(mTimingsFile, "\t\tMin \t\tMean \t\tMax\n");
    for (int ii = 0; ii < 5; ii++)
    {
      fprintf(mTimingsFile, "%s \t\t%.3g \t\t%.3g \t\t%.3g\n",
              lNames[ii].c_str(), lMins[ii], lMeans[ii], lMaxes[ii]);
    }
  }
#endif //NOMPI
}

LBM *SimulationMaster::GetLBM()
{
  return mLbm;
}

Net *SimulationMaster::GetNet()
{
  return mNet;
}
