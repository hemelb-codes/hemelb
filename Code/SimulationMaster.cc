#include "SimulationMaster.h"
#include "SimConfig.h"

#include "vis/ColourPalette.h"
#include "util/utilityFunctions.h"
#include "topology/TopologyReader.h"
#include "debug/Debugger.h"
#include "util/fileutils.h"

#include <limits>
#include <stdlib.h>

/**
 * Constructor for the SimulationMaster class
 *
 * Initialises member variables including the network topology
 * object.
 */
SimulationMaster::SimulationMaster(int iArgCount, char *iArgList[])
{
  // Initialise the network discovery. If this fails, abort.
  bool lTopologySuccess = true;
  mNetworkTopology
      = new hemelb::topology::NetworkTopology(&iArgCount, &iArgList, &lTopologySuccess);

  if (!lTopologySuccess)
  {
    printf("Couldn't get machine information for this network topology. Aborting.\n");
    Abort();
  }

  hemelb::debug::Debugger::Init(iArgList[0]);

  mLbm = NULL;
  mNet = NULL;
  mLocalLatDat = NULL;
  steeringCpt = NULL;

  mSimulationState.IsTerminating = 0;
  mSimulationState.DoRendering = 0;

  mLbTime = 0.0;
  mMPISendTime = 0.0;
  mMPIWaitTime = 0.0;
  mImagingTime = 0.0;
  mSnapshotTime = 0.0;

  mImagesWritten = 0;
  mSnapshotsWritten = 0;
}

/**
 * Destructor for the SimulationMaster class.
 *
 * Deallocates dynamically allocated memory to contained classes.
 */
SimulationMaster::~SimulationMaster()
{
  delete mNetworkTopology;

  if (mNet != NULL)
  {
    delete mNet;
  }
  if (mLbm != NULL)
  {
    delete mLbm;
  }
  if (mLocalLatDat != NULL)
  {
    delete mLocalLatDat;
  }

  if (steeringCpt != NULL)
  {
    delete steeringCpt;
  }
}

/**
 * Returns true if the current processor is the dedicated I/O
 * processor.
 */
bool SimulationMaster::IsCurrentProcTheIOProc()
{
  return mNetworkTopology->IsCurrentProcTheIOProc();
}

/**
 * Returns the number of processors involved in the simulation.
 */
int SimulationMaster::GetProcessorCount()
{
  return mNetworkTopology->GetProcessorCount();
}

/**
 * Initialises various elements of the simulation if necessary - steering,
 * domain decomposition, LBM and visualisation.
 */
void SimulationMaster::Initialise(hemelb::SimConfig *iSimConfig,
                                  unsigned int iImagesPerCycle,
                                  int iSteeringSessionid,
                                  FILE * bTimingsFile)
{
  mTimingsFile = bTimingsFile;

  // Initialise the Lbm.
  mLbm = new hemelb::lb::LBM(iSimConfig, mNetworkTopology);

  hemelb::topology::TopologyReader lTopologist;
  lTopologist.LoadAndDecompose(&mGlobLatDat, mLbm->total_fluid_sites, mLbm->siteMins,
                               mLbm->siteMaxes,
                               hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                               mNetworkTopology, mLbm->GetLbmParams(), iSimConfig, &mFileReadTime,
                               &mDomainDecompTime);

  // Initialise and begin the steering.
  if (mNetworkTopology->IsCurrentProcTheIOProc())
  {
    clientConnection = new hemelb::steering::ClientConnection(iSteeringSessionid);
  }
  else
  {
    clientConnection = NULL;
  }

  if (mNetworkTopology->IsCurrentProcTheIOProc())
  {
    imageSendCpt = new hemelb::steering::ImageSendComponent(mLbm, &mSimulationState,
                                                            mLbm->GetLbmParams(), clientConnection);
  }

  // Initialise the Net object and the Lbm.
  mLocalLatDat
      = new hemelb::lb::LocalLatticeData(
                                         mNetworkTopology->FluidSitesOnEachProcessor[mNetworkTopology->GetLocalRank()]);

  mNet = new hemelb::net::Net(mNetworkTopology);

  mStabilityTester = new hemelb::lb::StabilityTester(mLocalLatDat, mNet, mNetworkTopology,
                                                     &mSimulationState);

  double seconds = hemelb::util::myClock();
  int* lReceiveTranslator = mNet->Initialise(mGlobLatDat, mLocalLatDat);
  mNetInitialiseTime = hemelb::util::myClock() - seconds;

  mLbm->Initialise(lReceiveTranslator, mLocalLatDat);

  // Initialise the visualisation controller.
  hemelb::vis::controller
      = new hemelb::vis::Control(mLbm->GetLbmParams()->StressType, &mGlobLatDat);
  hemelb::vis::controller->initLayers(mNetworkTopology, &mGlobLatDat, mLocalLatDat);

  int images_period = (iImagesPerCycle == 0)
    ? 1e9
    : hemelb::util::NumericalFunctions::max<unsigned int>(1U, mLbm->period / iImagesPerCycle);

  steeringCpt = new hemelb::steering::SteeringComponent(images_period, clientConnection,
                                                        hemelb::vis::controller, mLbm, mNet,
                                                        mNetworkTopology, &mSimulationState);

  // Read in the visualisation parameters.
  mLbm->ReadVisParameters();

  hemelb::vis::controller->SetProjection(512, 512, iSimConfig->VisCentre.x,
                                         iSimConfig->VisCentre.y, iSimConfig->VisCentre.z,
                                         iSimConfig->VisLongitude, iSimConfig->VisLatitude,
                                         iSimConfig->VisZoom);
}

/**
 * Begin the simulation.
 */
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

  unsigned int snapshots_period = (lSnapshotsPerCycle == 0)
    ? 1e9
    : hemelb::util::NumericalFunctions::max<unsigned int>(1U, mLbm->period / lSnapshotsPerCycle);

  unsigned int images_period = (lImagesPerCycle == 0)
    ? 1e9
    : hemelb::util::NumericalFunctions::max<unsigned int>(1U, mLbm->period / lImagesPerCycle);

  bool is_finished = false;
  hemelb::lb::Stability stability = hemelb::lb::Stable;

  for (mSimulationState.CycleId = 1; mSimulationState.CycleId <= lSimulationConfig->NumCycles
      && !is_finished; mSimulationState.CycleId++)
  {
    bool restart = false;

    mLbm->InitMinMaxValues();

    for (mSimulationState.TimeStep = 1; mSimulationState.TimeStep <= mLbm->period; mSimulationState.TimeStep++)
    {
      ++total_time_steps;
      mSimulationState.IntraCycleTime = (PULSATILE_PERIOD * mSimulationState.TimeStep)
          / mLbm->period;

      bool write_snapshot_image = ( (mSimulationState.TimeStep % images_period) == 0)
        ? true
        : false;

      // Make sure we're rendering if we're writing this iteration.
      if (write_snapshot_image)
      {
        mSimulationState.DoRendering = true;
      }

      /* In the following two if blocks we do the core magic to ensure we only Render
       when (1) we are not sending a frame or (2) we need to output to disk */

      bool render_for_network_stream = false;
      if (mNetworkTopology->IsCurrentProcTheIOProc())
      {
        render_for_network_stream = imageSendCpt->ShouldRenderNewNetworkImage();
      }

      /* for debugging purposes we want to ensure we capture the variables in a single
       instant of time since variables might be altered by the thread half way through?
       This is to be done. */

      if (mNetworkTopology->IsCurrentProcTheIOProc() && mSimulationState.TimeStep % 100 == 0)
        printf("time step %i render_network_stream %i write_snapshot_image %i rendering %i\n",
               mSimulationState.TimeStep, render_for_network_stream, write_snapshot_image,
               mSimulationState.DoRendering);

      mLbm->UpdateBoundaryDensities(mSimulationState.CycleId, mSimulationState.TimeStep);

      // Cycle.
      {
        mLbm->RequestComms(mNet, mLocalLatDat);
        steeringCpt->RequestComms();
        mStabilityTester->RequestComms();

        mNet->Receive();

        {
          double lPrePreSend = MPI_Wtime();
          mLbm->PreSend(mLocalLatDat, mSimulationState.DoRendering);
          mLbTime += (MPI_Wtime() - lPrePreSend);
        }

        {
          double lPreSendTime = MPI_Wtime();
          mNet->Send();
          mMPISendTime += (MPI_Wtime() - lPreSendTime);
        }

        {
          double lPrePreReceive = MPI_Wtime();
          mLbm->PreReceive(mSimulationState.DoRendering, mLocalLatDat);
          mLbTime += (MPI_Wtime() - lPrePreReceive);
        }

        {
          double lPreWaitTime = MPI_Wtime();
          mNet->Wait(mLocalLatDat);
          mMPIWaitTime += (MPI_Wtime() - lPreWaitTime);
        }

        {
          double lPrePostStep = MPI_Wtime();
          mLbm->PostReceive(mLocalLatDat, mSimulationState.DoRendering);
          mLbTime += (MPI_Wtime() - lPrePostStep);
        }

        mStabilityTester->PostReceive();
        steeringCpt->PostReceive();

        mLbm->EndIteration(mLocalLatDat);
      }

      stability = hemelb::lb::Stable;

      restart = (mSimulationState.Stability == hemelb::lb::Unstable);
      if (restart)
      {
        break;
      }
      mLbm->UpdateInletVelocities(mSimulationState.TimeStep, *mLocalLatDat, mNet);

      double lPreImageTime = MPI_Wtime();

#ifndef NO_STREAKLINES
      hemelb::vis::controller->streaklines(mSimulationState.TimeStep, mLbm->period, &mGlobLatDat,
                                           mLocalLatDat);
#endif

      if (mSimulationState.DoRendering && !write_snapshot_image)
      {
        hemelb::vis::controller->render(RECV_BUFFER_A, &mGlobLatDat, mNetworkTopology);

        if (hemelb::vis::controller->mouse_x >= 0 && hemelb::vis::controller->mouse_y >= 0
            && steeringCpt->updatedMouseCoords)
        {
          for (int i = 0; i < hemelb::vis::controller->col_pixels_recv[RECV_BUFFER_A]; i++)
          {
            if (hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.isRt
                && int (hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.i)
                    == hemelb::vis::controller->mouse_x
                && int (hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i].i.j)
                    == hemelb::vis::controller->mouse_y)
            {
              double mouse_pressure, mouse_stress;
              mLbm->CalculateMouseFlowField(
                                            &hemelb::vis::controller->col_pixel_recv[RECV_BUFFER_A][i],
                                            mouse_pressure, mouse_stress,
                                            hemelb::vis::controller->density_threshold_min,
                                            hemelb::vis::controller->density_threshold_minmax_inv,
                                            hemelb::vis::controller->stress_threshold_max_inv);

              hemelb::vis::controller->setMouseParams(mouse_pressure, mouse_stress);

              break;
            }
          }
          steeringCpt->updatedMouseCoords = 0;
        }

        if (mNetworkTopology->IsCurrentProcTheIOProc())
        {
          imageSendCpt->isFrameReady = true;
        }
      }

      if (render_for_network_stream && mNetworkTopology->IsCurrentProcTheIOProc())
      {
        imageSendCpt->DoWork();
      }

      if (write_snapshot_image)
      {
        hemelb::vis::controller->render(RECV_BUFFER_B, &mGlobLatDat, mNetworkTopology);
        mImagesWritten++;

        if (mNetworkTopology->IsCurrentProcTheIOProc())
        {
          char image_filename[255];
          snprintf(image_filename, 255, "%08i.dat", mSimulationState.TimeStep);

          hemelb::vis::controller->writeImage(RECV_BUFFER_B, image_directory
              + std::string(image_filename), hemelb::vis::ColourPalette::pickColour);
        }
      }

      double lPreSnapshotTime = MPI_Wtime();
      mImagingTime += (lPreSnapshotTime - lPreImageTime);

      if (mSimulationState.TimeStep % snapshots_period == 0)
      {
        char snapshot_filename[255];
        snprintf(snapshot_filename, 255, "snapshot_%06i.dat", mSimulationState.TimeStep);

        mSnapshotsWritten++;
        mLbm->WriteConfigParallel(stability, snapshot_directory + std::string(snapshot_filename),
                                  mGlobLatDat, *mLocalLatDat);
      }

      mSnapshotTime += (MPI_Wtime() - lPreSnapshotTime);

      if (stability == hemelb::lb::StableAndConverged)
      {
        is_finished = true;
        break;
      }
      if (mSimulationState.IsTerminating)
      {
        is_finished = true;
        break;
      }
      if (mLbm->period > 400000)
      {
        is_unstable = true;
        break;
      }
    }

    if (restart)
    {
      hemelb::util::DeleteDirContents(snapshot_directory);
      hemelb::util::DeleteDirContents(image_directory);

      mLbm->Restart(mLocalLatDat);
      mStabilityTester->Reset();
      steeringCpt->Reset();

#ifndef NO_STREAKLINES
      hemelb::vis::controller->restart();
#endif
      if (mNetworkTopology->IsCurrentProcTheIOProc())
      {
        printf("restarting: period: %i\n", mLbm->period);
        fflush(0x0);
      }
      snapshots_period
          = (lSnapshotsPerCycle == 0)
            ? 1e9
            : hemelb::util::NumericalFunctions::max<unsigned int>(1U, mLbm->period
                / lSnapshotsPerCycle);

      images_period = (lImagesPerCycle == 0)
        ? 1e9
        : hemelb::util::NumericalFunctions::max<unsigned int>(1U, mLbm->period / lImagesPerCycle);

      mSimulationState.CycleId = 0;
      continue;
    }
    mLbm->CalculateFlowFieldValues();

    if (mNetworkTopology->IsCurrentProcTheIOProc())
    {
      fprintf(mTimingsFile, "cycle id: %li\n", mSimulationState.CycleId);
      printf("cycle id: %li\n", mSimulationState.CycleId);

      fflush(NULL);
    }
  }

  mSimulationState.CycleId
      = hemelb::util::NumericalFunctions::min<unsigned int>(mSimulationState.CycleId,
                                                            lSimulationConfig->NumCycles);

  PostSimulation(total_time_steps, hemelb::util::myClock() - simulation_time, is_unstable,
                 iStartTime);
}

/**
 * Called on error to abort the simulation and pull-down the MPI environment.
 */
void SimulationMaster::Abort()
{
  int err = MPI_Abort(MPI_COMM_WORLD, 1);

  // This gives us something to work from when we have an error - we get the rank
  // that calls abort, and we get a stack-trace from the exception having been thrown.
  fprintf(stderr, "Aborted by rank %d\n", mNetworkTopology->GetLocalRank());
  throw "SimulationMaster::Abort() called.";
}

/**
 * Steps that are taken when the simulation is complete.
 *
 * This function writes several bits of timing data to
 * the timing file.
 */
void SimulationMaster::PostSimulation(int iTotalTimeSteps,
                                      double iSimulationTime,
                                      bool iIsUnstable,
                                      double iStartTime)
{
  if (mNetworkTopology->IsCurrentProcTheIOProc())
  {
    fprintf(mTimingsFile, "\n");
    fprintf(mTimingsFile, "threads: %i, machines checked: %i\n\n",
            mNetworkTopology->GetProcessorCount(), mNetworkTopology->GetMachineCount());
    fprintf(mTimingsFile, "topology depths checked: %i\n\n", mNetworkTopology->GetDepths());
    fprintf(mTimingsFile, "fluid sites: %i\n\n", mLbm->total_fluid_sites);
    fprintf(mTimingsFile, "cycles and total time steps: %li, %i \n\n", mSimulationState.CycleId,
            iTotalTimeSteps);
    fprintf(mTimingsFile, "time steps per second: %.3f\n\n", iTotalTimeSteps / iSimulationTime);
  }

  if (iIsUnstable)
  {
    if (mNetworkTopology->IsCurrentProcTheIOProc())
    {
      fprintf(mTimingsFile, "Attention: simulation unstable with %i timesteps/cycle\n",
              mLbm->period);
      fprintf(mTimingsFile, "Simulation is terminated\n");
    }
  }
  else
  {
    if (mNetworkTopology->IsCurrentProcTheIOProc())
    {

      fprintf(mTimingsFile, "time steps per cycle: %i\n", mLbm->period);
      fprintf(mTimingsFile, "pressure min, max (mmHg): %e, %e\n", mLbm->GetMinPhysicalPressure(),
              mLbm->GetMaxPhysicalPressure());
      fprintf(mTimingsFile, "velocity min, max (m/s) : %e, %e\n", mLbm->GetMinPhysicalVelocity(),
              mLbm->GetMaxPhysicalVelocity());
      fprintf(mTimingsFile, "stress   min, max (Pa)  : %e, %e\n", mLbm->GetMinPhysicalStress(),
              mLbm->GetMaxPhysicalStress());
      fprintf(mTimingsFile, "\n");

      for (int n = 0; n < mLbm->inlets; n++)
      {
        fprintf(mTimingsFile, "inlet id: %i, average / peak velocity (m/s): %e / %e\n", n,
                mLbm->GetAverageInletVelocity(n), mLbm->GetPeakInletVelocity(n));
      }
      fprintf(mTimingsFile, "\n");

      fprintf(mTimingsFile, "\n");
      fprintf(mTimingsFile, "domain decomposition time (s):             %.3f\n", mDomainDecompTime);
      fprintf(mTimingsFile, "pre-processing buffer management time (s): %.3f\n", mNetInitialiseTime);
      fprintf(mTimingsFile, "input configuration reading time (s):      %.3f\n", mFileReadTime);

      fprintf(mTimingsFile, "total time (s):                            %.3f\n\n",
               (hemelb::util::myClock() - iStartTime));

      fprintf(mTimingsFile, "Sub-domains info:\n\n");

      for (unsigned int n = 0; n < mNetworkTopology->GetProcessorCount(); n++)
      {
        fprintf(mTimingsFile, "rank: %i, fluid sites: %i\n", n,
                mNetworkTopology->FluidSitesOnEachProcessor[n]);
      }
    }
  }

  PrintTimingData();
}

/**
 * Outputs a breakdown of the simulation time spent on different activities.
 */
void SimulationMaster::PrintTimingData()
{
  double lTimings[5] = { mLbTime, mMPISendTime, mMPIWaitTime, mImagingTime, mSnapshotTime };
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

  if (mNetworkTopology->IsCurrentProcTheIOProc())
  {
    for (int ii = 0; ii < 3; ii++)
      lTimings[ii] = 0.0;
  }

  MPI_Reduce(lTimings, lMaxes, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(lTimings, lMeans, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // Change the values for LBM and MPI on process 0 so they don't interfere with the min
  // operation (previously values were 0.0 so they won't affect max / mean
  // calc).
  if (mNetworkTopology->IsCurrentProcTheIOProc())
  {
    for (int ii = 0; ii < 3; ii++)
      lTimings[ii] = std::numeric_limits<double>::max();
  }
  MPI_Reduce(lTimings, lMins, 5, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if (mNetworkTopology->IsCurrentProcTheIOProc())
  {
    if (mNetworkTopology->GetProcessorCount() > 1)
    {
      for (int ii = 0; ii < 3; ii++)
        lMeans[ii] /= (double) (mNetworkTopology->GetProcessorCount() - 1);
      for (int ii = 3; ii < 5; ii++)
        lMeans[ii] /= (double) mNetworkTopology->GetProcessorCount();
    }

    fprintf(mTimingsFile, "\n\nPer-proc timing data (secs per [cycle,image,snapshot]): \n\n");
    fprintf(mTimingsFile, "\t\tMin \tMean \tMax\n");
    for (int ii = 0; ii < 5; ii++)
    {
      fprintf(mTimingsFile, "%s\t\t%.3g\t%.3g\t%.3g\n", lNames[ii].c_str(), lMins[ii], lMeans[ii],
              lMaxes[ii]);
    }
  }
}
