#include "SimulationMaster.h"
#include "configuration/SimConfig.h"

#include "io/XdrFileWriter.h"
#include "util/utilityFunctions.h"
#include "geometry/LatticeData.h"
#include "debug/Debugger.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "lb/HFunction.h"

#include "topology/NetworkTopology.h"

#include <map>
#include <limits>
#include <stdlib.h>

/**
 * Constructor for the SimulationMaster class
 *
 * Initialises member variables including the network topology
 * object.
 */
SimulationMaster::SimulationMaster(hemelb::configuration::CommandLine & options):
  timings()
{
  if (options.HasProblems()) {
    Abort();
  }

  timings[hemelb::reporting::Timers::total].Start();

  hemelb::debug::Debugger::Init(options.Arguments()[0]);

  mLatDat = NULL;

  mLbm = NULL;
  steeringCpt = NULL;
  mVisControl = NULL;
  mSimulationState = NULL;

  mImagesWritten = 0;
  mSnapshotsWritten = 0;

  snapshotsPerCycle = options.NumberOfSnapshotsPerCycle();
  imagesPerCycle = options.NumberOfImagesPerCycle();
  steeringSessionId = options.GetSteeringSessionId();;
  fileManager=new hemelb::reporting::FileManager(options,IsCurrentProcTheIOProc(),GetProcessorCount());
  if (fileManager->HasProblems())
  {
    Abort();
  }
  simConfig = hemelb::configuration::SimConfig::Load(fileManager->GetInputFile().c_str());
  //simConfig->Save(fileManager->GetInputFile()+"foo");
  fileManager->SaveConfiguration(simConfig);

  Initialise();

}

/**
 * Destructor for the SimulationMaster class.
 *
 * Deallocates dynamically allocated memory to contained classes.
 */
SimulationMaster::~SimulationMaster()
{


  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    delete imageSendCpt;
  }

  if (mLatDat != NULL)
  {
    delete mLatDat;
  }

  if (mLbm != NULL)
  {
    delete mLbm;
  }

  if (mInletValues != NULL)
  {
    delete mInletValues;
  }

  if (mOutletValues != NULL)
  {
    delete mOutletValues;
  }

  if (network != NULL)
  {
    delete network;
  }
  if (steeringCpt != NULL)
  {
    delete steeringCpt;
  }

  if (mVisControl != NULL)
  {
    delete mVisControl;
  }

  if (mStabilityTester != NULL)
  {
    delete mStabilityTester;
  }

  if (mEntropyTester != NULL)
  {
    delete mEntropyTester;
  }

  if (mSimulationState != NULL)
  {
    delete mSimulationState;
  }

  if (mUnits != NULL)
  {
    delete mUnits;
  }
  delete simConfig;
  delete fileManager;
}

/**
 * Returns true if the current processor is the dedicated I/O
 * processor.
 */
bool SimulationMaster::IsCurrentProcTheIOProc()
{
  return hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc();
}

/**
 * Returns the number of processors involved in the simulation.
 */
int SimulationMaster::GetProcessorCount()
{
  return hemelb::topology::NetworkTopology::Instance()->GetProcessorCount();
}

/**
 * Initialises various elements of the simulation if necessary - steering,
 * domain decomposition, LBM and visualisation.
 */
void SimulationMaster::Initialise()
{

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Beginning Initialisation.");

  mSimulationState = new hemelb::lb::SimulationState(simConfig->StepsPerCycle,
                                                     simConfig->NumCycles);

  hemelb::site_t mins[3], maxes[3];
  // TODO The way we initialise LbmParameters is not great.
  hemelb::lb::LbmParameters params(1000, 0.1);
  hemelb::site_t totalFluidSites;

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LatticeData.");
  mLatDat
  = new hemelb::geometry::LatticeData(hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                                      &totalFluidSites,
                                      mins,
                                      maxes,
                                      hemelb::topology::NetworkTopology::Instance()->FluidSitesOnEachProcessor,
                                      &params,
                                      simConfig,
                                      timings);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LBM.");
  mLbm = new hemelb::lb::LBM(simConfig, &mNet, mLatDat, mSimulationState, timings[hemelb::reporting::Timers::lb]);
  mLbm->SetSiteMinima(mins);
  mLbm->SetSiteMaxima(maxes);

  mLbm->GetLbmParams()->StressType = params.StressType;
  mLbm->SetTotalFluidSiteCount(totalFluidSites);

  // Initialise and begin the steering.
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    network = new hemelb::steering::Network(steeringSessionId);
  }
  else
  {
    network = NULL;
  }

  mStabilityTester = new hemelb::lb::StabilityTester(mLatDat, &mNet, mSimulationState);

  if (hemelb::log::Logger::ShouldDisplay<hemelb::log::Debug>())
  {
    int typesTested[1] = { 0 };
    mEntropyTester
    = new hemelb::lb::EntropyTester(typesTested, 1, mLatDat, &mNet, mSimulationState);
  }
  else
  {
    mEntropyTester = NULL;
  }

  timings[hemelb::reporting::Timers::netInitialise].Start();
  hemelb::site_t* lReceiveTranslator = mNet.Initialise(mLatDat);
  timings[hemelb::reporting::Timers::netInitialise].Stop();

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising visualisation controller.");
  mVisControl = new hemelb::vis::Control(mLbm->GetLbmParams()->StressType,
                                         &mNet,
                                         mSimulationState,
                                         mLatDat, timings[hemelb::reporting::Timers::visualisation]);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    imageSendCpt = new hemelb::steering::ImageSendComponent(mLbm,
                                                            mSimulationState,
                                                            mVisControl,
                                                            mLbm->GetLbmParams(),
                                                            network);
  }

  mUnits = new hemelb::util::UnitConverter(mLbm->GetLbmParams(),
                                           mSimulationState,
                                           mLatDat->GetVoxelSize());

  mInletValues
  = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::LatticeData::INLET_TYPE,
                                               mLatDat,
                                               simConfig->Inlets,
                                               mSimulationState,
                                               mUnits);

  mOutletValues
  = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::LatticeData::OUTLET_TYPE,
                                               mLatDat,
                                               simConfig->Outlets,
                                               mSimulationState,
                                               mUnits);

  mLbm->Initialise(lReceiveTranslator, mVisControl, mInletValues, mOutletValues, mUnits);

  steeringCpt = new hemelb::steering::SteeringComponent(network,
                                                        mVisControl,
                                                        mLbm,
                                                        &mNet,
                                                        mSimulationState,
                                                        simConfig,
                                                        mUnits);

  // Read in the visualisation parameters.
  mLbm->ReadVisParameters();
}

unsigned int SimulationMaster::OutputPeriod(unsigned int frequency){
  if (frequency==0){
    return 1000000000;
  }
  unsigned long roundedPeriod=mSimulationState->GetTimeStepsPerCycle()/frequency;
  return hemelb::util::NumericalFunctions::max(1U,(unsigned int) roundedPeriod);
}

void SimulationMaster::HandleActors()
{
  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
  != actors.end(); ++it)
  {
    (*it)->RequestComms();
  }

  mNet.Receive();
  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
  != actors.end(); ++it)
  {
    (*it)->PreSend();
  }
  timings[hemelb::reporting::Timers::mpiSend].Start();
  mNet.Send();
  timings[hemelb::reporting::Timers::mpiSend].Stop();

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
  != actors.end(); ++it)
  {
    (*it)->PreReceive();
  }

  timings[hemelb::reporting::Timers::mpiWait].Start();
  mNet.Wait();
  timings[hemelb::reporting::Timers::mpiWait].Stop();

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
  != actors.end(); ++it)
  {
    (*it)->PostReceive();
  }

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
  != actors.end(); ++it)
  {
    (*it)->EndIteration();
  }
}


void SimulationMaster::ResetUnstableSimulation(){
  fileManager->EmptyOutputDirectories();

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
  != actors.end(); ++it)
  {
    (*it)->Reset();
  }

#ifndef NO_STREAKLINES
  mVisControl->Reset();
#endif

  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("restarting: period: %i\n",
                                                                      mSimulationState->GetTimeStepsPerCycle());


  mSimulationState->Reset();
}

void SimulationMaster::WriteLocalImages(){
  for (std::multimap<unsigned long, unsigned long>::const_iterator it =
      snapshotsCompleted.find(mSimulationState->GetTimeStepsPassed()); it
      != snapshotsCompleted.end() && it->first == mSimulationState->GetTimeStepsPassed(); ++it)
  {
    mImagesWritten++;

    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {

      hemelb::io::XdrFileWriter * writer = fileManager->XdrImageWriter(
          1 + ( (it->second - 1) % mSimulationState->GetTimeStepsPerCycle()));

      const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result =
          mVisControl->GetResult(it->second);

      mVisControl ->WriteImage(writer,
                               *result,
                               mVisControl->mDomainStats,
                               mVisControl->mVisSettings);
    }
  }

  snapshotsCompleted.erase(mSimulationState->GetTimeStepsPassed());
}

void SimulationMaster::GenerateNetworkImages(){
  for (std::multimap<unsigned long, unsigned long>::const_iterator it =
      networkImagesCompleted.find(mSimulationState->GetTimeStepsPassed()); it
      != networkImagesCompleted.end() && it->first == mSimulationState->GetTimeStepsPassed(); ++it)
  {
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {

      const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result =
          mVisControl->GetResult(it->second);

      if (steeringCpt->updatedMouseCoords)
      {
        float density, stress;

        if (mVisControl->MouseIsOverPixel(result, &density, &stress))
        {
          double mouse_pressure, mouse_stress;
          mLbm->CalculateMouseFlowField(density,
                                        stress,
                                        mouse_pressure,
                                        mouse_stress,
                                        mVisControl->mDomainStats.density_threshold_min,
                                        mVisControl->mDomainStats.density_threshold_minmax_inv,
                                        mVisControl->mDomainStats.stress_threshold_max_inv);

          mVisControl->SetMouseParams(mouse_pressure, mouse_stress);
        }
        steeringCpt->updatedMouseCoords = false;
      }

      imageSendCpt->DoWork(result);

    }
  }

  networkImagesCompleted.erase(mSimulationState->GetTimeStepsPassed());
}

/**
 * Begin the simulation.
 */
 void SimulationMaster::RunSimulation()
{
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Beginning to run simulation.");

  timings[hemelb::reporting::Timers::simulation].Start();
  bool is_unstable = false;
  int total_time_steps = 0;

  unsigned int snapshots_period = OutputPeriod(snapshotsPerCycle);
  unsigned int images_period = OutputPeriod(imagesPerCycle);

  bool is_finished = false;
  hemelb::lb::Stability stability = hemelb::lb::Stable;


  actors.push_back(mLbm);
  actors.push_back(mInletValues);
  actors.push_back(mOutletValues);
  actors.push_back(steeringCpt);
  actors.push_back(mStabilityTester);
  if (mEntropyTester != NULL)
  {
    actors.push_back(mEntropyTester);
  }
  actors.push_back(mVisControl);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    actors.push_back(network);
  }

  for (; mSimulationState->GetTimeStepsPassed() <= mSimulationState->GetTotalTimeSteps()
  && !is_finished; mSimulationState->Increment())
  {
    ++total_time_steps;

    bool write_snapshot_image = ( (mSimulationState->GetTimeStep() % images_period) == 0)
                          ? true
                              : false;

    // Make sure we're rendering if we're writing this iteration.
    if (write_snapshot_image)
    {
      snapshotsCompleted.insert(std::pair<unsigned long, unsigned long>(mVisControl->Start(),
                                                                        mSimulationState->GetTimeStepsPassed()));
    }

    if (mSimulationState->GetDoRendering())
    {
      networkImagesCompleted.insert(std::pair<unsigned long, unsigned long>(mVisControl->Start(),
                                                                            mSimulationState->GetTimeStepsPassed()));
      mSimulationState->SetDoRendering(false);
    }

    /* In the following two if blocks we do the core magic to ensure we only Render
     when (1) we are not sending a frame or (2) we need to output to disk */

    /* for debugging purposes we want to ensure we capture the variables in a single
     instant of time since variables might be altered by the thread half way through?
     This is to be done. */

    bool render_for_network_stream = false;
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {
      render_for_network_stream = imageSendCpt->ShouldRenderNewNetworkImage();
      steeringCpt->readyForNextImage = render_for_network_stream;
    }

    if (mSimulationState->GetTimeStep() % 100 == 0)
    {
      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %i render_network_stream %i write_snapshot_image %i rendering %i",
                                                                          mSimulationState->GetTimeStep(),
                                                                          render_for_network_stream,
                                                                          write_snapshot_image,
                                                                          mSimulationState->GetDoRendering());
    }

    HandleActors();

    stability = hemelb::lb::Stable;

    if (mSimulationState->GetStability() == hemelb::lb::Unstable)
    {
      ResetUnstableSimulation();
      snapshots_period = OutputPeriod(snapshotsPerCycle);
      images_period = OutputPeriod(imagesPerCycle);
      continue;
    }
    mLbm->UpdateInletVelocities(mSimulationState->GetTimeStep());

#ifndef NO_STREAKLINES
    mVisControl->ProgressStreaklines(mSimulationState->GetTimeStep(),
                                     mSimulationState->GetTimeStepsPerCycle());
#endif

    if (snapshotsCompleted.count(mSimulationState->GetTimeStepsPassed()) > 0)
    {
      WriteLocalImages();

    }

    if (networkImagesCompleted.count(mSimulationState->GetTimeStepsPassed()) > 0)
    {
      GenerateNetworkImages();
    }

    timings[hemelb::reporting::Timers::snapshot].Start();

    if (mSimulationState->GetTimeStep() % snapshots_period == 0)
    {
      mSnapshotsWritten++;
      mLbm->WriteConfigParallel(stability, fileManager->SnapshotPath(mSimulationState->GetTimeStep()));
    }

    timings[hemelb::reporting::Timers::snapshot].Stop();

    if (stability == hemelb::lb::StableAndConverged)
    {
      is_finished = true;
      break;
    }
    if (mSimulationState->GetIsTerminating())
    {
      is_finished = true;
      break;
    }
    if (mSimulationState->GetTimeStepsPerCycle() > 400000)
    {
      is_unstable = true;
      break;
    }

    if (mSimulationState->GetTimeStep() == mSimulationState->GetTimeStepsPerCycle()
        && IsCurrentProcTheIOProc())
    {
      fileManager->Report()->Cycle(mSimulationState->GetCycleId());

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("cycle id: %li",
                                                                          mSimulationState->GetCycleId());

      fflush( NULL);
    }
  }
  timings[hemelb::reporting::Timers::simulation].Stop();
  PostSimulation(total_time_steps, is_unstable);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Finish running simulation.");
}

 /**
  * Called on error to abort the simulation and pull-down the MPI environment.
  */
 void SimulationMaster::Abort()
 {
   int err = MPI_Abort(MPI_COMM_WORLD, 1);

   // This gives us something to work from when we have an error - we get the rank
   // that calls abort, and we get a stack-trace from the exception having been thrown.
   hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Aborting");
   exit(1);
 }

 /**
  * Steps that are taken when the simulation is complete.
  *
  * This function writes several bits of timing data to
  * the timing file.
  */
 void SimulationMaster::PostSimulation(int iTotalTimeSteps, bool iIsUnstable)
 {
   timings[hemelb::reporting::Timers::total].Stop();
   timings.Reduce();
   if (IsCurrentProcTheIOProc())
   {
     fileManager->Report()->Phase1(mLbm->TotalFluidSiteCount(),
                               iTotalTimeSteps,
                               mSimulationState->GetCycleId(),
                               iIsUnstable, mSimulationState->GetTimeStepsPerCycle(), mImagesWritten, mSnapshotsWritten,
                               timings);
   }

 }

