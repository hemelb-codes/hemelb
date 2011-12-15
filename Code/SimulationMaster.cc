#include "SimulationMaster.h"
#include "configuration/SimConfig.h"

#include "io/writers/xdr/XdrFileWriter.h"
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
SimulationMaster::SimulationMaster(hemelb::configuration::CommandLine & options) :
  timings()
{
  if (options.HasProblems())
  {
    Abort();
  }

  timings[hemelb::reporting::Timers::total].Start();

  hemelb::debug::Debugger::Init(options.Arguments()[0]);

  latticeData = NULL;

  latticeBoltzmannModel = NULL;
  steeringCpt = NULL;
  visualisationControl = NULL;
  simulationState = NULL;

  snapshotsPerCycle = options.NumberOfSnapshotsPerCycle();
  imagesPerCycle = options.NumberOfImagesPerCycle();
  steeringSessionId = options.GetSteeringSessionId();
  ;
  fileManager = new hemelb::io::PathManager(options, IsCurrentProcTheIOProc(), GetProcessorCount());
  if (fileManager->HasProblems())
  {
    Abort();
  }
  simConfig = hemelb::configuration::SimConfig::Load(fileManager->GetInputFile().c_str());
  fileManager->SaveConfiguration(simConfig);
  Initialise();
  if (IsCurrentProcTheIOProc())
  {
    reporter = new hemelb::reporting::Reporter(fileManager->GetReportPath(),
                                               fileManager->GetInputFile(),
                                               latticeData->GetFluidSiteCountsOnEachProc(),
                                               latticeData->GetTotalFluidSites(),
                                               timings,
                                               *simulationState,
                                               *incompressibilityChecker);
  }
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
  delete latticeData;
  delete latticeBoltzmannModel;
  delete inletValues;
  delete outletValues;
  delete network;
  delete steeringCpt;
  delete visualisationControl;
  delete stabilityTester;
  delete entropyTester;
  delete simulationState;
  delete unitConvertor;
  delete incompressibilityChecker;

  delete simConfig;
  delete fileManager;
  if (IsCurrentProcTheIOProc())
  {
    delete reporter;
  }
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

  simulationState = new hemelb::lb::SimulationState(simConfig->StepsPerCycle, simConfig->NumCycles);

  // TODO The way we initialise LbmParameters is not great.
  hemelb::lb::LbmParameters params(1000, 0.1);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LatticeData.");

  timings[hemelb::reporting::Timers::netInitialise].Start();
  latticeData
      = hemelb::geometry::LatticeData::Load(hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                                            simConfig->DataFilePath,
                                            &params,
                                            timings);
  timings[hemelb::reporting::Timers::netInitialise].Stop();

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LBM.");
  latticeBoltzmannModel = new hemelb::lb::LBM(simConfig,
                                              &communicationNet,
                                              latticeData,
                                              simulationState,
                                              timings[hemelb::reporting::Timers::lb]);

  latticeBoltzmannModel->GetLbmParams()->StressType = params.StressType;

  // Initialise and begin the steering.
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    network = new hemelb::steering::Network(steeringSessionId);
  }
  else
  {
    network = NULL;
  }

  stabilityTester = new hemelb::lb::StabilityTester(latticeData,
                                                    &communicationNet,
                                                    simulationState,
                                                    timings);
  entropyTester = NULL;
  incompressibilityChecker = new hemelb::lb::IncompressibilityChecker<>(latticeData,
                                                                        &communicationNet,
                                                                        simulationState,
                                                                        timings);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising visualisation controller.");
  visualisationControl
      = new hemelb::vis::Control(latticeBoltzmannModel->GetLbmParams()->StressType,
                                 &communicationNet,
                                 simulationState,
                                 latticeData,
                                 timings[hemelb::reporting::Timers::visualisation]);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    imageSendCpt = new hemelb::steering::ImageSendComponent(latticeBoltzmannModel,
                                                            simulationState,
                                                            visualisationControl,
                                                            latticeBoltzmannModel->GetLbmParams(),
                                                            network);
  }

  unitConvertor = new hemelb::util::UnitConverter(latticeBoltzmannModel->GetLbmParams(),
                                                  simulationState,
                                                  latticeData->GetVoxelSize());

  inletValues
      = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::LatticeData::INLET_TYPE,
                                                   latticeData,
                                                   simConfig->Inlets,
                                                   simulationState,
                                                   unitConvertor);

  outletValues
      = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::LatticeData::OUTLET_TYPE,
                                                   latticeData,
                                                   simConfig->Outlets,
                                                   simulationState,
                                                   unitConvertor);

  latticeBoltzmannModel->Initialise(visualisationControl, inletValues, outletValues, unitConvertor);

  steeringCpt = new hemelb::steering::SteeringComponent(network,
                                                        visualisationControl,
                                                        latticeBoltzmannModel,
                                                        &communicationNet,
                                                        simulationState,
                                                        simConfig,
                                                        unitConvertor);

  // Read in the visualisation parameters.
  latticeBoltzmannModel->ReadVisParameters();
}

unsigned int SimulationMaster::OutputPeriod(unsigned int frequency)
{
  if (frequency == 0)
  {
    return 1000000000;
  }
  unsigned long roundedPeriod = simulationState->GetTimeStepsPerCycle() / frequency;
  return hemelb::util::NumericalFunctions::max(1U, (unsigned int) roundedPeriod);
}

void SimulationMaster::HandleActors()
{
  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it != actors.end(); ++it)
  {
    (*it)->RequestComms();
  }

  communicationNet.Receive();
  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it != actors.end(); ++it)
  {
    (*it)->PreSend();
  }
  timings[hemelb::reporting::Timers::mpiSend].Start();
  communicationNet.Send();
  timings[hemelb::reporting::Timers::mpiSend].Stop();

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it != actors.end(); ++it)
  {
    (*it)->PreReceive();
  }

  timings[hemelb::reporting::Timers::mpiWait].Start();
  communicationNet.Wait();
  timings[hemelb::reporting::Timers::mpiWait].Stop();

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it != actors.end(); ++it)
  {
    (*it)->PostReceive();
  }

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it != actors.end(); ++it)
  {
    (*it)->EndIteration();
  }
}

void SimulationMaster::ResetUnstableSimulation()
{
  fileManager->EmptyOutputDirectories();

  for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it != actors.end(); ++it)
  {
    (*it)->Reset();
  }

#ifndef NO_STREAKLINES
  visualisationControl->Reset();
#endif

  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("restarting: period: %i\n",
                                                                      simulationState->GetTimeStepsPerCycle());

  simulationState->Reset();
}

void SimulationMaster::WriteLocalImages()
{
  for (std::multimap<unsigned long, unsigned long>::const_iterator it =
      snapshotsCompleted.find(simulationState->GetTimeStepsPassed()); it
      != snapshotsCompleted.end() && it->first == simulationState->GetTimeStepsPassed(); ++it)
  {

    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {
      reporter->Image();
      hemelb::io::writers::xdr::XdrFileWriter * writer = fileManager->XdrImageWriter(1
          + ( (it->second - 1) % simulationState->GetTimeStepsPerCycle()));

      const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result =
          visualisationControl->GetResult(it->second);

      visualisationControl->WriteImage(writer,
                                       *result,
                                       visualisationControl->mDomainStats,
                                       visualisationControl->mVisSettings);

      delete writer;
    }
  }

  snapshotsCompleted.erase(simulationState->GetTimeStepsPassed());
}

void SimulationMaster::GenerateNetworkImages()
{
  for (std::multimap<unsigned long, unsigned long>::const_iterator it =
      networkImagesCompleted.find(simulationState->GetTimeStepsPassed()); it
      != networkImagesCompleted.end() && it->first == simulationState->GetTimeStepsPassed(); ++it)
  {
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {

      const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result =
          visualisationControl->GetResult(it->second);

      if (steeringCpt->updatedMouseCoords)
      {
        float density, stress;

        if (visualisationControl->MouseIsOverPixel(result, &density, &stress))
        {
          double mousePressure = 0.0, mouseStress = 0.0;
          latticeBoltzmannModel->CalculateMouseFlowField(density,
                                                         stress,
                                                         mousePressure,
                                                         mouseStress,
                                                         visualisationControl->mDomainStats.density_threshold_min,
                                                         visualisationControl->mDomainStats.density_threshold_minmax_inv,
                                                         visualisationControl->mDomainStats.stress_threshold_max_inv);

          visualisationControl->SetMouseParams(mousePressure, mouseStress);
        }
        steeringCpt->updatedMouseCoords = false;
      }

      imageSendCpt->DoWork(result);

    }
  }

  networkImagesCompleted.erase(simulationState->GetTimeStepsPassed());
}

/**
 * Begin the simulation.
 */
void SimulationMaster::RunSimulation()
{
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Beginning to run simulation.");

  timings[hemelb::reporting::Timers::simulation].Start();
  unsigned int snapshotsPeriod = OutputPeriod(snapshotsPerCycle);
  unsigned int imagesPeriod = OutputPeriod(imagesPerCycle);

  bool isFinished = false;
  hemelb::lb::Stability stability = hemelb::lb::Stable;

  actors.push_back(latticeBoltzmannModel);
  actors.push_back(inletValues);
  actors.push_back(outletValues);
  actors.push_back(steeringCpt);
  actors.push_back(stabilityTester);
  if (entropyTester != NULL)
  {
    actors.push_back(entropyTester);
  }
  actors.push_back(incompressibilityChecker);
  actors.push_back(visualisationControl);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    actors.push_back(network);
  }

  for (; simulationState->GetTimeStepsPassed() <= simulationState->GetTotalTimeSteps()
      && !isFinished; simulationState->Increment())
  {

    bool writeSnapshotImage = ( (simulationState->GetTimeStep() % imagesPeriod) == 0)
      ? true
      : false;

    // Make sure we're rendering if we're writing this iteration.
    if (writeSnapshotImage)
    {
      snapshotsCompleted.insert(std::pair<unsigned long, unsigned long>(visualisationControl->Start(),
                                                                        simulationState->GetTimeStepsPassed()));
    }

    if (simulationState->GetDoRendering())
    {
      networkImagesCompleted.insert(std::pair<unsigned long, unsigned long>(visualisationControl->Start(),
                                                                            simulationState->GetTimeStepsPassed()));
      simulationState->SetDoRendering(false);
    }

    /* In the following two if blocks we do the core magic to ensure we only Render
     when (1) we are not sending a frame or (2) we need to output to disk */

    /* for debugging purposes we want to ensure we capture the variables in a single
     instant of time since variables might be altered by the thread half way through?
     This is to be done. */

    bool renderForNetworkStream = false;
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {
      renderForNetworkStream = imageSendCpt->ShouldRenderNewNetworkImage();
      steeringCpt->readyForNextImage = renderForNetworkStream;
    }

    if (simulationState->GetTimeStep() % 100 == 0)
    {
      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %i render_network_stream %i write_snapshot_image %i rendering %i",
                                                                          simulationState->GetTimeStep(),
                                                                          renderForNetworkStream,
                                                                          writeSnapshotImage,
                                                                          simulationState->GetDoRendering());

    }

    HandleActors();

    stability = hemelb::lb::Stable;

    if (simulationState->GetStability() == hemelb::lb::Unstable)
    {
      ResetUnstableSimulation();
      snapshotsPeriod = OutputPeriod(snapshotsPerCycle);
      imagesPeriod = OutputPeriod(imagesPerCycle);
      continue;
    }
    latticeBoltzmannModel->UpdateInletVelocities(simulationState->GetTimeStep());

#ifndef NO_STREAKLINES
    visualisationControl->ProgressStreaklines(simulationState->GetTimeStep(),
                                              simulationState->GetTimeStepsPerCycle());
#endif

    if (snapshotsCompleted.count(simulationState->GetTimeStepsPassed()) > 0)
    {
      WriteLocalImages();

    }

    if (networkImagesCompleted.count(simulationState->GetTimeStepsPassed()) > 0)
    {
      GenerateNetworkImages();
    }

    timings[hemelb::reporting::Timers::snapshot].Start();

    if (simulationState->GetTimeStep() % snapshotsPeriod == 0)
    {
      if (IsCurrentProcTheIOProc())
      {
        reporter->Snapshot();
      }
      latticeBoltzmannModel->WriteConfigParallel(stability,
                                                 fileManager->SnapshotPath(simulationState->GetTimeStep()));
    }

    timings[hemelb::reporting::Timers::snapshot].Stop();

    if (stability == hemelb::lb::StableAndConverged)
    {
      isFinished = true;
      break;
    }
    if (simulationState->GetIsTerminating())
    {
      isFinished = true;
      break;
    }
    if (simulationState->GetTimeStepsPerCycle() > 400000)
    {
      if (IsCurrentProcTheIOProc())
      {
        reporter->Stability(false);
      }
      break;
    }

    if (simulationState->GetTimeStep() == simulationState->GetTimeStepsPerCycle()
        && IsCurrentProcTheIOProc())
    {

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("cycle id: %li",
                                                                          simulationState->GetCycleId());

      fflush(NULL);
    }
  }
  timings[hemelb::reporting::Timers::simulation].Stop();
  timings[hemelb::reporting::Timers::total].Stop();
  timings.Reduce();
  if (IsCurrentProcTheIOProc())
  {
    reporter->FillDictionary();
    reporter->Write();
  }
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Finish running simulation.");
}

/**
 * Called on error to abort the simulation and pull-down the MPI environment.
 */
void SimulationMaster::Abort()
{
  MPI_Abort(MPI_COMM_WORLD, 1);

  // This gives us something to work from when we have an error - we get the rank
  // that calls abort, and we get a stack-trace from the exception having been thrown.
  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Aborting");
  exit(1);
}

