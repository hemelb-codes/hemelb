#include "SimulationMaster.h"
#include "configuration/SimConfig.h"
#include "extraction/PropertyActor.h"
#include "extraction/LbDataSourceIterator.h"
#include "io/writers/xdr/XdrFileWriter.h"
#include "util/utilityFunctions.h"
#include "geometry/GeometryReader.h"
#include "geometry/LatticeData.h"
#include "debug/Debugger.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "lb/HFunction.h"
#include "colloids/ColloidController.h"

#include "topology/NetworkTopology.h"

#include <map>
#include <limits>
#include <cstdlib>

/**
 * Constructor for the SimulationMaster class
 *
 * Initialises member variables including the network topology
 * object.
 */
SimulationMaster::SimulationMaster(hemelb::configuration::CommandLine & options) :
    timings(), build_info()
{
  if (options.HasProblems())
  {
    Abort();
  }

  timings[hemelb::reporting::Timers::total].Start();

  hemelb::debug::Debugger::Init(options.Arguments()[0]);

  latticeData = NULL;

  colloidController = NULL;
  latticeBoltzmannModel = NULL;
  steeringCpt = NULL;
  propertyDataSource = NULL;
  visualisationControl = NULL;
  propertyExtractor = NULL;
  simulationState = NULL;
  snapshotsPerSimulation = options.NumberOfSnapshots();
  imagesPerSimulation = options.NumberOfImages();
  steeringSessionId = options.GetSteeringSessionId();

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
    reporter = new hemelb::reporting::Reporter(fileManager->GetReportPath(), fileManager->GetInputFile());
    reporter->AddReportable(&build_info);
    reporter->AddReportable(incompressibilityChecker);
    reporter->AddReportable(&timings);
    reporter->AddReportable(latticeData);
    reporter->AddReportable(simulationState);
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
  delete colloidController;
  delete latticeBoltzmannModel;
  delete inletValues;
  delete outletValues;
  delete network;
  delete steeringCpt;
  delete visualisationControl;
  delete propertyExtractor;
  delete propertyDataSource;
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

  simulationState = new hemelb::lb::SimulationState(simConfig->GetTimeStepLength(), simConfig->GetTotalTimeSteps());

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LatticeData.");

  timings[hemelb::reporting::Timers::latDatInitialise].Start();
// Use a reader to read in the file.
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Loading file and decomposing geometry.");

  hemelb::geometry::GeometryReader reader(
                   hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                   latticeType::GetLatticeInfo(),
                   timings);
                 
  hemelb::geometry::Geometry readGeometryData = reader.LoadAndDecompose(simConfig->GetDataFilePath());

  // Create a new lattice based on that info and return it.
  latticeData = new hemelb::geometry::LatticeData(
                   latticeType::GetLatticeInfo(),
                   readGeometryData);
  timings[hemelb::reporting::Timers::latDatInitialise].Stop();
 hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising Colloids.");
  colloidController = new hemelb::colloids::ColloidController(
                   &communicationNet,
                   latticeData,
                   &readGeometryData);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LBM.");
  latticeBoltzmannModel = new hemelb::lb::LBM<latticeType>(simConfig,
                                                           &communicationNet,
                                                           latticeData,
                                                           simulationState,
                                                           timings);

  // Initialise and begin the steering.
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    network = new hemelb::steering::Network(steeringSessionId,timings);
  }
  else
  {
    network = NULL;
  }

  stabilityTester = new hemelb::lb::StabilityTester<latticeType>(latticeData,
                                                                 &communicationNet,
                                                                 simulationState,
                                                                 timings);
  entropyTester = NULL;
  incompressibilityChecker =
      new hemelb::lb::IncompressibilityChecker<hemelb::net::PhasedBroadcastRegular<> >(latticeData,
                                                                                       &communicationNet,
                                                                                       simulationState,
                                                                                       latticeBoltzmannModel->GetPropertyCache(),
                                                                                       timings);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising visualisation controller.");
  visualisationControl = new hemelb::vis::Control(latticeBoltzmannModel->GetLbmParams()->StressType,
                                                  &communicationNet,
                                                  simulationState,
                                                  latticeBoltzmannModel->GetPropertyCache(),
                                                  latticeData,
                                                  timings[hemelb::reporting::Timers::visualisation]);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    imageSendCpt = new hemelb::steering::ImageSendComponent(simulationState,
                                                            visualisationControl,
                                                            latticeBoltzmannModel->GetLbmParams(),
                                                            network,
                                                            latticeBoltzmannModel->InletCount());

  } else {
    imageSendCpt=NULL;
  }

  unitConvertor = new hemelb::util::UnitConverter(latticeBoltzmannModel->GetLbmParams(),
                                                  simulationState,
                                                  latticeData->GetVoxelSize());

  inletValues = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::INLET_TYPE,
                                                           latticeData,
                                                           simConfig->GetInlets(),
                                                           simulationState,
                                                           unitConvertor);

  outletValues = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::OUTLET_TYPE,
                                                            latticeData,
                                                            simConfig->GetOutlets(),
                                                            simulationState,
                                                            unitConvertor);

  latticeBoltzmannModel->Initialise(visualisationControl, inletValues, outletValues, unitConvertor);

  steeringCpt = new hemelb::steering::SteeringComponent(network,
                                                        visualisationControl,
                                                        imageSendCpt,
                                                        &communicationNet,
                                                        simulationState,
                                                        simConfig,
                                                        unitConvertor);

  // Read in the visualisation parameters.
  latticeBoltzmannModel->ReadVisParameters();

  propertyDataSource = new hemelb::extraction::LbDataSourceIterator(latticeBoltzmannModel->GetPropertyCache(),
                                                                    *latticeData,
                                                                    *unitConvertor);

  if (simConfig->PropertyOutputCount() > 0)
  {

    for (unsigned outputNumber = 0; outputNumber < simConfig->PropertyOutputCount(); ++outputNumber)
    {
      simConfig->GetPropertyOutput(outputNumber)->filename = fileManager->GetDataExtractionPath()
          + simConfig->GetPropertyOutput(outputNumber)->filename;
    }

    propertyExtractor = new hemelb::extraction::PropertyActor(*simulationState,
                                                              simConfig->GetPropertyOutputs(),
                                                              *propertyDataSource);
  }

  imagesPeriod = OutputPeriod(imagesPerSimulation);

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
  if (propertyExtractor != NULL)
  {
    actors.push_back(propertyExtractor);
  }

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    actors.push_back(network);
  }
}

unsigned int SimulationMaster::OutputPeriod(unsigned int frequency)
{
  if (frequency == 0)
  {
    return 1000000000;
  }
  unsigned long roundedPeriod = simulationState->GetTotalTimeSteps() / frequency;
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

  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("restarting: time step length: %i\n",
                                                                      simulationState->GetTimeStepLength());

  simulationState->Reset();
}

void SimulationMaster::WriteLocalImages()
{
  /**
   * this map iteration iterates over all those image generation requests completing this step.
   * The map key (it->first) is the completion time step number.
   * The map value (it->second) is the initiation time step number.
   */
  for (MapType::const_iterator it =
      snapshotsCompleted.find(simulationState->GetTimeStep());
      it != snapshotsCompleted.end() && it->first == simulationState->GetTimeStep(); ++it)
  {

    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {
      reporter->Image();
      hemelb::io::writers::xdr::XdrFileWriter * writer = fileManager->XdrImageWriter(1
          + ( (it->second - 1) % simulationState->GetTimeStep()));

      const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result = visualisationControl->GetResult(it->second);

      visualisationControl->WriteImage(writer,
                                       *result,
                                       visualisationControl->domainStats,
                                       visualisationControl->visSettings);

      delete writer;
    }
  }

  snapshotsCompleted.erase(simulationState->GetTimeStep());
}

void SimulationMaster::GenerateNetworkImages()
{
  for (std::multimap<unsigned long, unsigned long>::const_iterator it =
      networkImagesCompleted.find(simulationState->GetTimeStep());
      it != networkImagesCompleted.end() && it->first == simulationState->GetTimeStep(); ++it)
  {
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {

      const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result = visualisationControl->GetResult(it->second);

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
                                                         visualisationControl->domainStats.density_threshold_min,
                                                         visualisationControl->domainStats.density_threshold_minmax_inv,
                                                         visualisationControl->domainStats.stress_threshold_max_inv);

          visualisationControl->SetMouseParams(mousePressure, mouseStress);
        }
        steeringCpt->updatedMouseCoords = false;
      }

      imageSendCpt->DoWork(result);

    }
  }

  networkImagesCompleted.erase(simulationState->GetTimeStep());
}

/**
 * Begin the simulation.
 */
void SimulationMaster::RunSimulation()
{
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Beginning to run simulation.");
  timings[hemelb::reporting::Timers::simulation].Start();

  while (simulationState->GetTimeStep() <= simulationState->GetTotalTimeSteps())
  {
    DoTimeStep();
    if (simulationState->GetTimeStep() > MAX_TIME_STEPS)
    {
      simulationState->SetStability(hemelb::lb::Unstable);
      break;
    }
    if (simulationState->IsTerminating())
    {
      break;
    }
  }

  timings[hemelb::reporting::Timers::simulation].Stop();
  Finalise();
}

void SimulationMaster::Finalise()
{
  timings[hemelb::reporting::Timers::total].Stop();
  timings.Reduce();
  if (IsCurrentProcTheIOProc())
  {
    reporter->FillDictionary();
    reporter->Write();
  }
  // DTMP: Logging output on communication as debug output for now.
  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("sync points: %lld, bytes sent: %lld",
                                                                        communicationNet.SyncPointsCounted, 
                                                                        communicationNet.BytesSent);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Finish running simulation.");
}

void SimulationMaster::DoTimeStep()
{
  bool writeSnapshotImage = ( (simulationState->GetTimeStep() % imagesPeriod) == 0) ?
    true :
    false;

  // Make sure we're rendering if we're writing this iteration.
  if (writeSnapshotImage)
  {
    /***
     * snapShotsCompleted and networkImagesCompleted are multimaps.
     * The keys are the iterations on which production of an image will complete, and should be written or sent over the network.
     * The values are the iterations on which the image creation began.
     */
    snapshotsCompleted.insert(std::pair<unsigned long, unsigned long>(visualisationControl->Start(),
                                                                      simulationState->GetTimeStep()));
  }

  if (simulationState->IsRendering())
  {
    // Here, Start() actually triggers the render.
    networkImagesCompleted.insert(std::pair<unsigned long, unsigned long>(visualisationControl->Start(),
                                                                          simulationState->GetTimeStep()));
    hemelb::log::Logger::Log<hemelb::log::Debug,
      hemelb::log::Singleton>("%d images currently being composited for the steering client", networkImagesCompleted.size());
    simulationState->SetIsRendering(false);
  }

  /* In the following two if blocks we do the core magic to ensure we only Render
   when (1) we are not sending a frame or (2) we need to output to disk */

  /* TODO for debugging purposes we want to ensure we capture the variables in a single
   instant of time since variables might be altered by the thread half way through?
   This is to be done. */

  bool renderForNetworkStream = false;
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc()&&!steeringCpt->readyForNextImage)
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
                                                                        simulationState->IsRendering());

  }

  RecalculatePropertyRequirements();

  HandleActors();

  if (simulationState->GetStability() == hemelb::lb::Unstable)
  {
    ResetUnstableSimulation();
    imagesPeriod = OutputPeriod(imagesPerSimulation);
    return;
  }

#ifndef NO_STREAKLINES
  visualisationControl->ProgressStreaklines(simulationState->GetTimeStep(), simulationState->GetTotalTimeSteps());
#endif

  if (snapshotsCompleted.count(simulationState->GetTimeStep()) > 0)
  {
    WriteLocalImages();

  }

  if (networkImagesCompleted.count(simulationState->GetTimeStep()) > 0)
  {
    GenerateNetworkImages();
  }

  timings[hemelb::reporting::Timers::snapshot].Start();

  if (ShouldWriteSnapshot())
  {
    if (IsCurrentProcTheIOProc())
    {
      reporter->Snapshot();
    }
    latticeBoltzmannModel->WriteConfigParallel(simulationState->GetStability(),
                                               fileManager->SnapshotPath(simulationState->GetTimeStep()));
  }

  timings[hemelb::reporting::Timers::snapshot].Stop();

  if (simulationState->GetTimeStep() % FORCE_FLUSH_PERIOD == 0 && IsCurrentProcTheIOProc())
  {
    fflush(NULL);
  }
  simulationState->Increment();
}

void SimulationMaster::RecalculatePropertyRequirements()
{
  // Get the property cache & reset its list of properties to get.
  hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

  propertyCache.ResetRequirements();

  // Check whether we're rendering images or snapshotting on this iteration.
  if (visualisationControl->IsRendering() || ShouldWriteSnapshot())
  {
    propertyCache.densityCache.SetRefreshFlag();
    propertyCache.velocityCache.SetRefreshFlag();

    if (simConfig->GetStressType() == hemelb::lb::ShearStress)
    {
      propertyCache.shearStressCache.SetRefreshFlag();
    }
    else if (simConfig->GetStressType() == hemelb::lb::VonMises)
    {
      propertyCache.vonMisesStressCache.SetRefreshFlag();
    }
  }

  if (incompressibilityChecker != NULL)
  {
    propertyCache.densityCache.SetRefreshFlag();
  }

  // If extracting property results, check what's required by them.
  if (propertyExtractor != NULL)
  {
    propertyExtractor->SetRequiredProperties(propertyCache);
  }

  // If using streaklines, the velocity will be needed.
#ifndef NO_STREAKLINES
  propertyCache.velocityCache.SetRefreshFlag();
#endif
}

bool SimulationMaster::ShouldWriteSnapshot()
{
  return simulationState->GetTimeStep() % OutputPeriod(snapshotsPerSimulation) == 0;
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

