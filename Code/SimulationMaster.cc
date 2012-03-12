#include "SimulationMaster.h"
#include "configuration/SimConfig.h"
#include "extraction/PropertyActor.h"
#include "extraction/LbDataSourceIterator.h"
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

  simulationState = new hemelb::lb::SimulationState(simConfig->TimeStepLength,simConfig->TotalTimeSteps);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LatticeData.");

  timings[hemelb::reporting::Timers::netInitialise].Start();
  latticeData = hemelb::geometry::LatticeData::Load(hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                                                    latticeType::GetLatticeInfo(),
                                                    simConfig->DataFilePath,
                                                    timings);
  timings[hemelb::reporting::Timers::netInitialise].Stop();

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LBM.");
  latticeBoltzmannModel = new hemelb::lb::LBM<latticeType>(simConfig,
                                                           &communicationNet,
                                                           latticeData,
                                                           simulationState,
                                                           timings[hemelb::reporting::Timers::lb]);

  // Initialise and begin the steering.
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    network = new hemelb::steering::Network(steeringSessionId);
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
      new hemelb::lb::IncompressibilityChecker<hemelb::net::PhasedBroadcastRegular<>, latticeType>(latticeData,
                                                                                                   &communicationNet,
                                                                                                   simulationState,
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
  }

  unitConvertor = new hemelb::util::UnitConverter(latticeBoltzmannModel->GetLbmParams(),
                                                  simulationState,
                                                  latticeData->GetVoxelSize());

  inletValues = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::INLET_TYPE,
                                                           latticeData,
                                                           simConfig->Inlets,
                                                           simulationState,
                                                           unitConvertor);

  outletValues = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::OUTLET_TYPE,
                                                            latticeData,
                                                            simConfig->Outlets,
                                                            simulationState,
                                                            unitConvertor);

  latticeBoltzmannModel->Initialise(visualisationControl, inletValues, outletValues, unitConvertor);

  steeringCpt = new hemelb::steering::SteeringComponent(network,
                                                        visualisationControl,
                                                        &communicationNet,
                                                        simulationState,
                                                        simConfig,
                                                        unitConvertor);

  // Read in the visualisation parameters.
  latticeBoltzmannModel->ReadVisParameters();

  propertyDataSource = new hemelb::extraction::LbDataSourceIterator(latticeBoltzmannModel->GetPropertyCache(),
                                                                    *latticeData,
                                                                    *unitConvertor);

  if (simConfig->propertyOutputs.size() > 0)
  {

    for (unsigned outputNumber = 0; outputNumber < simConfig->propertyOutputs.size(); ++outputNumber)
    {
      simConfig->propertyOutputs[outputNumber]->filename = fileManager->GetDataExtractionPath()
          + simConfig->propertyOutputs[outputNumber]->filename;
    }

    propertyExtractor = new hemelb::extraction::PropertyActor(*simulationState,
                                                              simConfig->propertyOutputs,
                                                              *propertyDataSource);
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
  for (std::multimap<unsigned long, unsigned long>::const_iterator it =
      snapshotsCompleted.find(simulationState->GetTimeStepsPassed());
      it != snapshotsCompleted.end() && it->first == simulationState->GetTimeStepsPassed(); ++it)
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

  snapshotsCompleted.erase(simulationState->GetTimeStepsPassed());
}

void SimulationMaster::GenerateNetworkImages()
{
  for (std::multimap<unsigned long, unsigned long>::const_iterator it =
      networkImagesCompleted.find(simulationState->GetTimeStepsPassed());
      it != networkImagesCompleted.end() && it->first == simulationState->GetTimeStepsPassed(); ++it)
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

  networkImagesCompleted.erase(simulationState->GetTimeStepsPassed());
}

/**
 * Begin the simulation.
 */
void SimulationMaster::RunSimulation()
{
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Beginning to run simulation.");

  timings[hemelb::reporting::Timers::simulation].Start();
  unsigned int imagesPeriod = OutputPeriod(imagesPerSimulation);

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
  if (propertyExtractor != NULL)
  {
    actors.push_back(propertyExtractor);
  }

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    actors.push_back(network);
  }

  for (; simulationState->GetTimeStepsPassed() <= simulationState->GetTotalTimeSteps() && !isFinished;
      simulationState->Increment())
  {

    bool writeSnapshotImage = ( (simulationState->GetTimeStep() % imagesPeriod) == 0) ?
      true :
      false;

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

    RecalculatePropertyRequirements();

    HandleActors();

    stability = hemelb::lb::Stable;

    if (simulationState->GetStability() == hemelb::lb::Unstable)
    {
      ResetUnstableSimulation();
      imagesPeriod = OutputPeriod(imagesPerSimulation);
      continue;
    }

#ifndef NO_STREAKLINES
    visualisationControl->ProgressStreaklines(simulationState->GetTimeStep(), simulationState->GetTotalTimeSteps());
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

    if (IsSnapshotting())
    {
      if (IsCurrentProcTheIOProc())
      {
        reporter->Snapshot();
      }
      latticeBoltzmannModel->WriteConfigParallel(stability, fileManager->SnapshotPath(simulationState->GetTimeStep()));
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
    if (simulationState->GetTotalTimeSteps() > 400000)
    {
      if (IsCurrentProcTheIOProc())
      {
        reporter->Stability(false);
      }
      break;
    }

    if (simulationState->GetTimeStep()%1000==0 && IsCurrentProcTheIOProc())
    {
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

void SimulationMaster::RecalculatePropertyRequirements()
{
  // Get the property cache & reset its list of properties to get.
  hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

  propertyCache.ResetRequirements();

  // Check whether we're rendering images or snapshotting on this iteration.
  if (visualisationControl->IsRendering() || IsSnapshotting())
  {
    propertyCache.densityCache.SetRefreshFlag();
    propertyCache.velocityCache.SetRefreshFlag();

    if (simConfig->StressType == hemelb::lb::ShearStress)
    {
      propertyCache.shearStressCache.SetRefreshFlag();
    }
    else if (simConfig->StressType == hemelb::lb::VonMises)
    {
      propertyCache.vonMisesStressCache.SetRefreshFlag();
    }
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

bool SimulationMaster::IsSnapshotting()
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

