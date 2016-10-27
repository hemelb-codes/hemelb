
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "SimulationMaster.h"
#include "configuration/SimConfig.h"
#include "extraction/PropertyActor.h"
#include "extraction/LbDataSourceIterator.h"
#include "io/writers/xdr/XdrFileWriter.h"
#include "util/utilityFunctions.h"
#include "geometry/GeometryReader.h"
#include "geometry/LatticeData.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "lb/HFunction.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "colloids/ColloidController.h"
#include "net/BuildInfo.h"
#include "colloids/BodyForces.h"
#include "colloids/BoundaryConditions.h"
#include "comm/MpiEnvironment.h"

#include <map>
#include <limits>
#include <cstdlib>

/**
 * Constructor for the SimulationMaster class
 *
 * Initialises member variables including the network topology
 * object.
 */
SimulationMaster::SimulationMaster(hemelb::configuration::CommandLine & options, hemelb::comm::Communicator::ConstPtr ioComm) :
  ioComms(ioComm), timings(ioComm), build_info(), communicationNet(ioComm), asyncCommQ(hemelb::comm::Async::New(ioComm))
{
  timings[hemelb::reporting::Timers::total].Start();

  latticeData = NULL;

  colloidController = NULL;
  latticeBoltzmannModel = NULL;
  steeringCpt = NULL;
  propertyDataSource = NULL;
  propertyExtractor = NULL;
  simulationState = NULL;
  stepManager = NULL;
  netConcern = NULL;
  neighbouringDataManager = NULL;
  steeringSessionId = options.GetSteeringSessionId();

  fileManager = new hemelb::io::PathManager(options, IsCurrentProcTheIOProc(), GetProcessorCount());
  simConfig = hemelb::configuration::SimConfig::New(fileManager->GetInputFile());
  unitConverter = &simConfig->GetUnitConverter();
  monitoringConfig = simConfig->GetMonitoringConfiguration();

  fileManager->SaveConfiguration(simConfig);
  Initialise();
  if (IsCurrentProcTheIOProc())
  {
    reporter = new hemelb::reporting::Reporter(fileManager->GetReportPath(),
                                               fileManager->GetInputFile());
    reporter->AddReportable(&build_info);
    if (monitoringConfig->doIncompressibilityCheck)
    {
      reporter->AddReportable(incompressibilityChecker);
    }
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

  delete latticeData;
  delete colloidController;
  delete latticeBoltzmannModel;
  delete inletValues;
  delete outletValues;
  delete network;
  delete steeringCpt;
  delete propertyExtractor;
  delete propertyDataSource;
  delete stabilityTester;
  delete entropyTester;
  delete simulationState;
  delete incompressibilityChecker;
  delete neighbouringDataManager;

  delete simConfig;
  delete fileManager;
  if (IsCurrentProcTheIOProc())
  {
    delete reporter;
  }
  delete stepManager;
  delete netConcern;
}

/**
 * Returns true if the current processor is the dedicated I/O
 * processor.
 */
bool SimulationMaster::IsCurrentProcTheIOProc()
{
  return ioComms->OnIORank();
}

/**
 * Returns the number of processors involved in the simulation.
 */
int SimulationMaster::GetProcessorCount()
{
  return ioComms->Size();
}

/**
 * Initialises various elements of the simulation if necessary - steering,
 * domain decomposition, and LBM.
 */
void SimulationMaster::Initialise()
{

  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Beginning Initialisation.");

  simulationState = new hemelb::lb::SimulationState(simConfig->GetTimeStepLength(),
                                                    simConfig->GetTotalTimeSteps());

  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Initialising LatticeData.");

  timings[hemelb::reporting::Timers::latDatInitialise].Start();
  // Use a reader to read in the file.
  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Loading file and decomposing geometry.");

  hemelb::geometry::GeometryReader reader(hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                                          latticeType::GetLatticeInfo(),
                                          timings, ioComms);
  hemelb::geometry::Geometry readGeometryData =
      reader.LoadAndDecompose(simConfig->GetDataFilePath());

  // Create a new lattice based on that info and return it.
  latticeData = new hemelb::geometry::LatticeData(latticeType::GetLatticeInfo(), readGeometryData, ioComms);

  timings[hemelb::reporting::Timers::latDatInitialise].Stop();

  neighbouringDataManager =
      new hemelb::geometry::neighbouring::NeighbouringDataManager(*latticeData,
                                                                  latticeData->GetNeighbouringData(),
                                                                  asyncCommQ);
  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Initialising LBM.");
  latticeBoltzmannModel = new hemelb::lb::LBM<latticeType>(simConfig,
                                                           asyncCommQ,
                                                           latticeData,
                                                           simulationState,
                                                           timings,
                                                           neighbouringDataManager);

  hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

  if (simConfig->HasColloidSection())
  {
    timings[hemelb::reporting::Timers::colloidInitialisation].Start();
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Loading Colloid config.");
    std::string colloidConfigPath = simConfig->GetColloidConfigPath();
    hemelb::io::xml::Document xml(colloidConfigPath);

    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Creating Body Forces.");
    hemelb::colloids::BodyForces::InitBodyForces(xml);

    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Creating Boundary Conditions.");
    hemelb::colloids::BoundaryConditions::InitBoundaryConditions(latticeData, xml);

    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Initialising Colloids.");
    colloidController =
        new hemelb::colloids::ColloidController(*latticeData,
                                                *simulationState,
                                                readGeometryData,
                                                xml,
                                                propertyCache,
                                                latticeBoltzmannModel->GetLbmParams(),
                                                fileManager->GetColloidPath(),
                                                ioComms,
                                                timings);
  }
  timings[hemelb::reporting::Timers::colloidInitialisation].Stop();

  // Initialise and begin the steering.
  if (ioComms->OnIORank())
  {
    network = new hemelb::steering::Network(steeringSessionId, timings);
  }
  else
  {
    network = NULL;
  }

  stabilityTester = new hemelb::lb::StabilityTester<latticeType>(latticeData,
                                                                 ioComms,
                                                                 simulationState,
                                                                 timings,
                                                   monitoringConfig->doConvergenceCheck,
                                                   monitoringConfig->convergenceRelativeTolerance);
  entropyTester = NULL;

  if (monitoringConfig->doIncompressibilityCheck)
  {
    incompressibilityChecker = new hemelb::lb::IncompressibilityChecker(latticeData,
                                                ioComms,
                                                simulationState,
                                                latticeBoltzmannModel->GetPropertyCache(),
                                                timings);
  }
  else
  {
    incompressibilityChecker = NULL;
  }

  inletValues = new hemelb::lb::iolets::BoundaryValues(hemelb::geometry::INLET_TYPE,
                                                       latticeData,
                                                       simConfig->GetInlets(),
                                                       simulationState,
                                                       ioComms,
                                                       *unitConverter);

  outletValues = new hemelb::lb::iolets::BoundaryValues(hemelb::geometry::OUTLET_TYPE,
                                                        latticeData,
                                                        simConfig->GetOutlets(),
                                                        simulationState,
                                                        ioComms,
                                                        *unitConverter);

  latticeBoltzmannModel->Initialise(inletValues, outletValues, unitConverter);
  neighbouringDataManager->ShareNeeds();
  neighbouringDataManager->TransferNonFieldDependentInformation();

  steeringCpt = new hemelb::steering::SteeringComponent(network,
                                                        &communicationNet,
                                                        simulationState,
                                                        simConfig,
                                                        unitConverter,
                                                        timings);

  propertyDataSource =
      new hemelb::extraction::LbDataSourceIterator(latticeBoltzmannModel->GetPropertyCache(),
                                                   *latticeData,
                                                   ioComms->Rank(),
                                                   *unitConverter);

  if (simConfig->PropertyOutputCount() > 0)
  {

    for (unsigned outputNumber = 0; outputNumber < simConfig->PropertyOutputCount(); ++outputNumber)
    {
      simConfig->GetPropertyOutput(outputNumber)->filename = fileManager->GetDataExtractionPath()
          + simConfig->GetPropertyOutput(outputNumber)->filename;
    }

    propertyExtractor = new hemelb::extraction::PropertyActor(*simulationState,
                                                              simConfig->GetPropertyOutputs(),
                                                              *propertyDataSource,
                                                              timings, ioComms);
  }

  stepManager = new hemelb::net::phased::StepManager(2,
                                                     &timings,
                                                     hemelb::net::separate_communications);
  netConcern = new hemelb::net::phased::NetConcern(communicationNet);
  stepManager->RegisterIteratedActorSteps(*neighbouringDataManager, 0);
  if (colloidController != NULL)
  {
    stepManager->RegisterIteratedActorSteps(*colloidController, 1);
  }
  stepManager->RegisterIteratedActorSteps(*latticeBoltzmannModel, 1);

  stepManager->RegisterIteratedActorSteps(*inletValues, 1);
  stepManager->RegisterIteratedActorSteps(*outletValues, 1);
  stepManager->RegisterIteratedActorSteps(*steeringCpt, 1);
  stepManager->RegisterIteratedActorSteps(*stabilityTester, 1);
  stepManager->RegisterCommsSteps(*stabilityTester, 1);
  if (entropyTester != NULL)
  {
    stepManager->RegisterIteratedActorSteps(*entropyTester, 1);
  }

  if (monitoringConfig->doIncompressibilityCheck)
  {
    stepManager->RegisterIteratedActorSteps(*incompressibilityChecker, 1);
    stepManager->RegisterCommsSteps(*incompressibilityChecker, 1);
  }

  if (propertyExtractor != NULL)
  {
    stepManager->RegisterIteratedActorSteps(*propertyExtractor, 1);
  }

  if (ioComms->OnIORank())
  {
    stepManager->RegisterIteratedActorSteps(*network, 1);
  }
  stepManager->RegisterCommsForAllPhases(*netConcern);
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
  stepManager->CallActions();
}

void SimulationMaster::OnUnstableSimulation()
{
  LogStabilityReport();
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Aborting: time step length: %f\n",
                                                                         simulationState->GetTimeStepLength());
  Finalise();
  Abort();
}


/**
 * Begin the simulation.
 */
void SimulationMaster::RunSimulation()
{
  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Beginning to run simulation.");
  timings[hemelb::reporting::Timers::simulation].Start();

  int LastTS = simulationState->GetTotalTimeSteps();
  while (simulationState->GetTimeStep() <= LastTS)
  {
    // We need to keep the master rank in sync with the workers
    // Here, each rank notifies that it is beginning a time step
    auto syncReq = ioComms->Ibarrier();
    
    DoTimeStep();
    
    syncReq->Wait();
    // Now all ranks have at least started this time step
    
    if (simulationState->IsTerminating())
    {
      LastTS = std::min(simulationState->GetTimeStep(), simulationState->GetTotalTimeSteps());
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

  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Finish running simulation.");
}

void SimulationMaster::DoTimeStep()
{
  // If the simulation is finishing, all running collective actions
  // must wait this time step
  if (simulationState->IsTerminating())
  {
    if (stabilityTester)
      stabilityTester->MustFinishThisTimeStep();
    if (entropyTester)
      entropyTester->MustFinishThisTimeStep();
    if (incompressibilityChecker)
      incompressibilityChecker->MustFinishThisTimeStep();
    if (steeringCpt)
      steeringCpt->MustFinishThisTimeStep();
  }
  
  if (simulationState->GetTimeStep() % 100 == 0)
  {
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %i",
                                                                        simulationState->GetTimeStep());
    LogStabilityReport();
  }

  RecalculatePropertyRequirements();

  HandleActors();

  if (simulationState->GetStability() == hemelb::lb::Unstable)
  {
    OnUnstableSimulation();
  }

  // If the user requested to terminate converged steady flow simulations, mark
  // simulation to be finished ASAP.
  if ( (simulationState->GetStability() == hemelb::lb::StableAndConverged)
      && monitoringConfig->convergenceTerminate)
  {
    LogStabilityReport();
    simulationState->SetIsTerminating(true);
  }

  if ( (simulationState->GetTimeStep() % 500 == 0) && colloidController != NULL)
    colloidController->OutputInformation(simulationState->GetTimeStep());

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

  if (monitoringConfig->doIncompressibilityCheck)
  {
    propertyCache.densityCache.SetRefreshFlag();
    propertyCache.velocityCache.SetRefreshFlag();
  }

  // If extracting property results, check what's required by them.
  if (propertyExtractor != NULL)
  {
    propertyExtractor->SetRequiredProperties(propertyCache);
  }
}

/**
 * Called on error to abort the simulation and pull-down the MPI environment.
 */
void SimulationMaster::Abort()
{
  // This gives us something to work from when we have an error - we get the rank
  // that calls abort, and we get a stack-trace from the exception having been thrown.
  hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::OnePerCore>("Aborting");
  hemelb::comm::MpiEnvironment::Abort(1);

  exit(1);
}

void SimulationMaster::LogStabilityReport()
{
  if (monitoringConfig->doIncompressibilityCheck)
  {
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %i, tau %.6f, max_relative_press_diff %.3f, Ma %.3f, max_vel_phys %e",
                                                                        simulationState->GetTimeStep(),
                                                                        latticeBoltzmannModel->GetLbmParams()->GetTau(),
                                                                        incompressibilityChecker->GetMaxRelativeDensityDifference(),
                                                                        incompressibilityChecker->GetGlobalLargestVelocityMagnitude()
                                                                            / hemelb::Cs,
                                                                        unitConverter->ConvertVelocityToPhysicalUnits(incompressibilityChecker->GetGlobalLargestVelocityMagnitude()));
  }

  if (simulationState->GetStability() == hemelb::lb::StableAndConverged)
  {
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %i, steady flow simulation converged.",
                                                                        simulationState->GetTimeStep());
  }
}

