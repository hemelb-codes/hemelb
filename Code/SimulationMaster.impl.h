// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_SIMULATIONMASTER_IMPL_H
#define HEMELB_SIMULATIONMASTER_IMPL_H

#include "SimulationMaster.h"

#include <map>
#include <limits>
#include <cstdlib>
#include <boost/uuid/uuid_io.hpp>

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
#include "net/IOCommunicator.h"
#include "colloids/BodyForces.h"
#include "colloids/BoundaryConditions.h"

#ifdef HEMELB_BUILD_RBC
#  include "redblood/CellController.h"
#  include "redblood/FaderCell.h"
#  include "redblood/MeshIO.h"
#  include "redblood/RBCConfig.h"
#endif

namespace hemelb
{
  /**
   * Constructor for the SimulationMaster class
   *
   * Initialises member variables including the network topology
   * object.
   */
  template<class TRAITS>
  SimulationMaster<TRAITS>::SimulationMaster(hemelb::configuration::CommandLine & options,
                                             const hemelb::net::IOCommunicator& ioComm) :
      ioComms(ioComm.Duplicate()), timings(ioComms), build_info(), cellController(nullptr),
          communicationNet(ioComms)
  {
    timings[hemelb::reporting::Timers::total].Start();

    imagesPerSimulation = options.NumberOfImages();
    steeringSessionId = options.GetSteeringSessionId();

    fileManager = std::make_shared<hemelb::io::PathManager>(options,
                                                            IsCurrentProcTheIOProc(),
                                                            GetProcessorCount());
    simConfig.reset(hemelb::configuration::SimConfig::New(fileManager->GetInputFile()));
    unitConverter = &simConfig->GetUnitConverter();
    monitoringConfig = simConfig->GetMonitoringConfiguration();

    fileManager->SaveConfiguration(simConfig.get());
    Initialise();
    if (IsCurrentProcTheIOProc())
    {
      reporter = std::make_shared<hemelb::reporting::Reporter>(fileManager->GetReportPath(),
                                                               fileManager->GetInputFile());
      reporter->AddReportable(&build_info);
      if (monitoringConfig->doIncompressibilityCheck)
      {
        reporter->AddReportable(incompressibilityChecker.get());
      }
      reporter->AddReportable(&timings);
      reporter->AddReportable(latticeData.get());
      reporter->AddReportable(simulationState.get());
    }
  }

  /**
   * Destructor for the SimulationMaster class.
   *
   * Deallocates dynamically allocated memory to contained classes.
   */
  template<class TRAITS> SimulationMaster<TRAITS>::~SimulationMaster()
  {
  }

  /**
   * Returns true if the current processor is the dedicated I/O
   * processor.
   */
  template<class TRAITS>
  bool SimulationMaster<TRAITS>::IsCurrentProcTheIOProc()
  {
    return ioComms.OnIORank();
  }

  /**
   * Returns the number of processors involved in the simulation.
   */
  template<class TRAITS>
  int SimulationMaster<TRAITS>::GetProcessorCount()
  {
    return ioComms.Size();
  }

  /**
   * Initialises various elements of the simulation if necessary - steering,
   * domain decomposition, LBM and visualisation.
   */
  template<class TRAITS>
  void SimulationMaster<TRAITS>::Initialise()
  {

    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Beginning Initialisation.");

    simulationState = std::make_shared<hemelb::lb::SimulationState>(simConfig->GetTimeStepLength(),
                                                                    simConfig->GetTotalTimeSteps());

    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Initialising LatticeData.");

    timings[hemelb::reporting::Timers::latDatInitialise].Start();
    // Use a reader to read in the file.
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Loading file and decomposing geometry.");

    hemelb::geometry::GeometryReader reader(hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                                            latticeType::GetLatticeInfo(),
                                            timings,
                                            ioComms);
    hemelb::geometry::Geometry readGeometryData =
        reader.LoadAndDecompose(simConfig->GetDataFilePath());

    // Create a new lattice based on that info and return it.
    latticeData = std::make_shared<hemelb::geometry::LatticeData>(latticeType::GetLatticeInfo(),
                                                                  readGeometryData,
                                                                  ioComms);

    timings[hemelb::reporting::Timers::latDatInitialise].Stop();

    neighbouringDataManager = std::make_shared<
        hemelb::geometry::neighbouring::NeighbouringDataManager>(*latticeData,
                                                                 latticeData->GetNeighbouringData(),
                                                                 communicationNet);
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Initialising LBM.");
    latticeBoltzmannModel =
        std::make_shared<hemelb::lb::LBM<Traits>>(simConfig.get(),
                                                  &communicationNet,
                                                  latticeData.get(),
                                                  simulationState.get(),
                                                  timings,
                                                  neighbouringDataManager.get());

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
      hemelb::colloids::BoundaryConditions::InitBoundaryConditions(latticeData.get(), xml);

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Initialising Colloids.");
      colloidController =
          std::make_shared<hemelb::colloids::ColloidController>(*latticeData,
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

    if (simConfig->HasRBCSection())
    {
#ifdef HEMELB_BUILD_RBC
      auto rbcConfig = simConfig->GetRBCConfig();
      hemelb::redblood::CellContainer cells;
      typedef hemelb::redblood::CellController<Traits> Controller;
      auto const controller = std::make_shared<Controller>(*latticeData,
                                                           cells,
                                                           rbcConfig->GetRBCMeshes(),
                                                           timings,
                                                           rbcConfig->GetBoxSize(),
                                                           rbcConfig->GetCell2Cell(),
                                                           rbcConfig->GetCell2Wall(),
                                                           ioComms);
      controller->SetCellInsertion(rbcConfig->GetInserter());
      controller->SetOutlets(*rbcConfig->GetRBCOutlets());
      cellController = std::static_pointer_cast<hemelb::net::IteratedAction>(controller);

      auto output_callback =
          [this](const hemelb::redblood::CellContainer & cells)
          {
	    auto rbcConfig = simConfig->GetRBCConfig();
            auto timestep = simulationState->Get0IndexedTimeStep();
            if ((timestep % rbcConfig->GetRBCOutputPeriod()) == 0)
            {
              log::Logger::Log<log::Info, log::OnePerCore>("printstep %d, num cells %d", timestep, cells.size());

              // Create output directory for current writing step. Requires syncing to
              // ensure no process goes ahead before directory is created.
              std::string rbcOutputDir;
              try
              {
                rbcOutputDir = fileManager->GetRBCOutputPathWithSubdir(std::to_string(timestep));
              }
              catch(Exception& e)
              {
                std::stringstream message;
                message << e.what() << std::endl
                        << "Error " << errno << ": " << std::strerror(errno);
                log::Logger::Log<log::Critical, log::OnePerCore>(message.str());
                ioComms.Abort(-1);
                exit(-1);
              }
              ioComms.Barrier();

              for (auto cell : cells)
              {
                std::stringstream filename;
                filename << rbcOutputDir << cell->GetTag() << "_t_" << timestep << ".vtp";

                std::shared_ptr<redblood::CellBase> cell_base;
                auto fader_cell_cast = std::dynamic_pointer_cast<redblood::FaderCell>(cell);
                if(fader_cell_cast)
                {
                  cell_base = fader_cell_cast->GetWrapeeCell();
                }
                else
                {
                  cell_base = cell;
                }
                auto cell_cast = std::dynamic_pointer_cast<redblood::Cell>(cell_base);
                assert(cell_cast);
		auto meshio = redblood::VTKMeshIO{};
		meshio.writeFile(filename.str(), *cell_cast, *unitConverter);
              }
            }
          };
      controller->AddCellChangeListener(output_callback);
#else
      throw hemelb::Exception() << "Trying to create red blood cell controller with HEMELB_BUILD_RBC=OFF";
#endif
    }

    // Initialise and begin the steering.
    if (ioComms.OnIORank())
    {
      network = std::make_shared<hemelb::steering::Network>(steeringSessionId, timings);
    }

    stabilityTester =
        std::make_shared<hemelb::lb::StabilityTester<latticeType>>(latticeData.get(),
                                                                   &communicationNet,
                                                                   simulationState.get(),
                                                                   timings,
                                                                   monitoringConfig);
    if (monitoringConfig->doIncompressibilityCheck)
    {
      incompressibilityChecker =
          std::make_shared<
              hemelb::lb::IncompressibilityChecker<hemelb::net::PhasedBroadcastRegular<>>>(latticeData.get(),
                                                                                           &communicationNet,
                                                                                           simulationState.get(),
                                                                                           latticeBoltzmannModel->GetPropertyCache(),
                                                                                           timings);
    }

    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Initialising visualisation controller.");
    visualisationControl =
        std::make_shared<hemelb::vis::Control>(latticeBoltzmannModel->GetLbmParams()->StressType,
                                               &communicationNet,
                                               simulationState.get(),
                                               latticeBoltzmannModel->GetPropertyCache(),
                                               latticeData.get(),
                                               timings[hemelb::reporting::Timers::visualisation]);

    if (ioComms.OnIORank())
    {
      imageSendCpt =
          std::make_shared<hemelb::steering::ImageSendComponent>(simulationState.get(),
                                                                 visualisationControl.get(),
                                                                 latticeBoltzmannModel->GetLbmParams(),
                                                                 network.get(),
                                                                 latticeBoltzmannModel->InletCount());

    }

    inletValues = std::make_shared<hemelb::lb::iolets::BoundaryValues>(hemelb::geometry::INLET_TYPE,
                                                                       latticeData.get(),
                                                                       simConfig->GetInlets(),
                                                                       simulationState.get(),
                                                                       ioComms,
                                                                       *unitConverter);

    outletValues =
        std::make_shared<hemelb::lb::iolets::BoundaryValues>(hemelb::geometry::OUTLET_TYPE,
                                                             latticeData.get(),
                                                             simConfig->GetOutlets(),
                                                             simulationState.get(),
                                                             ioComms,
                                                             *unitConverter);

    latticeBoltzmannModel->Initialise(visualisationControl.get(),
                                      inletValues.get(),
                                      outletValues.get(),
                                      unitConverter);
    latticeBoltzmannModel->SetInitialConditions(ioComms);
    neighbouringDataManager->ShareNeeds();
    neighbouringDataManager->TransferNonFieldDependentInformation();

    steeringCpt = std::make_shared<hemelb::steering::SteeringComponent>(network.get(),
                                                                        visualisationControl.get(),
                                                                        imageSendCpt.get(),
                                                                        &communicationNet,
                                                                        simulationState.get(),
                                                                        simConfig.get(),
                                                                        unitConverter);

    // Read in the visualisation parameters.
    latticeBoltzmannModel->ReadVisParameters();

    propertyDataSource =
        std::make_shared<hemelb::extraction::LbDataSourceIterator>(latticeBoltzmannModel->GetPropertyCache(),
                                                                   *latticeData,
                                                                   ioComms.Rank(),
                                                                   *unitConverter);

    if (simConfig->PropertyOutputCount() > 0)
    {

      for (unsigned outputNumber = 0; outputNumber < simConfig->PropertyOutputCount();
          ++outputNumber)
      {
        simConfig->GetPropertyOutput(outputNumber)->filename = fileManager->GetDataExtractionPath()
            + simConfig->GetPropertyOutput(outputNumber)->filename;
      }

      propertyExtractor =
          std::make_shared<hemelb::extraction::PropertyActor>(*simulationState,
                                                              simConfig->GetPropertyOutputs(),
                                                              *propertyDataSource,
                                                              timings,
                                                              ioComms);
    }

    imagesPeriod = OutputPeriod(imagesPerSimulation);

    stepManager =
        std::make_shared<hemelb::net::phased::StepManager>(2,
                                                           &timings,
                                                           hemelb::net::separate_communications);
    netConcern = std::make_shared<hemelb::net::phased::NetConcern>(communicationNet);
    stepManager->RegisterIteratedActorSteps(*neighbouringDataManager.get(), 0);
    if (colloidController)
    {
      stepManager->RegisterIteratedActorSteps(*colloidController, 1);
    }
    if (cellController)
    {
      stepManager->RegisterIteratedActorSteps(*cellController, 1);
    }
    stepManager->RegisterIteratedActorSteps(*latticeBoltzmannModel, 1);

    stepManager->RegisterIteratedActorSteps(*inletValues, 1);
    stepManager->RegisterIteratedActorSteps(*outletValues, 1);
    stepManager->RegisterIteratedActorSteps(*steeringCpt, 1);
    stepManager->RegisterIteratedActorSteps(*stabilityTester, 1);
    if (entropyTester)
    {
      stepManager->RegisterIteratedActorSteps(*entropyTester, 1);
    }

    if (monitoringConfig->doIncompressibilityCheck)
    {
      stepManager->RegisterIteratedActorSteps(*incompressibilityChecker, 1);
    }
    stepManager->RegisterIteratedActorSteps(*visualisationControl, 1);
    if (propertyExtractor)
    {
      stepManager->RegisterIteratedActorSteps(*propertyExtractor, 1);
    }

    if (ioComms.OnIORank())
    {
      stepManager->RegisterIteratedActorSteps(*network, 1);
    }
    stepManager->RegisterCommsForAllPhases(*netConcern);
  }

  template<class TRAITS>
  unsigned int SimulationMaster<TRAITS>::OutputPeriod(unsigned int frequency)
  {
    if (frequency == 0)
    {
      return 1000000000;
    }
    unsigned long roundedPeriod = simulationState->GetTotalTimeSteps() / frequency;
    return hemelb::util::NumericalFunctions::max(1U, (unsigned int) roundedPeriod);
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::HandleActors()
  {
    stepManager->CallActions();
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::OnUnstableSimulation()
  {
    LogStabilityReport();
    hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Aborting: time step length: %f\n",
                                                                           simulationState->GetTimeStepLength());
    Finalise();
    Abort();
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::WriteLocalImages()
  {
    /**
     * this map iteration iterates over all those image generation requests completing this step.
     * The map key (it->first) is the completion time step number.
     * The map value (it->second) is the initiation time step number.
     */
    for (MapType::const_iterator it = writtenImagesCompleted.find(simulationState->GetTimeStep());
        it != writtenImagesCompleted.end() && it->first == simulationState->GetTimeStep(); ++it)
    {

      if (ioComms.OnIORank())
      {
        reporter->Image();
        std::unique_ptr<hemelb::io::writers::Writer> writer(fileManager->XdrImageWriter(1
            + ( (it->second - 1) % simulationState->GetTimeStep())));

        const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result =
            visualisationControl->GetResult(it->second);

        visualisationControl->WriteImage(writer.get(),
                                         *result,
                                         visualisationControl->domainStats,
                                         visualisationControl->visSettings);

      }
    }

    writtenImagesCompleted.erase(simulationState->GetTimeStep());
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::GenerateNetworkImages()
  {
    for (std::multimap<unsigned long, unsigned long>::const_iterator it =
        networkImagesCompleted.find(simulationState->GetTimeStep());
        it != networkImagesCompleted.end() && it->first == simulationState->GetTimeStep(); ++it)
    {
      if (ioComms.OnIORank())
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
  template<class TRAITS>
  void SimulationMaster<TRAITS>::RunSimulation()
  {
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Beginning to run simulation.");
    timings[hemelb::reporting::Timers::simulation].Start();

    while (simulationState->GetTimeStep() <= simulationState->GetTotalTimeSteps())
    {
      DoTimeStep();
      if (simulationState->IsTerminating())
      {
        break;
      }
    }

    timings[hemelb::reporting::Timers::simulation].Stop();
    Finalise();
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::Finalise()
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

  template<class TRAITS>
  void SimulationMaster<TRAITS>::DoTimeStep()
  {
    log::Logger::Log<log::Debug, log::OnePerCore>("Current LB time: %e",
                                                  simulationState->GetTime());
    bool writeImage = ( (simulationState->GetTimeStep() % imagesPeriod) == 0) ?
      true :
      false;

    // Make sure we're rendering if we're writing this iteration.
    if (writeImage)
    {
      /***
       * writtenImagesCompleted and networkImagesCompleted are multimaps.
       * The keys are the iterations on which production of an image will complete, and should be written or sent over the network.
       * The values are the iterations on which the image creation began.
       */
      writtenImagesCompleted.insert(std::pair<unsigned long, unsigned long>(visualisationControl->Start(),
                                                                            simulationState->GetTimeStep()));
    }

    if (simulationState->IsRendering())
    {
      // Here, Start() actually triggers the render.
      networkImagesCompleted.insert(std::pair<unsigned long, unsigned long>(visualisationControl->Start(),
                                                                            simulationState->GetTimeStep()));
      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("%d images currently being composited for the steering client",
                                                                           networkImagesCompleted.size());
      simulationState->SetIsRendering(false);
    }

    /* In the following two if blocks we do the core magic to ensure we only Render
     when (1) we are not sending a frame or (2) we need to output to disk */

    /* TODO for debugging purposes we want to ensure we capture the variables in a single
     instant of time since variables might be altered by the thread half way through?
     This is to be done. */

    bool renderForNetworkStream = false;
    if (ioComms.OnIORank() && !steeringCpt->readyForNextImage)
    {
      renderForNetworkStream = imageSendCpt->ShouldRenderNewNetworkImage();
      steeringCpt->readyForNextImage = renderForNetworkStream;
    }

    if (simulationState->GetTimeStep() % 100 == 0)
    {
      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %i render_network_stream %i write_image_to_disk %i rendering %i",
                                                                          simulationState->GetTimeStep(),
                                                                          renderForNetworkStream,
                                                                          writeImage,
                                                                          simulationState->IsRendering());
      LogStabilityReport();
    }

    RecalculatePropertyRequirements();

    HandleActors();

    if (simulationState->GetStability() == hemelb::lb::Unstable)
    {
      OnUnstableSimulation();
    }

    // If the user requested to terminate converged steady flow simulations, mark
    // simulation to be finished at the end of the current timestep.
    if ( (simulationState->GetStability() == hemelb::lb::StableAndConverged)
        && monitoringConfig->convergenceTerminate)
    {
      LogStabilityReport();
      simulationState->SetIsTerminating(true);
    }

    if ( (simulationState->GetTimeStep() % 500 == 0) && colloidController)
      colloidController->OutputInformation(simulationState->GetTimeStep());

#ifndef NO_STREAKLINES
    visualisationControl->ProgressStreaklines(simulationState->GetTimeStep(),
                                              simulationState->GetTotalTimeSteps());
#endif

    if (writtenImagesCompleted.count(simulationState->GetTimeStep()) > 0)
    {
      WriteLocalImages();

    }

    if (networkImagesCompleted.count(simulationState->GetTimeStep()) > 0)
    {
      GenerateNetworkImages();
    }

    if (simulationState->GetTimeStep() % FORCE_FLUSH_PERIOD == 0 && IsCurrentProcTheIOProc())
    {
      fflush(nullptr);
    }
    simulationState->Increment();
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::RecalculatePropertyRequirements()
  {
    // Get the property cache & reset its list of properties to get.
    hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

    propertyCache.ResetRequirements();

    // Check whether we're rendering images on this iteration.
    if (visualisationControl->IsRendering())
    {
      propertyCache.densityCache.SetRefreshFlag();
      propertyCache.velocityCache.SetRefreshFlag();

      if (simConfig->GetStressType() == hemelb::lb::ShearStress)
      {
        propertyCache.wallShearStressMagnitudeCache.SetRefreshFlag();
      }
      else if (simConfig->GetStressType() == hemelb::lb::VonMises)
      {
        propertyCache.vonMisesStressCache.SetRefreshFlag();
      }
    }

    if (monitoringConfig->doIncompressibilityCheck)
    {
      propertyCache.densityCache.SetRefreshFlag();
      propertyCache.velocityCache.SetRefreshFlag();
    }

    // If extracting property results, check what's required by them.
    if (propertyExtractor)
    {
      propertyExtractor->SetRequiredProperties(propertyCache);
    }

    // If using streaklines, the velocity will be needed.
#ifndef NO_STREAKLINES
    propertyCache.velocityCache.SetRefreshFlag();
#endif
  }

  /**
   * Called on error to abort the simulation and pull-down the MPI environment.
   */
  template<class TRAITS>
  void SimulationMaster<TRAITS>::Abort()
  {
    // This gives us something to work from when we have an error - we get the rank
    // that calls abort, and we get a stack-trace from the exception having been thrown.
    hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::OnePerCore>("Aborting");
    hemelb::net::MpiEnvironment::Abort(1);

    exit(1);
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::LogStabilityReport()
  {
    if (monitoringConfig->doIncompressibilityCheck
        && incompressibilityChecker->AreDensitiesAvailable())
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

  template<class TRAITS>
  const hemelb::util::UnitConverter& SimulationMaster<TRAITS>::GetUnitConverter() const
  {
    return *unitConverter;
  }
}

#endif
