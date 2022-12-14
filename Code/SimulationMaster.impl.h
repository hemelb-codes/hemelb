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
#include "configuration/SimBuilder.h"
#include "extraction/PropertyActor.h"
#include "extraction/LbDataSourceIterator.h"
#include "io/writers/XdrFileWriter.h"
#include "util/utilityFunctions.h"
#include "geometry/Domain.h"
#include "log/Logger.h"
#include "lb/HFunction.h"
#include "io/xml.h"
#include "net/BuildInfo.h"
#include "net/IOCommunicator.h"
#include "colloids/BodyForces.h"
#include "colloids/BoundaryConditions.h"

#ifdef HEMELB_BUILD_COLLOIDS
#  include "colloids/ColloidController.h"
#endif

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
    SimulationMaster<TRAITS>::SimulationMaster(configuration::CommandLine & options,
                                               const net::IOCommunicator& ioComm) :
            ioComms(ioComm.Duplicate()), timings(ioComms), build_info(),
            communicationNet(ioComms)
    {
        // Start the main timer!
        timings[reporting::Timers::total].Start();

        fileManager = std::make_shared<io::PathManager>(options,
                                                                IsCurrentProcTheIOProc(),
                                                                GetProcessorCount());
        log::Logger::Log<log::Info, log::Singleton>("Reading configuration from %s", fileManager->GetInputFile().c_str());
        // Convert XML to configuration
        simConfig = configuration::SimConfig::New(fileManager->GetInputFile());
        // Use it to initialise self
        auto builder = configuration::SimBuilder(*simConfig);
        log::Logger::Log<log::Info, log::Singleton>("Beginning Initialisation.");
        builder(*this);
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
   * Initialises various elements of the simulation if necessary:
   * domain decomposition, LBM and visualisation.
   */
  template<class TRAITS>
  void SimulationMaster<TRAITS>::Initialise()
  {

  }

  template<class TRAITS>
  unsigned int SimulationMaster<TRAITS>::OutputPeriod(unsigned int frequency)
  {
    if (frequency == 0)
    {
      return 1000000000;
    }
    unsigned long roundedPeriod = simulationState->GetTotalTimeSteps() / frequency;
    return util::NumericalFunctions::max(1U, (unsigned int) roundedPeriod);
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
    log::Logger::Log<log::Warning, log::Singleton>("Aborting: time step length: %f\n",
                                                                           simulationState->GetTimeStepLength());
    Finalise();
    Abort();
  }

  /**
   * Begin the simulation.
   */
  template<class TRAITS>
  void SimulationMaster<TRAITS>::RunSimulation()
  {
    log::Logger::Log<log::Info, log::Singleton>("Beginning to run simulation.");
    timings[reporting::Timers::simulation].Start();

    while (simulationState->GetTimeStep() <= simulationState->GetTotalTimeSteps())
    {
      DoTimeStep();
      if (simulationState->IsTerminating())
      {
        break;
      }
    }

    timings[reporting::Timers::simulation].Stop();
    Finalise();
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::Finalise()
  {
    timings[reporting::Timers::total].Stop();
    timings.Reduce();
    if (IsCurrentProcTheIOProc())
    {
      reporter->FillDictionary();
      reporter->Write();
    }
    // DTMP: Logging output on communication as debug output for now.
    log::Logger::Log<log::Debug, log::OnePerCore>("sync points: %lld, bytes sent: %lld",
                                                                          communicationNet.SyncPointsCounted,
                                                                          communicationNet.BytesSent);

    log::Logger::Log<log::Info, log::Singleton>("Finish running simulation.");
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::DoTimeStep()
  {
    log::Logger::Log<log::Debug, log::OnePerCore>("Current LB time: %e",
                                                  simulationState->GetTime());

    /* In the following two if blocks we do the core magic to ensure we only Render
     when (1) we are not sending a frame or (2) we need to output to disk */

    /* TODO for debugging purposes we want to ensure we capture the variables in a single
     instant of time since variables might be altered by the thread half way through?
     This is to be done. */

    if (simulationState->GetTimeStep() % 100 == 0)
    {
      log::Logger::Log<log::Info, log::Singleton>("time step %i",
                                                                          simulationState->GetTimeStep());
      LogStabilityReport();
    }

    RecalculatePropertyRequirements();

    HandleActors();

    if (simulationState->GetStability() == lb::Unstable)
    {
      OnUnstableSimulation();
    }

    // If the user requested to terminate converged steady flow simulations, mark
    // simulation to be finished at the end of the current timestep.
    if ( (simulationState->GetStability() == lb::StableAndConverged)
        && stabilityTester->ShouldTerminateWhenConverged())
    {
      LogStabilityReport();
      simulationState->SetIsTerminating(true);
    }

#ifdef HEMELB_BUILD_COLLOIDS
    if (simulationState->GetTimeStep() % 500 == 0) {
        if (auto c = std::dynamic_pointer_cast<colloids::ColloidController>(colloidController)) {
            c->OutputInformation(simulationState->GetTimeStep());
        }
    }
#endif
    if (simulationState->GetTimeStep() % FORCE_FLUSH_PERIOD == 0 && IsCurrentProcTheIOProc())
    {
      fflush(nullptr);
    }

    fieldData->SwapOldAndNew();
    simulationState->Increment();
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::RecalculatePropertyRequirements()
  {
    // Get the property cache & reset its list of properties to get.
    lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

    propertyCache.ResetRequirements();

    if (incompressibilityChecker)
    {
      propertyCache.densityCache.SetRefreshFlag();
      propertyCache.velocityCache.SetRefreshFlag();
    }

    // If extracting property results, check what's required by them.
    if (propertyExtractor)
    {
      propertyExtractor->SetRequiredProperties(propertyCache);
    }
  }

  /**
   * Called on error to abort the simulation and pull-down the MPI environment.
   */
  template<class TRAITS>
  void SimulationMaster<TRAITS>::Abort()
  {
    // This gives us something to work from when we have an error - we get the rank
    // that calls abort, and we get a stack-trace from the exception having been thrown.
    log::Logger::Log<log::Critical, log::OnePerCore>("Aborting");
    net::MpiEnvironment::Abort(1);

    exit(1);
  }

  template<class TRAITS>
  void SimulationMaster<TRAITS>::LogStabilityReport()
  {
    if (incompressibilityChecker && incompressibilityChecker->AreDensitiesAvailable())
    {
      log::Logger::Log<log::Info, log::Singleton>("time step %i, tau %.6f, max_relative_press_diff %.3f, Ma %.3f, max_vel_phys %e",
                                                                          simulationState->GetTimeStep(),
                                                                          latticeBoltzmannModel->GetLbmParams()->GetTau(),
                                                                          incompressibilityChecker->GetMaxRelativeDensityDifference(),
                                                                          incompressibilityChecker->GetGlobalLargestVelocityMagnitude()
                                                                              / Cs,
                                                                          unitConverter->ConvertVelocityToPhysicalUnits(incompressibilityChecker->GetGlobalLargestVelocityMagnitude()));
    }

    if (simulationState->GetStability() == lb::StableAndConverged)
    {
      log::Logger::Log<log::Info, log::Singleton>("time step %i, steady flow simulation converged.",
                                                                          simulationState->GetTimeStep());
    }
  }

  template<class TRAITS>
  const util::UnitConverter& SimulationMaster<TRAITS>::GetUnitConverter() const
  {
    return *unitConverter;
  }
}

#endif
