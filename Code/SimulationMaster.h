// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H
#include <memory>

#include "lb/lattices/Lattices.h"
#include "lb/lb.hpp"
#include "lb/StabilityTester.h"
#include "net/net.h"
#include "lb/EntropyTester.h"
#include "lb/iolets/BoundaryValues.h"
#include "util/UnitConverter.h"
#include "configuration/CommandLine.h"
#include "io/PathManager.h"
#include "reporting/Reporter.h"
#include "reporting/Timers.h"
#include "reporting/BuildInfo.h"
#include "lb/IncompressibilityChecker.hpp"
#include "net/phased/StepManager.h"
#include "net/phased/NetConcern.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "Traits.h"

namespace hemelb
{
    namespace configuration { class SimBuilder; }
    namespace extraction { class PropertyActor; }

    template<class TRAITS = Traits<>>
  class SimulationMaster
  {
      friend class configuration::SimBuilder;
    public:
      using Traits = TRAITS;

      SimulationMaster(configuration::CommandLine &options,
                       const net::IOCommunicator& ioComms);
      virtual ~SimulationMaster();

      void Abort();

      bool IsCurrentProcTheIOProc();

      int GetProcessorCount();

      void RunSimulation();
      lb::SimulationState const& GetState() const
      {
        return *simulationState;
      }
      util::UnitConverter const & GetUnitConverter() const;
      void Finalise();
#     ifdef HEMELB_DOING_UNITTESTS
      //! Makes it easy to add cell controller without messy input files
      void RegisterActor(net::phased::Concern &concern, net::phased::StepManager::Phase phase)
      {
        stepManager->RegisterIteratedActorSteps(concern, phase);
      }
      //! Access to lattice data for debugging
      geometry::FieldData& GetFieldData()
      {
        assert(fieldData);
        return *fieldData;
      }
      std::shared_ptr<net::IteratedAction> GetCellController() {
        return cellController;
      }
#     endif
    protected:

      std::shared_ptr<lb::BoundaryValues> inletValues;
      std::shared_ptr<lb::BoundaryValues> outletValues;
      virtual void DoTimeStep();

      /* The next quantities are protected because they are used by MultiscaleSimulationMaster */
      // Set the lattice type via a build parameter
      using latticeType = typename Traits::Lattice;
      std::shared_ptr<geometry::Domain> domainData;
      std::shared_ptr<geometry::FieldData> fieldData;
      std::shared_ptr<lb::LBM<Traits>> latticeBoltzmannModel;
      std::shared_ptr<geometry::neighbouring::NeighbouringDataManager>
        neighbouringDataManager;
      net::IOCommunicator ioComms;
      std::shared_ptr<configuration::SimConfig> simConfig;

    private:
      void Initialise();
      void SetupReporting(); // set up the reporting file
      unsigned int OutputPeriod(unsigned int frequency);
      void HandleActors();
      void OnUnstableSimulation();
      /**
       * Updates the property caches record of which properties need to be calculated
       * and cached on this iteration.
       */
      void RecalculatePropertyRequirements();

      /**
       * Helper method to log simulation parameters related to stability and accuracy
       */
      void LogStabilityReport();

      std::shared_ptr<io::PathManager> fileManager;
      reporting::Timers timings;
      std::shared_ptr<reporting::Reporter> reporter;
      reporting::BuildInfo build_info;

      std::shared_ptr<lb::SimulationState> simulationState;

      /** Struct containing the configuration of various checkers/testers */
      std::shared_ptr<lb::StabilityTester<latticeType>> stabilityTester;
      std::shared_ptr<lb::EntropyTester<latticeType>> entropyTester;
      /** Actor in charge of checking the maximum density difference across the domain */
      std::shared_ptr<lb::IncompressibilityChecker<net::PhasedBroadcastRegular<> >>
        incompressibilityChecker;

      std::shared_ptr<net::IteratedAction> cellController;
      std::shared_ptr<net::IteratedAction> colloidController;
      net::Net communicationNet;

      std::shared_ptr<util::UnitConverter> unitConverter;

      std::shared_ptr<extraction::IterableDataSource> propertyDataSource;
      std::shared_ptr<extraction::PropertyActor> propertyExtractor;

      std::shared_ptr<net::phased::StepManager> stepManager;
      std::shared_ptr<net::phased::NetConcern> netConcern;

      static constexpr LatticeTimeStep FORCE_FLUSH_PERIOD = 1000;
  };
}

#include "SimulationMaster.impl.h"
#endif /* HEMELB_SIMULATIONMASTER_H */
