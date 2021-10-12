// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>
#include <random>
#include <memory>
#include <iterator>

#include <catch2/catch.hpp>

#include "redblood/CellController.h"
#include "redblood/parallel/IntegrateVelocities.h"
#include "redblood/parallel/CellParallelization.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/redblood/parallel/ParallelFixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;
    //! Parallel and Sequential move in lock step
    class MPILockStepTests : public helpers::FolderTestFixture
    {
    public:
      MPILockStepTests();

      template<class STENCIL> void testIntegration()
      {
	Check<STENCIL>();
      }

    protected:
      std::shared_ptr<hemelb::configuration::CommandLine> options;

      //! Meta-function to create simulation type
      template<class STENCIL>
      struct MasterSim
      {
	typedef ::hemelb::Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type LBTraits;
	typedef typename LBTraits::ChangeStencil<STENCIL>::Type Traits;
	typedef OpenedSimulationMaster<Traits> Type;
      };

      //! Creates a master simulation
      template<class STENCIL>
      std::shared_ptr<typename MasterSim<STENCIL>::Type> CreateMasterSim(
									 net::MpiCommunicator const &comm) const
      {
	typedef typename MasterSim<STENCIL>::Type MasterSim;
	return std::make_shared<MasterSim>(*options, comm);
      }

      template<class STENCIL> void Check();
    };


    MPILockStepTests::MPILockStepTests(): FolderTestFixture() {
      // Have everything ready to creates simulations
      if (Comms().Rank() == 0) {
	CopyResourceToTempdir("cyl_l100_r5.xml");
	CopyResourceToTempdir("cyl_l100_r5.gmy");
	CopyResourceToTempdir("rbc_ico_2880.msh");

	// This simulation duration is sufficient to pick up the
	// original force spreading issue that motivated the test. Run
	// the test for longer in order to check other aspects of the
	// parallel implementation against a sequential run.
	ModifyXMLInput("cyl_l100_r5.xml", { "simulation", "steps", "value" }, 1000);
      }
      Comms().Barrier();

      auto const result = Comms().Rank() == 0 ?
	"root_result" :
	"others";
      options = std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
	    "-in",
	    "cyl_l100_r5.xml",
	    "-i",
	    "1",
	    "-ss",
	    "1111",
	    "-out",
	    result });
    }

    void checkBarycenter(
			 net::MpiCommunicator const &world, hemelb::redblood::CellContainer const &cells) {

      std::vector<uint64_t> uuids;
      std::vector<LatticePosition> barycenters;
      for(auto const &cell: cells) {
	uuids.push_back(*static_cast<uint64_t const*>(static_cast<void const*>(&cell->GetTag())));
	barycenters.push_back(cell->GetBarycenter());
      }
      uuids = world.Gather(uuids, 0);
      barycenters = world.Gather(barycenters, 0);
      if(world.Rank() == 0) {
	REQUIRE(uuids.size() == 2 * cells.size());
	std::map<uint64_t, LatticePosition> serial;
	for(std::size_t i(0); i < cells.size(); ++i) {
	  serial[uuids[i]] = barycenters[i];
	}
	for(std::size_t i(serial.size()); i < uuids.size(); ++i) {
	  REQUIRE(std::size_t(1) == serial.count(uuids[i]));
	  REQUIRE(ApproxV(barycenters[i]).Margin(1e-11) == serial[uuids[i]]);
	}
      }
      world.Barrier();
    }

    void checkNonZeroForceSites(
          net::MpiCommunicator const &world,
          geometry::LatticeData const &latDat,
          std::size_t &nbtests,
          hemelb::redblood::CellContainer const &cells) {
      std::vector<LatticeVector> positions;
      std::vector<LatticeForceVector> forces;
      for(site_t i(0); i < latDat.GetLocalFluidSiteCount(); ++i) {
	auto const site = latDat.GetSite(i);
	if(site.GetForce().GetMagnitude() > 1e-12) {
	  positions.push_back(site.GetGlobalSiteCoords());
	  forces.push_back(site.GetForce());
	}
      }
      if(cells.size() == 0 and world.Rank() == 0) {
	REQUIRE(std::size_t(0) == positions.size());
      }
      if(world.Rank() == 0) {
	auto const parallel_positions = world.Gather(std::vector<LatticeVector>{}, 0);
	auto const parallel_forces = world.Gather(std::vector<LatticeForceVector>{}, 0);
	REQUIRE(positions.size() == parallel_positions.size());
	for(std::size_t i(0); i < positions.size(); ++i) {
	  auto const i_found = std::find(
					 parallel_positions.begin(), parallel_positions.end(), positions[i]);
	  REQUIRE(i_found != parallel_positions.end());
	  auto const actual_force = parallel_forces[i_found - parallel_positions.begin()];
	  REQUIRE(ApproxV(forces[i]).Margin(1e-11) == actual_force);
	  ++nbtests;
	}
      } else {
	world.Gather(positions, 0);
	world.Gather(forces, 0);
      }
      world.Barrier();
    }

    template<class STENCIL>
    void MPILockStepTests::Check()
    {
      using hemelb::redblood::CellContainer;
      using hemelb::redblood::TemplateCellContainer;
      using hemelb::redblood::CellInserter;
      typedef typename MasterSim<STENCIL>::Type::Traits Traits;

      auto const world = Comms();
      auto const color = world.Rank() == 0;
      auto const split = world.Split(color);
      if(world.Size() < 3) {
	log::Logger::Log<log::Debug, log::Singleton>(
						     "Lock step tests of no interest if fewer than 3 processors");
	return;
      }

      auto master = CreateMasterSim<STENCIL>(split);
      REQUIRE(master);

      auto controller = std::static_pointer_cast<hemelb::redblood::CellController<Traits>>(
											   master->GetCellController());
      auto const originalCellInserter = controller->GetCellInsertion();
      auto rank = world.Rank();
      auto cellInserter = [&originalCellInserter, rank](CellInserter const &adder) {
	auto const transformCell = [adder, originalCellInserter, rank](CellContainer::value_type cell) {
	  static boost::uuids::uuid nbCells = {{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
						0x0, 0x0, 0x0, 0x0, 0x0, static_cast<uint8_t>(rank)}};
	  cell->SetTag(nbCells);
	  ++*static_cast<int64_t*>(static_cast<void*>(&nbCells));
	  *cell += LatticePosition(0,0,3.6);
	  adder(cell);
	};
	originalCellInserter(transformCell);
      };
      controller->SetCellInsertion(cellInserter);
      std::size_t nbtests = 0;
      controller->AddCellChangeListener(
					std::bind(
						  checkNonZeroForceSites,
						  std::cref(world), std::cref(master->GetLatticeData()), std::ref(nbtests),
						  std::placeholders::_1
						  )
					);
      controller->AddCellChangeListener(
					std::bind(
						  checkBarycenter,
						  std::cref(world),
						  std::placeholders::_1
						  )
					);

      // run the simulation
      master->RunSimulation();
      master->Finalise();
      REQUIRE((nbtests > 0 or world.Rank() != 0));
    }

    METHOD_AS_TEST_CASE(MPILockStepTests::testIntegration<hemelb::redblood::stencil::FourPoint>,
			"Test running in lockstep with four point stencil",
			"[redblood][.long]");
    METHOD_AS_TEST_CASE(MPILockStepTests::testIntegration<hemelb::redblood::stencil::ThreePoint>,
			"Test running in lockstep with three point stencil",
			"[redblood][.long]");
    METHOD_AS_TEST_CASE(MPILockStepTests::testIntegration<hemelb::redblood::stencil::CosineApprox>,
			"Test running in lockstep with cosine stencil",
			"[redblood][.long]");
    METHOD_AS_TEST_CASE(MPILockStepTests::testIntegration<hemelb::redblood::stencil::TwoPoint>,
			"Test running in lockstep with two point stencil",
			"[redblood][.long]");
  }
}

