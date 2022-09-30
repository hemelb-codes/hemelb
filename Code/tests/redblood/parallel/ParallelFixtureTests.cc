// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>
#include <random>
#include <memory>
#include <iterator>

#include <catch2/catch.hpp>

#include "redblood/parallel/IntegrateVelocities.h"
#include "redblood/parallel/NodeCharacterizer.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/redblood/parallel/ParallelFixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;
    //! Parallel and Sequential move in lock step
    class ParallelFixtureTests : public helpers::FolderTestFixture
    {
    public:
      ParallelFixtureTests();

      //! if owner procs thinks position affects proc i, then proc i knows it as well
      void testTransitiveOwnership();
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
    };


    ParallelFixtureTests::ParallelFixtureTests() : FolderTestFixture()
    {
      // Have everything ready to creates simulations
      if (Comms().Rank() == 0)
	{
	  CopyResourceToTempdir("large_cylinder_rbc.xml");
	  CopyResourceToTempdir("large_cylinder.gmy");
	  CopyResourceToTempdir("red_blood_cell.txt");

	  ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 20000);
	}
      Comms().Barrier();

      options = std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
	    "-in",
	    "large_cylinder_rbc.xml",
	    "-out",
	    "" });
    }

    void ParallelFixtureTests::testTransitiveOwnership()
    {
      typedef hemelb::redblood::stencil::FourPoint Stencil;

      auto const world = Comms();
      if(world.Size() < 2)
        {
          return;
        }

      auto master = CreateMasterSim<Stencil>(world);
      REQUIRE(master);

      auto& fd = master->GetFieldData();
      auto& dom = fd.GetDomain();
      helpers::ZeroOutForces(fd);
      auto const nmid = 20;
      auto const nedges = 20;
      auto const positions = GatherSpecialPositions(dom, nmid, nedges, world);

      auto graphComm =
	world.Graph(hemelb::redblood::parallel::ComputeProcessorNeighbourhood(world,
                                                                          dom,
									      2e-6 / master->GetSimConfig()->GetVoxelSize()));

      auto const& globalCoordsToProcMap = hemelb::redblood::parallel::ComputeGlobalCoordsToProcMap(graphComm, dom);

      for(std::size_t i(0); i < positions.size(); ++i)
        {
          auto const procs = hemelb::redblood::parallel::details::positionAffectsProcs<Stencil>(
												globalCoordsToProcMap, positions[i].as<LatticeDistance>());

          // Send set of affected procs as known by owner proc
          decltype(world.Rank()) positions_are_from_this_proc = i / (nmid + nedges) + 1;
          int N(procs.size());
          world.Broadcast(N, positions_are_from_this_proc);
          std::vector<proc_t> expected(N);
          if(world.Rank() == positions_are_from_this_proc)
	    {
	      std::copy(procs.begin(), procs.end(), expected.begin());
	    }
          world.Broadcast(expected, positions_are_from_this_proc);

          std::set<proc_t> expected_set(expected.begin(), expected.end());

          // Owner knows thinks current position affects this proc
          if(expected_set.count(world.Rank()))
	    {
	      // so this proc must know that current positions affects it
	      REQUIRE(procs.count(world.Rank()));
	      // and affects owner proc
	      REQUIRE(procs.count(positions_are_from_this_proc));
	    }
          else
	    {
	      // otherwise, this proc should not think it is affected
	      REQUIRE(not procs.count(world.Rank()));
	    }

          // similarly, if this proc thinks current position affects owner proc
          if(not procs.count(positions_are_from_this_proc))
	    {
	      // then owner proc must think current position affects this proc
	      REQUIRE(not expected_set.count(world.Rank()));
	    }

        }
    }

    METHOD_AS_TEST_CASE(ParallelFixtureTests::testTransitiveOwnership,
			"Test transitive ownership",
			"[redblood]");
  }
}

