// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "redblood/RBCConfig.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/redblood/Fixtures.h"
#include "tests/redblood/parallel/ParallelFixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    class GraphCommsTests : public helpers::FolderTestFixture
    {
    public:
      GraphCommsTests();

      void testLexicographicalOrderingOfVectors();
      void testDumbGraphCommunicator();
      void testGraphCommunicator();
      void testComputeCellsEffectiveSize();
      void testComputeGlobalCoordsToProcMap();

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
      template<class STENCIL = hemelb::redblood::stencil::FourPoint>
      std::shared_ptr<typename MasterSim<STENCIL>::Type> CreateMasterSim(
              net::MpiCommunicator const &comm) const
      {
	typedef typename MasterSim<STENCIL>::Type MasterSim;
	return std::make_shared<MasterSim>(*options, comm);
      }

    };

    GraphCommsTests::GraphCommsTests() : FolderTestFixture() {
      // Have everything ready to creates simulations
      if (Comms().Rank() == 0) {
	CopyResourceToTempdir("large_cylinder_rbc.xml");
	CopyResourceToTempdir("large_cylinder.gmy");
	CopyResourceToTempdir("red_blood_cell.txt");

	// This simulation duration is sufficient to pick up the original force spreading
	// issue that motivated the test. Run the test for longer in order to check other
	// aspects of the parallel implementation against a sequential run.
	ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 1000);
      }
      Comms().Barrier();

      options =
	std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
	      "-in",
	      "large_cylinder_rbc.xml",
    });
    }

    void GraphCommsTests::testLexicographicalOrderingOfVectors() {
      auto cmp = hemelb::redblood::parallel::VectorLexicographicalOrdering<util::Vector3D<int>>{};

      util::Vector3D<int> foo{0, 0, 0};
      util::Vector3D<int> bar{1, 0, 0};

      // // Catch implementation
      // CHECK(cmp(foo, bar));

      // bar = {0, 1, 0};
      // CHECK(cmp(foo, bar));

      // bar = {0, 0, 1};
      // CHECK(cmp(foo, bar));

      // bar = {0, 0, 0};
      // CHECK(!cmp(foo, bar));

      // foo = {0, 0, 1};
      // CHECK(!cmp(foo, bar));
      REQUIRE(cmp(foo, bar));

      bar = {0, 1, 0};
      REQUIRE(cmp(foo, bar));

      bar = {0, 0, 1};
      REQUIRE(cmp(foo, bar));

      bar = {0, 0, 0};
      REQUIRE(!cmp(foo, bar));

      foo = {0, 0, 1};
      REQUIRE(!cmp(foo, bar));
    }

    void GraphCommsTests::testDumbGraphCommunicator()
    {
      using hemelb::redblood::parallel::ComputeProcessorNeighbourhood;

      auto comms = Comms();
      auto neighbourhoods = ComputeProcessorNeighbourhood(comms);

      for (auto procid_neighbours : util::enumerate(neighbourhoods)) {
	// All other processes must be neighbours
	REQUIRE(int(procid_neighbours.value.size()) == comms.Size() - 1);

	for (auto neighbour = 0; neighbour < comms.Size(); ++neighbour) {
	  if (static_cast<unsigned>(neighbour) != procid_neighbours.index) {
	    // All other process indices must be in the std::vector
	    REQUIRE(find(procid_neighbours.value.begin(),
			 procid_neighbours.value.end(),
			 neighbour) != procid_neighbours.value.end());
	  }
	}
      }
    }

    void GraphCommsTests::testGraphCommunicator()
    {
      using hemelb::redblood::parallel::ComputeProcessorNeighbourhood;

      // Test only makes sense if run with 4 cores
      auto comms = Comms();
      REQUIRE(4 == comms.Size());

      // Setup simulation with cylinder
      auto master = CreateMasterSim(comms);
      REQUIRE(master);
      auto &latticeData = master->GetLatticeData();

      // Compute neighbourhoods (cylinder is 4.8e-5 long, a cell
      // effective size of 2e-6 won't let cells span across more than
      // two subdomains)
      auto neighbourhoods =
	ComputeProcessorNeighbourhood(comms,
				      latticeData,
				      2e-6 / master->GetSimConfig()->GetVoxelSize());

      // Parmetis divides the cylinder in four consecutive cylindrical
      // subdomains with interfaces roughly parallel to the iolets.
      // The ranks are ordered 0,3,1,2 along the positive direction of
      // the z axis.
      std::vector<std::vector<int>> expected_neighbourhoods;
      expected_neighbourhoods.push_back( { 3 });
      expected_neighbourhoods.push_back( { 2, 3 });
      expected_neighbourhoods.push_back( { 1 });
      expected_neighbourhoods.push_back( { 0, 1 });

      // Compare computed vs expected neighbourhoods
      for (auto procid_neighbours : util::enumerate(neighbourhoods)) {
	REQUIRE(procid_neighbours.value
		== expected_neighbourhoods[procid_neighbours.index]);
      }
    }

    void GraphCommsTests::testComputeCellsEffectiveSize() {
      using hemelb::redblood::parallel::ComputeCellsEffectiveSize;

      // Setup simulation with cylinder
      auto comms = Comms();
      auto master = CreateMasterSim(comms);
      REQUIRE(master);

      auto simConf = master->GetSimConfig();
      REQUIRE(simConf->HasRBCSection());
      auto rbcConf = simConf->GetRBCConfig();
      REQUIRE(rbcConf);
      // Biggest cell radius in lattice units times a tolerance
      REQUIRE(
        Approx(
	  parallel::MAXIMUM_SIZE_TO_RADIUS_RATIO * (8e-06 / simConf->GetVoxelSize())
	).margin(1e-9) ==
	ComputeCellsEffectiveSize(rbcConf->GetRBCMeshes())
      );
    }

    void GraphCommsTests::testComputeGlobalCoordsToProcMap()
    {
      using hemelb::redblood::parallel::ComputeGlobalCoordsToProcMap;
      using hemelb::redblood::parallel::ComputeProcessorNeighbourhood;

      // Test only makes sense if run with 4 cores
      auto comms = Comms();
      REQUIRE(4 == comms.Size());

      // Setup simulation with cylinder
      auto master = CreateMasterSim(comms);
      REQUIRE(master);
      auto &latticeData = master->GetLatticeData();
      auto graphComm = comms.Graph(hemelb::redblood::parallel::ComputeProcessorNeighbourhood(comms,
											     latticeData,
											     2e-6 / master->GetSimConfig()->GetVoxelSize()));
      auto const& globalCoordsToProcMap = hemelb::redblood::parallel::ComputeGlobalCoordsToProcMap(graphComm, latticeData);

      // The first lattice site for each rank is
      std::vector<LatticeVector> lattice_coords;
      lattice_coords.push_back({2,6,2});
      lattice_coords.push_back({2,6,19});
      lattice_coords.push_back({2,6,29});
      lattice_coords.push_back({2,6,9});

      // Parmetis divides the cylinder in four consecutive cylindrical subdomains.
      // The ranks are ordered 0,3,1,2 along the positive direction of the z axis.
      // The cell size is such that only ranks to the left and right are neighbours.
      std::vector<std::vector<proc_t>> neighs_to_know_about;
      neighs_to_know_about.push_back({0, 3});
      neighs_to_know_about.push_back({1, 2, 3});
      neighs_to_know_about.push_back({1, 2});
      neighs_to_know_about.push_back({0, 1, 3});

      // Each process should know who is the owner of a local lattice and a lattice belonging to each neighbour
      for (auto neigh : neighs_to_know_about[comms.Rank()]) {
	// Use map::at to throw when the lattice coordinates are not known
	REQUIRE(globalCoordsToProcMap.at(lattice_coords[neigh]) == neigh);
      }
    }

    METHOD_AS_TEST_CASE(GraphCommsTests::testLexicographicalOrderingOfVectors,
			"testLexicographicalOrderingOfVectors",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testDumbGraphCommunicator,
			"testDumbGraphCommunicator",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testGraphCommunicator,
			"testGraphCommunicator",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testComputeCellsEffectiveSize,
			"testComputeCellsEffectiveSize",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testComputeGlobalCoordsToProcMap,
			"testComputeGlobalCoordsToProcMap",
			"[redblood]");
  }
}
