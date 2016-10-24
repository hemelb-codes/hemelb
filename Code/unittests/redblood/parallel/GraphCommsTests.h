//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_GRAPHCOMMSTEST_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_GRAPHCOMMSTEST_H

#include "unittests/helpers/FolderTestFixture.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/redblood/parallel/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class GraphCommsTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (GraphCommsTests);
          CPPUNIT_TEST (testDumbGraphCommunicator);
          CPPUNIT_TEST (testGraphCommunicator);
          CPPUNIT_TEST (testComputeCellsEffectiveSize);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();
          void testDumbGraphCommunicator();
          void testGraphCommunicator();
          void testComputeCellsEffectiveSize();

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

      void GraphCommsTests::setUp()
      {
        FolderTestFixture::setUp();

        // Have everything ready to creates simulations
        if (Comms().Rank() == 0)
        {
          CopyResourceToTempdir("large_cylinder_rbc.xml");
          CopyResourceToTempdir("large_cylinder.gmy");
          CopyResourceToTempdir("red_blood_cell.txt");

          // This simulation duration is sufficient to pick up the original force spreading
          // issue that motivated the test. Run the test for longer in order to check other
          // aspects of the parallel implementation against a sequential run.
          ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 1000);
        }
        HEMELB_MPI_CALL(MPI_Barrier, (Comms()));

        options =
            std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
                                                                                      "-in",
                                                                                      "large_cylinder_rbc.xml",
                                                                                      "-i",
                                                                                      "1",
                                                                                      "-ss",
                                                                                      "1111" });
      }

      void GraphCommsTests::testDumbGraphCommunicator()
      {
        using hemelb::redblood::ComputeProcessorNeighbourhood;

        auto comms = Comms();
        auto neighbourhoods = ComputeProcessorNeighbourhood(comms);

        for (auto procid_neighbours : util::enumerate(neighbourhoods))
        {
          // All other processes must be neighbours
          CPPUNIT_ASSERT_EQUAL(int(procid_neighbours.value.size()), comms.Size() - 1);

          for (auto neighbour = 0; neighbour < comms.Size(); ++neighbour)
          {
            if (static_cast<unsigned>(neighbour) != procid_neighbours.index)
            {
              // All other process indices must be in the std::vector
              CPPUNIT_ASSERT(find(procid_neighbours.value.begin(),
                                  procid_neighbours.value.end(),
                                  neighbour) != procid_neighbours.value.end());
            }
          }
        }
      }

      void GraphCommsTests::testGraphCommunicator()
      {
        using hemelb::redblood::ComputeProcessorNeighbourhood;

        // Test only makes sense if run with 4 cores
        auto comms = Comms();
        CPPUNIT_ASSERT_EQUAL(4, comms.Size());

        // Setup simulation with cylinder
        auto master = CreateMasterSim(comms);
        CPPUNIT_ASSERT(master);
        auto &latticeData = master->GetLatticeData();

        // Compute neighbourhoods (cylinder is 4.8e-5 long, a cell effective size of 2e-6 won't let cells span across more than two subdomains)
        auto neighbourhoods =
            ComputeProcessorNeighbourhood(comms,
                                          latticeData,
                                          2e-6 / master->GetSimConfig()->GetVoxelSize());

        // Parmetis divides the cylinder in four consecutive cylindrical subdomains with interfaces roughly parallel to the iolets.
        // The ranks are ordered 0,3,1,2 along the positive direction of the z axis.
        std::vector<std::vector<int>> expected_neighbourhoods;
        expected_neighbourhoods.push_back( { 3 });
        expected_neighbourhoods.push_back( { 2, 3 });
        expected_neighbourhoods.push_back( { 1 });
        expected_neighbourhoods.push_back( { 0, 1 });

        // Compare computed vs expected neighbourhoods
        for (auto procid_neighbours : util::enumerate(neighbourhoods))
        {
          CPPUNIT_ASSERT(procid_neighbours.value
              == expected_neighbourhoods[procid_neighbours.index]);
        }
      }

      void GraphCommsTests::testComputeCellsEffectiveSize()
      {
        using hemelb::redblood::ComputeCellsEffectiveSize;

        // Setup simulation with cylinder
        auto comms = Comms();
        auto master = CreateMasterSim(comms);
        CPPUNIT_ASSERT(master);

        // Biggest cell radius in lattice units times a tolerance
        CPPUNIT_ASSERT_DOUBLES_EQUAL(hemelb::redblood::EFFECTIVE_SIZE_TO_RADIUS_RATIO
                                         * (8e-06 / master->GetSimConfig()->GetVoxelSize()),
                                     ComputeCellsEffectiveSize(master->GetSimConfig()->GetRBCMeshes()),
                                     1e-9);

      }

      CPPUNIT_TEST_SUITE_REGISTRATION (GraphCommsTests);
    }
  }
}

#endif
