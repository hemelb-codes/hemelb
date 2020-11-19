// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPIPARALLELINTEGRATIONTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPIPARALLELINTEGRATIONTESTS_H

#include <cppunit/TestFixture.h>

#include <algorithm>
#include <random>
#include <memory>

#include "redblood/CellController.h"
#include "redblood/parallel/IntegrateVelocities.h"
#include "redblood/parallel/CellParallelization.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/LatticeDataAccess.h"
#include "unittests/helpers/FolderTestFixture.h"
#include "unittests/redblood/parallel/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class MPIParallelIntegrationTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (MPIParallelIntegrationTests);
          CPPUNIT_TEST (testIntegration<hemelb::redblood::stencil::FourPoint>);
//          CPPUNIT_TEST (testIntegration<hemelb::redblood::stencil::ThreePoint>);
//          CPPUNIT_TEST (testIntegration<hemelb::redblood::stencil::CosineApprox>);
//          CPPUNIT_TEST (testIntegration<hemelb::redblood::stencil::TwoPoint>);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();

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

      void MPIParallelIntegrationTests::setUp()
      {
        FolderTestFixture::setUp();

        // Have everything ready to creates simulations
        if (net::MpiCommunicator::World().Rank() == 0)
        {
          CopyResourceToTempdir("large_cylinder_rbc.xml");
          CopyResourceToTempdir("large_cylinder.gmy");
          CopyResourceToTempdir("red_blood_cell.txt");

          // Get the RBC drop off point closer to a subdomain border so that we can test for cell communication
          // related issues with fewer timesteps.
          ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 5000);
          ModifyXMLInput("large_cylinder_rbc.xml", { "inlets", "inlet", "flowextension", "length", "value" }, 20e-6);

        }
        HEMELB_MPI_CALL(MPI_Barrier, (net::MpiCommunicator::World()));

        auto const result = Comms().Rank() == 0 ?
          "root_result" :
          "others";
        options = std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
                                                              "-in",
                                                              "large_cylinder_rbc.xml",
                                                              "-i",
                                                              "1",
                                                              "-ss",
                                                              "1111",
                                                              "-out",
                                                              result });
      }


      template<class STENCIL>
      void MPIParallelIntegrationTests::Check()
      {
        using hemelb::redblood::CellContainer;
        using hemelb::redblood::TemplateCellContainer;
        typedef typename MasterSim<STENCIL>::Traits Traits;

        auto const world = net::MpiCommunicator::World();
        auto const color = world.Rank() == 0;
        auto const split = world.Split(color);

        auto master = CreateMasterSim<STENCIL>(split);
        CPPUNIT_ASSERT(master);

        unsigned num_cells;
        auto checkNumCells = [&num_cells]( const hemelb::redblood::CellContainer & cells)
        {
          num_cells = cells.size();
        };

        auto controller = std::static_pointer_cast<hemelb::redblood::CellController<Traits>>(master->GetCellController());
        CPPUNIT_ASSERT(controller);
        controller->AddCellChangeListener(checkNumCells);

        // run the simulation
        master->RunSimulation();
        master->Finalise();

        // check that both simulations have one cell in the domain
        unsigned num_cells_sequential;
        unsigned num_cells_parallel;
        if (color)
        {
          num_cells_sequential = num_cells;
          num_cells_parallel = 0;
        }
        else
        {
          num_cells_parallel = num_cells;
        }
        world.Broadcast(num_cells_sequential, 0);
        num_cells_parallel = world.AllReduce(num_cells_parallel, MPI_SUM);
        CPPUNIT_ASSERT_EQUAL(num_cells_sequential, 1u);
        CPPUNIT_ASSERT_EQUAL(num_cells_parallel, 1u);

        // Find out number of vertices in the RBC mesh
        unsigned num_vertices = color ? (*controller->GetCells().begin())->GetVertices().size() : 0;
        world.Broadcast(num_vertices, 0);

        // Find out the location of the mesh vertices in the sequential run
        hemelb::redblood::MeshData::Vertices sequential_vertices(num_vertices);
        if (color)
        {
          sequential_vertices = (*controller->GetCells().begin())->GetVertices();
        }
        HEMELB_MPI_CALL(MPI_Bcast, (&sequential_vertices[0], 3*sequential_vertices.size(), net::MpiDataType<LatticeDistance>(), 0, net::MpiCommunicator::World()));

        // Find out the location of the mesh vertices in the parallel run
        hemelb::redblood::MeshData::Vertices vertices_location_reduction(num_vertices, {0,0,0});
        if (!color && (controller->GetCells().size()==1))
        {
          // If I'm the process in the parallel run than owns the only existing cell
          vertices_location_reduction = (*controller->GetCells().begin())->GetVertices();
        }
        hemelb::redblood::MeshData::Vertices parallel_vertices(num_vertices);
        HEMELB_MPI_CALL(MPI_Allreduce,
                        (net::MpiConstCast(&vertices_location_reduction[0]), &parallel_vertices[0], 3*vertices_location_reduction.size(), net::MpiDataType<LatticeDistance>(), MPI_SUM, net::MpiCommunicator::World()));

        // Compare locations in sequential and parallel runs
        for (auto const item : util::zip(sequential_vertices, parallel_vertices))
        {
          auto delta = 1e-12;
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item)[0], std::get<1>(item)[0], delta);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item)[1], std::get<1>(item)[1], delta);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item)[2], std::get<1>(item)[2], delta);
        }

      }

      CPPUNIT_TEST_SUITE_REGISTRATION (MPIParallelIntegrationTests);
    }
  }
}

#endif  // ONCE
