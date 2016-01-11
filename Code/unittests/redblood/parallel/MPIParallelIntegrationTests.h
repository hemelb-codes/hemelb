//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPIPARALLELINTEGRATIONTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPIPARALLELINTEGRATIONTESTS_H

#include <cppunit/TestFixture.h>

#include <algorithm>
#include <random>
#include <memory>

#include "redblood/parallel/IntegrateVelocities.h"
#include "redblood/parallel/CellParallelization.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/LatticeDataAccess.h"
#include "unittests/helpers/FolderTestFixture.h"
#include "unittests/redblood/parallel/Fixtures.h"
#include "redblood/CellController.h"

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

          ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 5000);
          ModifyXMLInput("large_cylinder_rbc.xml",
                         { "redbloodcells", "controller", "stencil" },
                         "two");
          ModifyXMLInput("large_cylinder_rbc.xml", { "outlets",
                                                     "outlet",
                                                     "flowextension",
                                                     "length",
                                                     "value" },
                         40e-6);
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

        unsigned num_cells_sequential;
        unsigned num_cells_parallel;
        if (world.Rank() == 0)
        {
          num_cells_sequential = num_cells;
          num_cells_parallel = 0;
        }
        else
        {
          num_cells_parallel = num_cells;
        }

        world.Broadcast(num_cells_sequential, 0);
        world.AllReduce(num_cells_parallel, MPI_SUM);

        CPPUNIT_ASSERT_EQUAL(num_cells_sequential, num_cells_parallel);

      }

      CPPUNIT_TEST_SUITE_REGISTRATION (MPIParallelIntegrationTests);
    }
  }
}

#endif  // ONCE
