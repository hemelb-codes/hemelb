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
#include <iterator>

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
      //! Parallel and Sequential move in lock step
      class MPILockStepTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (MPILockStepTests);
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


      void MPILockStepTests::setUp()
      {
        FolderTestFixture::setUp();

        // Have everything ready to creates simulations
        if (net::MpiCommunicator::World().Rank() == 0)
        {
          CopyResourceToTempdir("large_cylinder_rbc.xml");
          CopyResourceToTempdir("large_cylinder.gmy");
          CopyResourceToTempdir("red_blood_cell.txt");

          ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 20000);
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

      void checkBarycenter(
          net::MpiCommunicator const &world, hemelb::redblood::CellContainer const &cells) {

        std::vector<uint64_t> uuids;
        std::vector<LatticePosition> barycenters;
        for(auto const &cell: cells)
        {
          uuids.push_back(*static_cast<uint64_t const*>(static_cast<void const*>(&cell->GetTag())));
          barycenters.push_back(cell->GetBarycenter());
        }
        uuids = world.Gather(uuids, 0);
        barycenters = world.Gather(barycenters, 0);
        if(world.Rank() == 0)
        {
          CPPUNIT_ASSERT_EQUAL(uuids.size(), 2 * cells.size());
          std::map<uint64_t, LatticePosition> serial;
          for(std::size_t i(0); i < cells.size(); ++i)
          {
            serial[uuids[i]] = barycenters[i];
          }
          for(std::size_t i(serial.size()); i < uuids.size(); ++i)
          {
            CPPUNIT_ASSERT_EQUAL(std::size_t(1), serial.count(uuids[i]));
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].x, serial[uuids[i]].x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].y, serial[uuids[i]].y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].z, serial[uuids[i]].z, 1e-8);
          }
        }
      }

      void checkCellNodes(
          net::MpiCommunicator const &world, hemelb::redblood::CellContainer const &cells) {

        std::vector<uint64_t> uuids;
        std::vector<unsigned> nNodes;
        std::vector<LatticePosition> nodes;
        for(auto const &cell: cells)
        {
          uuids.push_back(*static_cast<uint64_t const*>(static_cast<void const*>(&cell->GetTag())));
          nNodes.push_back(cells.size());
          for(auto const& cell: cells)
          {
            std::copy(cell->GetVertices().begin(), cell->GetVertices().end(),
                std::back_inserter(nodes));
          }
        }
        uuids = world.Gather(uuids, 0);
        nNodes = world.Gather(nNodes, 0);
        nodes = world.Gather(nodes, 0);
        if(world.Rank() == 0)
        {
          CPPUNIT_ASSERT_EQUAL(uuids.size(), 2 * cells.size());
          std::map<uint64_t, unsigned> serial;
          for(std::size_t i(0); i < cells.size(); ++i)
          {
            serial[uuids[i]] = barycenters[i];
          }
          for(std::size_t i(serial.size()); i < uuids.size(); ++i)
          {
          std::cout << "YEAH I AM THERE " << barycenters[i] << std::endl;
            CPPUNIT_ASSERT_EQUAL(std::size_t(1), serial.count(uuids[i]));
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].x, serial[uuids[i]].x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].y, serial[uuids[i]].y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].z, serial[uuids[i]].z, 1e-8);
          }
        }
      }


      template<class STENCIL>
      void MPILockStepTests::Check()
      {
        using hemelb::redblood::CellContainer;
        using hemelb::redblood::TemplateCellContainer;
        using hemelb::redblood::CellInserter;
        typedef typename MasterSim<STENCIL>::Type::Traits Traits;

        auto const world = net::MpiCommunicator::World();
        auto const color = world.Rank() == 0;
        auto const split = world.Split(color);
        if(world.Size() < 3)
        {
          log::Logger::Log<log::Debug, log::Singleton>(
              "Lock step tests of no interest if fewer than 3 processors");
          return;
        }

        auto master = CreateMasterSim<STENCIL>(split);
        CPPUNIT_ASSERT(master);

        auto controller = std::static_pointer_cast<hemelb::redblood::CellController<Traits>>(
            master->GetCellController());
        auto const originalCellInserter = controller->GetCellInsertion();
        auto cellInserter = [&originalCellInserter](CellInserter const &adder) {
          auto const transformCell = [adder, originalCellInserter](CellContainer::value_type cell) {
            static boost::uuids::uuid nbCells = {{0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
              0x0, 0x0, 0x0, 0x0, 0x0, 0x0}};
            cell->SetTag(nbCells);
            ++*static_cast<int64_t*>(static_cast<void*>(&nbCells));
            *cell += LatticePosition(0,0,3);
            adder(cell);
          };
          originalCellInserter(transformCell);
        };
        controller->SetCellInsertion(cellInserter);
        controller->AddCellChangeListener(std::bind(checkBarycenter, world, std::placeholders::_1));
        auto const checkLent = [controller](CellContainer const &) {
          if(controller->GetLentCells().size() > 0)
          {
            std::cout << "doest have lent cells" << std::endl;
          }
        };
        controller->AddCellChangeListener(checkLent);

        // run the simulation
        master->RunSimulation();
        master->Finalise();
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (MPILockStepTests);
    }
  }
}

#endif  // ONCE
