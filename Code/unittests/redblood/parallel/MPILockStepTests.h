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
          CPPUNIT_TEST (testIntegration<hemelb::redblood::stencil::ThreePoint>);
          CPPUNIT_TEST (testIntegration<hemelb::redblood::stencil::CosineApprox>);
          CPPUNIT_TEST (testIntegration<hemelb::redblood::stencil::TwoPoint>);
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
        if (Comms().Rank() == 0)
        {
          CopyResourceToTempdir("large_cylinder_rbc.xml");
          CopyResourceToTempdir("large_cylinder.gmy");
          CopyResourceToTempdir("red_blood_cell.txt");

          // This simulation duration is sufficient to pick up the original force exchange
          // issue that motivated the test. Run the test for longer in order to check other
          // aspects of the parallel implementation against a sequential run.
          ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 1000);
        }
        HEMELB_MPI_CALL(MPI_Barrier, (Comms()));

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
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].x, serial[uuids[i]].x, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].y, serial[uuids[i]].y, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenters[i].z, serial[uuids[i]].z, 1e-12);
          }
        }
        world.Barrier();
      }

      void checkNonZeroForceSites(
          net::MpiCommunicator const &world,
          geometry::LatticeData const &latDat,
          std::size_t &nbtests,
          hemelb::redblood::CellContainer const &cells) {
        std::vector<LatticePosition> positions;
        std::vector<LatticeForceVector> forces;
        for(site_t i(0); i < latDat.GetLocalFluidSiteCount(); ++i)
        {
          auto const site = latDat.GetSite(i);
          if(site.GetForce().GetMagnitude() > 1e-12)
          {
            positions.push_back(site.GetGlobalSiteCoords());
            forces.push_back(site.GetForce());
          }
        }
        if(cells.size() == 0 and world.Rank() == 0)
        {
          CPPUNIT_ASSERT_EQUAL(std::size_t(0), positions.size());
        }
        if(world.Rank() == 0)
        {
          auto const parallel_positions = world.Gather(std::vector<LatticePosition>{}, 0);
          auto const parallel_forces = world.Gather(std::vector<LatticeForceVector>{}, 0);
          CPPUNIT_ASSERT_EQUAL(positions.size(), parallel_positions.size());
          for(std::size_t i(0); i < positions.size(); ++i)
          {
            auto const i_found = std::find(
                parallel_positions.begin(), parallel_positions.end(), positions[i]);
            CPPUNIT_ASSERT(i_found != parallel_positions.end());
            auto const actual_force = parallel_forces[i_found - parallel_positions.begin()];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(forces[i].x, actual_force.x, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(forces[i].y, actual_force.y, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(forces[i].z, actual_force.z, 1e-12);
            ++nbtests;
          }
        }
        else
        {
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

        // run the simulation
        master->RunSimulation();
        master->Finalise();
        CPPUNIT_ASSERT(nbtests > 0 or world.Rank() != 0);
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (MPILockStepTests);
    }
  }
}

#endif  // ONCE
