//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPI_INTEGRATE_VELOCITIES_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPI_INTEGRATE_VELOCITIES_TESTS_H

#include <cppunit/TestFixture.h>

#include<algorithm>
#include <random>

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
      class MPIIntegrateVelocitiesTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (MPIIntegrateVelocitiesTests);
          CPPUNIT_TEST (testMidRegion<hemelb::redblood::stencil::FourPoint> );
          CPPUNIT_TEST (testEdgeRegion<hemelb::redblood::stencil::FourPoint> );
          CPPUNIT_TEST (testAll<hemelb::redblood::stencil::FourPoint> );
          CPPUNIT_TEST (testMidRegion<hemelb::redblood::stencil::ThreePoint> );
          CPPUNIT_TEST (testEdgeRegion<hemelb::redblood::stencil::ThreePoint> );
          CPPUNIT_TEST (testAll<hemelb::redblood::stencil::ThreePoint> );
          CPPUNIT_TEST (testMidRegion<hemelb::redblood::stencil::CosineApprox> );
          CPPUNIT_TEST (testEdgeRegion<hemelb::redblood::stencil::CosineApprox> );
          CPPUNIT_TEST (testAll<hemelb::redblood::stencil::CosineApprox> );
          CPPUNIT_TEST (testMidRegion<hemelb::redblood::stencil::TwoPoint> );
          CPPUNIT_TEST (testEdgeRegion<hemelb::redblood::stencil::TwoPoint> );
          CPPUNIT_TEST (testAll<hemelb::redblood::stencil::TwoPoint> );CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();

          template<class STENCIL> void testMidRegion()
          {
            Check<STENCIL>(2, 0, 1);
          }
          template<class STENCIL> void testEdgeRegion()
          {
            Check<STENCIL>(0, 2, 1);
          }
          template<class STENCIL> void testAll()
          {
            Check<STENCIL>(3, 3, 2);
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

          template<class STENCIL> void Check(size_t mid, size_t edges, size_t nCells);
      };

      void MPIIntegrateVelocitiesTests::setUp()
      {
        using hemelb::configuration::CommandLine;
        FolderTestFixture::setUp();

        // Have everything ready to creates simulations
        if (net::MpiCommunicator::World().Rank() == 0)
        {
          CopyResourceToTempdir("red_blood_cell.txt");
          TiXmlDocument doc(resources::Resource("large_cylinder.xml").Path());
          CopyResourceToTempdir("large_cylinder.xml");
          ModifyXMLInput("large_cylinder.xml", { "simulation", "steps", "value" }, 2);
          CopyResourceToTempdir("large_cylinder.gmy");
        }
        HEMELB_MPI_CALL(MPI_Barrier, (net::MpiCommunicator::World()));
        auto const result = Comms().Rank() == 0 ?
          "root_result" :
          "others";
        options = std::make_shared<CommandLine>(CommandLine { "hemelb",
                                                              "-in",
                                                              "large_cylinder.xml",
                                                              "-i",
                                                              "1",
                                                              "-ss",
                                                              "1111",
                                                              "-out",
                                                              result });
      }

      template<class STENCIL>
      void MPIIntegrateVelocitiesTests::Check(size_t mid, size_t edges, size_t nCells)
      {
        using hemelb::redblood::CellContainer;
        using hemelb::redblood::TemplateCellContainer;
        typedef typename MasterSim<STENCIL>::Traits Traits;
        auto const world = net::MpiCommunicator::World();
        auto const color = world.Rank() == 0;
        auto const split = world.Split(color);
        auto master = CreateMasterSim<STENCIL>(split);
        auto &latDat = master->GetLatticeData();
        helpers::ZeroOutForces(latDat);

        // Setup same distribution for all instances
        for (Direction i = 0; i < Traits::Lattice::NUMVECTORS; ++i)
        {
          auto func = [i](LatticePosition const &x)
          {
            double const fac = (1e0 + static_cast<double>(i) * 0.1);
            return x.Dot(x) * fac - x.Dot(LatticePosition(1, 1, 1));
          };
          helpers::setUpDistribution<typename Traits::Lattice>(&latDat, i, func);
        }

        // Figure out positions to use for cell nodes
        auto const cells = CreateCellsFromSpecialPositions(latDat, mid, edges, world, nCells);
        auto const checkMoves = cells[0]->GetVertices()[0];
        auto owned =
            split.Size() != 1 ?
              CellContainer { cells.begin() + split.Rank() * nCells, cells.begin()
                                  + (1 + split.Rank()) * nCells } :
              CellContainer { cells.begin(), cells.end() };
        auto const distributions = hemelb::redblood::parallel::nodeDistributions(latDat, owned);

        // Goes through "ExchangeCells" to figure out who owns/lends what.
        // Ownership is pre-determined here: first nCells got to 0, second nCells to 2, etc...
        TemplateCellContainer const templates =
            { { cells[0]->GetTemplateName(), cells[0]->clone() } };
        auto ownership = [&cells, nCells](CellContainer::const_reference cell)
        {
          size_t i(0);
          for(auto const& c: cells)
          {
            if(cell->GetTag() == c->GetTag())
            {
              break;
            }
            ++i;
          }
          proc_t result = i / nCells;
          return result;
        };
        hemelb::redblood::parallel::ExchangeCells xchange(CreateDumbGraphComm(split), split);
        xchange.PostCellMessageLength(distributions, owned, ownership);
        xchange.PostCells(distributions, owned, ownership);
        auto const distCells = xchange.ReceiveCells(templates);
        auto const &lentCells = std::get<2>(distCells);

        // Actually perform velocity integration
        hemelb::redblood::parallel::IntegrateVelocities integrator(CreateDumbGraphComm(split));
        integrator.PostMessageLength(lentCells);
        integrator.ComputeLocalVelocitiesAndUpdatePositions<Traits>(latDat, owned);
        integrator.PostVelocities<Traits>(latDat, lentCells);
        integrator.UpdatePositionsNonLocal(distributions, owned);

        if (world.Rank() == 0)
        {
          CPPUNIT_ASSERT( (checkMoves - cells[0]->GetVertices()[0]).GetMagnitude() > 1e-8);
        }

        for (auto const &cell : cells)
        {
          typedef std::vector<LatticePosition> Positions;
          Positions positions = world.Rank() == 0 ?
            cell->GetVertices() :
            Positions(cell->GetNumberOfNodes(), 0);
          world.Broadcast(positions, 0);
          if (world.Rank() != 0 and owned.find(cell) != owned.end())
          {
            for (auto const item : util::zip(cell->GetVertices(), positions))
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<1>(item).x, std::get<0>(item).x, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<1>(item).y, std::get<0>(item).y, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<1>(item).z, std::get<0>(item).z, 1e-8);
            }
          }
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (MPIIntegrateVelocitiesTests);
    }
  }
}

#endif  // ONCE
