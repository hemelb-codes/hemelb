//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPI_SPREAD_FORCES_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPI_SPREAD_FORCES_TESTS_H

#include <cppunit/TestFixture.h>

#include<algorithm>
#include <random>

#include "redblood/parallel/SpreadForces.h"
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
      class MPISpreadForcesTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (MPISpreadForcesTests);
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

      void MPISpreadForcesTests::setUp()
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
      void MPISpreadForcesTests::Check(size_t mid, size_t edges, size_t nCells)
      {
        using hemelb::redblood::CellContainer;
        using hemelb::redblood::parallel::nodeDistributions;
        auto const world = net::MpiCommunicator::World();
        if (world.Size() == 1)
        {
          return;
        }
        auto const color = world.Rank() == 0;
        auto const split = world.Split(color);
        auto master = CreateMasterSim<STENCIL>(split);
        auto &latDat = master->GetLatticeData();
        helpers::ZeroOutForces(latDat);

        // Figure out positions to use for cell nodes
        auto const cells = CreateCellsFromSpecialPositions(latDat, mid, edges, world, nCells);
        auto const owned =
            split.Size() != 1 ?
              CellContainer { cells.begin() + split.Rank() * nCells, cells.begin()
                                  + (1 + split.Rank()) * nCells } :
              CellContainer { cells.begin(), cells.end() };
        auto const distributions = nodeDistributions(latDat, owned);

        hemelb::redblood::parallel::SpreadForces mpi_spreader(CreateDumbGraphComm(split));
        mpi_spreader.PostMessageLength(distributions, owned);
        mpi_spreader.ComputeForces(owned);
        mpi_spreader.PostForcesAndNodes(distributions, owned);
        mpi_spreader.SpreadLocalForces<typename MasterSim<STENCIL>::Traits>(latDat, owned);
        mpi_spreader.SpreadNonLocalForces<typename MasterSim<STENCIL>::Traits>(latDat);

        std::vector<LatticeVector> indices;
        std::vector<LatticeForceVector> forces;
        if (color)
        {
          for (std::size_t i(0); i < std::size_t(latDat.GetLocalFluidSiteCount()); ++i)
          {
            auto const site = latDat.GetSite(i);
            if (site.GetForce().GetMagnitudeSquared() > 1e-8)
            {
              indices.push_back(site.GetGlobalSiteCoords());
              forces.push_back(site.GetForce());
            }
          }
        }
        int nIndices = indices.size();
        world.Broadcast(nIndices, 0);
        CPPUNIT_ASSERT(nIndices > 0);
        indices.resize(nIndices);
        world.Broadcast(indices, 0);

        if (not color)
        {
          for (auto coords : indices)
          {
            auto const id = latDat.GetProcIdFromGlobalCoords(coords);
            if (id == split.Rank())
            {
              auto const site = latDat.GetSite(coords);
              forces.push_back(-site.GetForce());
            }
            else
            {
              forces.push_back(0);
            }
          }
        }
        CPPUNIT_ASSERT_EQUAL(size_t(nIndices), size_t(forces.size()));
        std::vector<LatticeForceVector> summed(indices.size());
        HEMELB_MPI_CALL(MPI_Allreduce,
                        (forces.data(), summed.data(), forces.size() * 3, net::MpiDataType<double>(), MPI_SUM, world));
        for (auto const force : summed)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(0, force.x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(0, force.y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(0, force.z, 1e-8);
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (MPISpreadForcesTests);
    }
  }
}

#endif  // ONCE
