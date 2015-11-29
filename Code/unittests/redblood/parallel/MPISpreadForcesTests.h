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

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      using hemelb::redblood::CellContainer;
      using hemelb::redblood::parallel::nodeDistributions;
      //! Make some functionality available
      template<class TRAITS>
      class OpenedSimulationMaster : public SimulationMaster<TRAITS> 
      {
        public:
          using SimulationMaster<TRAITS>::SimulationMaster;
          using SimulationMaster<TRAITS>::Finalise;

          void DoTimeStep()
          {
            SimulationMaster<TRAITS>::DoTimeStep();
          }
      };

      //! Open SpreadForces object to introspection
      class OpenedSpreadForces : public hemelb::redblood::parallel::SpreadForces
      {
        public:
          using hemelb::redblood::parallel::SpreadForces::SpreadForces;
#         define HEMELB_MACRO(Name, name, TYPE)   \
            TYPE & Get ## Name()                  \
            {                                     \
              return name;                        \
            }

            HEMELB_MACRO(SendNodeCount, sendNodeCount, net::INeighborAllToAll<int>);
            HEMELB_MACRO(SendPositions, sendPositions, net::INeighborAllToAllV<LatticePosition>);
            HEMELB_MACRO(SendForces, sendForces, net::INeighborAllToAllV<LatticeForceVector>);
            HEMELB_MACRO(
                CellForces, cellForces, hemelb::redblood::parallel::SpreadForces::CellForces);
#         undef HEMELB_MACRO
      };

      class DummyCell : public NodeCell
      {
        public:
          LatticeForce force;

          DummyCell(
              LatticePosition const&position, LatticeForce f = 0e0,
              std::string const &templateName = "nope")
            : DummyCell(std::vector<LatticePosition>{position}, f, templateName)
          {
          }
          DummyCell(
              std::vector<LatticePosition> const &positions, LatticeForce f = 0e0,
              std::string const &templateName = "nope")
            : NodeCell(positions, templateName), force(f)
          {
          }
          template<class ITER>
          DummyCell(ITER first, ITER last,
              LatticeForce f = 0e0, std::string const &templateName = "nope")
            : NodeCell(first, last, templateName), force(f)
          {
          }

          LatticeEnergy operator()(std::vector<LatticeForceVector> &f) const override
          {
            f.resize(GetNumberOfNodes());
            std::fill(f.begin(), f.end(), force);
            return 0e0;
          }
          std::unique_ptr<CellBase> cloneImpl() const override
          {
            return std::unique_ptr<DummyCell>{new DummyCell(*this)};
          }
      };

      class MPISpreadForcesTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (MPISpreadForcesTests);
          CPPUNIT_TEST (testMidregion<hemelb::redblood::stencil::FourPoint>);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();

          //! Tests that forces from a single node is spread to other processes
          template<class STENCIL> void testMidregion();


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

          //! \brief gathers mid-domain and egde positions from all procs
          //! \details If there are insufficient number of edges, mid-domains are used instead.
          //! erase removes the components from the first process.
          std::vector<LatticePosition> GetPositions(
            geometry::LatticeData const & latDat,
            size_t mid, size_t edges, net::MpiCommunicator const &c, bool erase=true);
          //! Creates list of cells for each set of positions from each process
          std::vector<CellContainer::value_type> GetCells(
            geometry::LatticeData const & latDat,
            size_t mid, size_t edges, net::MpiCommunicator const &c, bool erase=true);

          template<class STENCIL> void Check(size_t mid, size_t edges);

          net::MpiCommunicator GetGraphComm(net::MpiCommunicator const &comm);
      };

      void MPISpreadForcesTests::setUp()
      {
        using hemelb::configuration::CommandLine;
        FolderTestFixture::setUp();

        // Have everything ready to creates simulations
        if(net::MpiCommunicator::World().Rank() == 0)
        {
          CopyResourceToTempdir("red_blood_cell.txt");
          TiXmlDocument doc(resources::Resource("large_cylinder.xml").Path());
          CopyResourceToTempdir("large_cylinder.xml");
          ModifyXMLInput("large_cylinder.xml", {"simulation", "steps", "value"}, 2);
          CopyResourceToTempdir("large_cylinder.gmy");
        }
        HEMELB_MPI_CALL(MPI_Barrier, (net::MpiCommunicator::World()));
        auto const result = Comms().Rank() == 0 ? "root_result": "others";
        options = std::make_shared<CommandLine>(
          CommandLine{"hemelb", "-in", "large_cylinder.xml",
            "-i", "1", "-ss", "1111", "-out", result});
      }

      net::MpiCommunicator MPISpreadForcesTests::GetGraphComm(net::MpiCommunicator const &comm)
      {
        if(comm.Size() == 1)
        {
        }
        // setups a graph communicator that in-practice is all-to-all
        // Simpler than setting up something realistic
        std::vector<std::vector<int>> vertices;
        for(size_t i(0); i < comm.Size(); ++i)
        {
          vertices.push_back(std::vector<int>());
          for(size_t j(0); j < comm.Size(); ++j)
          {
            if(j != i)
            {
              vertices[i].push_back(j);
            }
          }
        }
        return comm.Graph(vertices);
      }

      std::vector<LatticePosition> MPISpreadForcesTests::GetPositions(
          geometry::LatticeData const & latDat,
          size_t mid, size_t edges, net::MpiCommunicator const &c, bool erase)
      {
        std::random_device rd;
        std::mt19937 g(rd());

        int const nMids = latDat.GetMidDomainSiteCount();
        int const nEdges = latDat.GetDomainEdgeCollisionCount(0);
        std::vector<LatticePosition> positions(c.Size() * (mid + edges));
        std::vector<int> shuf(nMids);
        std::iota(shuf.begin(), shuf.end(), 0);
        std::shuffle(shuf.begin(), shuf.end(), g);
        mid += std::max(0, static_cast<int>(edges) - nEdges);
        edges = std::min(edges, static_cast<size_t>(nEdges));
        for(size_t i(0); i < mid; ++i)
        {
          auto const site = latDat.GetSite(shuf[i]);
          positions[c.Rank() * (mid +  edges) + i] = site.GetGlobalSiteCoords();
        }
        shuf.resize(nEdges);
        std::iota(shuf.begin(), shuf.end(), 0);
        std::shuffle(shuf.begin(), shuf.end(), g);
        for(size_t i(0); i < edges; ++i)
        {
          auto const site = latDat.GetSite(nMids + shuf[i]);
          positions[c.Rank() * (mid +  edges) + i + mid] = site.GetGlobalSiteCoords();
        }

        auto const sendType = net::MpiDataType<LatticePosition>();
        HEMELB_MPI_CALL(
            MPI_Allgather,
            (MPI_IN_PLACE, mid + edges, sendType, positions.data(), mid + edges, sendType, c)
        );

        if(erase)
        {
          for(size_t i(0); i < positions.size() - mid - edges; ++i)
          {
            positions[i] = positions[i + mid + edges];
          }
          positions.resize(positions.size() - mid - edges);
        }

        return positions;
      }

      std::vector<CellContainer::value_type> MPISpreadForcesTests::GetCells(
          geometry::LatticeData const & latDat,
          size_t mid, size_t edges, net::MpiCommunicator const &c, bool erase)
      {
        auto const positions = GetPositions(latDat, mid, edges, c, erase);
        auto const n = mid + edges;
        std::vector<CellContainer::value_type> cells;
        for(auto i_first = positions.begin(); i_first != positions.end(); i_first += n)
        {
          cells.push_back(std::make_shared<DummyCell>(
                std::vector<LatticePosition>{i_first, i_first + n}, 1e0));
        }
        return cells;
      }

      template<class STENCIL> void MPISpreadForcesTests::testMidregion()
      {
        Check<STENCIL>(2, 0);
      }

      template<class STENCIL> void MPISpreadForcesTests::Check(size_t mid, size_t edges)
      {
        if(net::MpiCommunicator::World().Size() == 1)
        {
          return;
        }
        auto const world = net::MpiCommunicator::World();
        auto const color = world.Rank() == 0;
        auto const split = world.Split(color);
        auto master = CreateMasterSim<STENCIL>(split);
        auto &latDat = master->GetLatticeData();
        helpers::ZeroOutForces(latDat);

        // Figure out positions to use for cell nodes
        auto const cells = GetCells(latDat, mid, edges, world, true);
        auto const owned = split.Size() != 1 ?
          CellContainer{cells[split.Rank()]}:
          CellContainer{cells.begin(), cells.end()};
        auto const distributions = nodeDistributions(latDat, owned);

        OpenedSpreadForces mpi_spreader(GetGraphComm(split));
        mpi_spreader.PostMessageLength(distributions, owned);
        mpi_spreader.ComputeForces(owned);
        mpi_spreader.PostForcesAndNodes(distributions, owned);
        mpi_spreader.SpreadLocalForces<typename MasterSim<STENCIL>::Traits>(latDat, owned);
        mpi_spreader.SpreadNonLocalForces(latDat);

        std::vector<LatticeVector> indices;
        std::vector<LatticeForceVector> forces;
        if(color)
        {
          for(size_t i(0); i < latDat.GetLocalFluidSiteCount(); ++i)
          {
            auto const site = latDat.GetSite(i);
            if(site.GetForce().GetMagnitudeSquared() > 1e-8)
            {
              indices.push_back(site.GetGlobalSiteCoords());
              forces.push_back(site.GetForce());
            }
          }
        }
        int n = indices.size();
        world.Broadcast(n, 0);
        CPPUNIT_ASSERT(n > 0);
        indices.resize(n);
        world.Broadcast(indices, 0);

        if(not color)
        {
          for(auto coords: indices)
          {
            auto const id = latDat.GetProcIdFromGlobalCoords(coords);
            if(id == split.Rank())
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
        CPPUNIT_ASSERT_EQUAL(size_t(n), size_t(forces.size()));
        std::vector<LatticeForceVector> summed(indices.size());
        HEMELB_MPI_CALL(MPI_Allreduce,
          (forces.data(), summed.data(), forces.size() * 3,
           net::MpiDataType<double>(), MPI_SUM, world));
        for(auto const force: summed)
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
