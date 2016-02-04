//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURES_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURES_H

#include "redblood/Cell.h"
#include "redblood/parallel/CellParallelization.h"
#include "SimulationMaster.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      //! \brief gathers mid-domain and egde positions from all procs
      //! \details If there are insufficient number of edges, mid-domains are used instead.
      //! erase removes the components from the first process.
      std::vector<LatticePosition> GatherSpecialPositions(geometry::LatticeData const & latDat,
                                                          size_t mid, size_t edges,
                                                          net::MpiCommunicator const &c)
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
        for (size_t i(0); i < mid; ++i)
        {
          auto const site = latDat.GetSite(shuf[i]);
          positions[c.Rank() * (mid + edges) + i] = site.GetGlobalSiteCoords();
        }
        shuf.resize(nEdges);
        std::iota(shuf.begin(), shuf.end(), 0);
        std::shuffle(shuf.begin(), shuf.end(), g);
        for (size_t i(0); i < edges; ++i)
        {
          auto const site = latDat.GetSite(nMids + shuf[i]);
          positions[c.Rank() * (mid + edges) + i + mid] = site.GetGlobalSiteCoords();
        }

        auto const sendType = net::MpiDataType<LatticePosition>();
        HEMELB_MPI_CALL(MPI_Allgather,
                        (MPI_IN_PLACE, mid + edges, sendType, positions.data(), mid + edges, sendType, c));

        // Erase contribution from first proc (acting as serial oracle)
        for (size_t i(0); i < positions.size() - mid - edges; ++i)
        {
          positions[i] = positions[i + mid + edges];
        }
        positions.resize(positions.size() - mid - edges);

        return positions;
      }

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

      class DummyCell : public NodeCell
      {
        public:
          LatticeForce force;

          DummyCell(LatticePosition const&position, LatticeForce f = 0e0,
                    std::string const &templateName = "nope") :
              DummyCell(std::vector<LatticePosition> { position }, f, templateName)
          {
          }
          DummyCell(std::vector<LatticePosition> const &positions, LatticeForce f = 0e0,
                    std::string const &templateName = "nope") :
              NodeCell(positions, templateName), force(f)
          {
          }
          template<class ITER>
          DummyCell(ITER first, ITER last, LatticeForce f = 0e0, std::string const &templateName =
                        "nope") :
              NodeCell(first, last, templateName), force(f)
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
            return std::unique_ptr<DummyCell> { new DummyCell(*this) };
          }
      };

      //! Creates list of cells for each set of positions from each process
      std::vector<hemelb::redblood::CellContainer::value_type> CreateCellsFromSpecialPositions(
          geometry::LatticeData const & latDat, size_t mid, size_t edges,
          net::MpiCommunicator const &c, size_t nCells)
      {
        using hemelb::redblood::CellContainer;
        auto const positions = GatherSpecialPositions(latDat, mid * nCells, edges * nCells, c);
        std::vector<CellContainer::value_type> cells;
        for (auto i_first = positions.begin(); i_first != positions.end(); i_first += mid + edges)
        {
          cells.push_back(std::make_shared<DummyCell>(std::vector<LatticePosition> { i_first,
                                                                                     i_first + mid
                                                                                         + edges },
                                                      1e0));
        }
        return cells;
      }

      net::MpiCommunicator CreateDumbGraphComm(net::MpiCommunicator const &comm)
      {
        if (comm.Size() == 1)
        {
        }
        // setups a graph communicator that in-practice is all-to-all
        // Simpler than setting up something realistic
        std::vector<std::vector<int>> vertices;
        for (std::size_t i(0); i < std::size_t(comm.Size()); ++i)
        {
          vertices.push_back(std::vector<int>());
          for (std::size_t j(0); j < std::size_t(comm.Size()); ++j)
          {
            if (j != i)
            {
              vertices[i].push_back(j);
            }
          }
        }
        return comm.Graph(vertices);
      }

    }
  }
}

#endif  // ONCE
