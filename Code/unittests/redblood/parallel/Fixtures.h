// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURES_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURES_H

#include "redblood/Cell.h"
#include "SimulationMaster.h"
#include <random>

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      //! \brief gathers mid-domain and egde positions from all procs
      //! \details If there are insufficient number of edges, mid-domains are used instead.
      //! erase removes the components from the first process.
      std::vector<LatticeVector> GatherSpecialPositions(geometry::LatticeData const & latDat,
                                                          size_t mid, size_t edges,
                                                          net::MpiCommunicator const &c)
      {
        std::random_device rd;
        std::mt19937 g(rd());

        int const nMids = latDat.GetMidDomainCollisionCount(0);
        int const nEdges = latDat.GetDomainEdgeCollisionCount(0);
        std::vector<LatticeVector> positions(c.Size() * (mid + edges));
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
          auto const site = latDat.GetSite(latDat.GetMidDomainSiteCount() + shuf[i]);
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

          std::shared_ptr<hemelb::configuration::SimConfig> GetSimConfig()
          {
            return this->simConfig;
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
            std::fill(f.begin(), f.end(), LatticeForceVector{force});
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
        return comm.Graph(hemelb::redblood::parallel::ComputeProcessorNeighbourhood(comm));
      }
    }
  }
}

#endif  // ONCE
