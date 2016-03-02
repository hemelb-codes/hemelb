//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_REDBLOODE_PARALLELIZATION_SPREAD_FORCES_H
#define HEMELB_REDBLOODE_PARALLELIZATION_SPREAD_FORCES_H

#include <map>
#include <vector>

#include <boost/uuid/uuid.hpp>

#include "redblood/parallel/NodeCharacterizer.h"
#include "redblood/parallel/CellParallelization.h"
#include "redblood/Cell.h"
#include "redblood/GridAndCell.h"

#include "net/MpiCommunicator.h"
#include "net/INeighborAllToAll.h"
#include "net/INeighborAllToAllV.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      class SpreadForces
      {
        public:
          //! Forces for a single cell
          typedef std::vector<LatticeForceVector> Forces;
          //! \brief Forces for owned cells
          typedef std::map<boost::uuids::uuid, Forces> CellForces;
          //! Type of the object holding distributions
          typedef CellParallelization::NodeDistributions NodeDistributions;
          //! Container holding cells lent by other processes
          typedef CellParallelization::LentCells LentCells;

          SpreadForces(net::MpiCommunicator const &graphComm) :
              sendNodeCount(graphComm), sendPositions(graphComm), sendForces(graphComm)
          {
          }

          //! \brief Computes and posts length of message when sending forces to other procs
          //| \note This function must be called before PostForcesAndNode
          //! \param[in] distributions: Node distributions of the cells owned by this process
          //! \param[in] owned: Cells currently owned by this process
          void PostMessageLength(NodeDistributions const& distributions,
                                 CellContainer const &owned);
          //! Computes and caches forces
          //! \note This function must be called prior to PostForcesAndNodes and SpreadLocal. It
          //! can be called before or after PostMessageLength.
          //! \return Sum of energies over all cells
          virtual LatticeEnergy ComputeForces(CellContainer const &owned);
          //! \brief Post non-local forces
          //! \param[in] distributions tells us for each proc the list of nodes it requires
          //! \param[in] cells a container of cells owned and managed by this process
          void PostForcesAndNodes(NodeDistributions const &distributions,
                                  CellContainer const &owned);
          //! \brief Spreads local forces
          //! \tparam TRAITS defines TRAITS::Stencil needed to actually do the spreading.
          template<class TRAITS = Traits<>>
          void SpreadLocalForces(geometry::LatticeData & latticeData,
                                 CellContainer const &owned) const;
          //! \brief Receive and spread forces from other procs
          template<class TRAITS = Traits<>>
          void SpreadNonLocalForces(geometry::LatticeData & latticeData);

        protected:
          //! Sends total number of shared nodes
          net::INeighborAllToAll<int> sendNodeCount;
          //! Sends positions for the nodes
          net::INeighborAllToAllV<LatticePosition> sendPositions;
          //! Sends forces
          net::INeighborAllToAllV<LatticeForceVector> sendForces;
          //! Holds forces for each cell
          CellForces cellForces;

          //! Helper function to iterate over cells that do need sending
          template<class FUNCTOR>
          void IterateOverMessageCells(NodeDistributions const& distributions,
                                       CellContainer const &owned, FUNCTOR functor);
      };

      template<class FUNCTOR>
      void SpreadForces::IterateOverMessageCells(NodeDistributions const& distributions,
                                                 CellContainer const &owned, FUNCTOR functor)
      {
        auto const neighbors = sendNodeCount.GetCommunicator().GetNeighbors();
        // Count the number of vertices and cells
        for (auto const & cell : owned)
        {
          assert(distributions.count(cell->GetTag()) == 1);
          auto const dist = distributions.find(cell->GetTag())->second;
          for (auto const neighbor : util::enumerate(neighbors))
          {
            auto const nVertices = dist.CountNodes(neighbor.value);
            if (nVertices > 0)
            {
              functor(neighbor.index, neighbor.value, dist[neighbor.value], cell);
            }
          }
        }
      }

      template<class TRAITS>
      void SpreadForces::SpreadLocalForces(geometry::LatticeData & latticeData,
                                           CellContainer const &owned) const
      {
        namespace hrd = hemelb::redblood::details;
        typedef typename TRAITS::Stencil Stencil;
        for (auto const cell : owned)
        {
          assert(cellForces.count(cell->GetTag()) == 1);
          auto const& forces = cellForces.find(cell->GetTag())->second;
          hrd::spreadForce2Grid<hrd::SpreadForces, Stencil>(cell,
                                                            hrd::SpreadForces(forces, latticeData));
        }
      }

      template<class TRAITS>
      void SpreadForces::SpreadNonLocalForces(geometry::LatticeData &latticeData)
      {
        namespace hrd = hemelb::redblood::details;
        typedef typename TRAITS::Stencil Stencil;

        sendPositions.receive();
        sendForces.receive();
        assert(sendPositions.GetReceiveBuffer().size() == sendForces.GetReceiveBuffer().size());
#       ifndef NDEBUG
        if(sendPositions.GetReceiveBuffer().size() > 0)
        {
          log::Logger::Log<log::Debug, log::OnePerCore>(
              "Number of forces spread from other procs: %i",
              sendPositions.GetReceiveBuffer().size()
          );
        }
#       endif
        hrd::spreadForce2Grid<hrd::SpreadForces, Stencil>(sendPositions.GetReceiveBuffer(),
                                                          hrd::SpreadForces(sendForces.GetReceiveBuffer(),
                                                                            latticeData));
      }
    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
