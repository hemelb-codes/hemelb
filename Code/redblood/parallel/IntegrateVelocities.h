// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_REDBLOOD_PARALLEL_INTEGRATEVELOCITIES_H
#define HEMELB_REDBLOOD_PARALLEL_INTEGRATEVELOCITIES_H

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
      //! \brief Integrates velocities across processes
      //! \details Uses the fact that the velocities are additive. We update the positions from the
      //! local knowledge in one step, and then from the non-local knowledge in another.
      class IntegrateVelocities
      {
        public:
          //! Forces for a single cell
          typedef std::vector<LatticeVelocity> Velocities;

          IntegrateVelocities(net::MpiCommunicator const &graphComm) :
              sendNodeCount(graphComm), sendVelocities(graphComm)
          {
          }

          //! \brief Computes and posts length of message when sending velocities to other procs
          //! \param[in] lent: Cells lent to this proc
          //! \note Must be called before PostVelocities
          void PostMessageLength(LentCells const &lent);
          //! Computes and caches velocities
          //! \tparam TRAITS defines TRAITS::Stencil and TRAITS::Kernel needed to actually integrate
          //! the velocities.
          //! \brief After this call, the node-positions are in an indeterminate state: the nodes
          //! have been update with portions of the velocities that this proc knows about. However,
          //! the velocities from other procs are not integrated until UpdatePositionsNonLocal.
          //! \note Can be called at anytime.
          template<class TRAITS = Traits<>>
          void ComputeLocalVelocitiesAndUpdatePositions(geometry::LatticeData const &latDat,
                                                        CellContainer &owned);
          //! \brief Post non-local velocities
          //! \param[in] distributions tells us for each proc the list of nodes it requires
          //! \param[in] cells a container of cells owned and managed by this process
          //! \note Must be called after PostMessageLength and before UpdatePositionsLocal
          template<class TRAITS = Traits<>>
          void PostVelocities(geometry::LatticeData const &latDat, LentCells const &lent);
          //! \brief Gathers velocities from other procs and upates positions
          //! \note Must be called after PostVelocities
          void UpdatePositionsNonLocal(NodeDistributions const& distributions,
                                       CellContainer &owned);

        protected:
          //! Sends total number of shared nodes
          net::INeighborAllToAll<int> sendNodeCount;
          //! Sends positions for the nodes
          net::INeighborAllToAllV<LatticeVelocity> sendVelocities;
      };

      template<class TRAITS>
      void IntegrateVelocities::ComputeLocalVelocitiesAndUpdatePositions(
          geometry::LatticeData const &latticeData, CellContainer &owned)
      {
        typedef typename TRAITS::Kernel Kernel;
        typedef typename TRAITS::Stencil Stencil;
        std::vector<LatticeVelocity> velocities;
        for (auto const &cell : owned)
        {
          velocities.resize(cell->GetNumberOfNodes());
          std::fill(velocities.begin(), velocities.end(), LatticeVelocity { 0, 0, 0 });
          velocitiesOnMesh<Kernel, Stencil>(cell, latticeData, velocities);
          *cell += velocities;
        }
      }

      template<class TRAITS>
      void IntegrateVelocities::PostVelocities(geometry::LatticeData const &latDat,
                                               LentCells const &lent)
      {
        typedef typename TRAITS::Kernel Kernel;
        typedef typename TRAITS::Stencil Stencil;
        sendVelocities.SetSendCounts(sendNodeCount.GetSendBuffer());
        std::map<int, int> offsets;
        for (auto const neighbor : sendVelocities.GetCommunicator().GetNeighbors())
        {
          offsets[neighbor] = 0;
        }
        std::vector<LatticeVelocity> velocities;
        for (auto const lentCells : lent)
        {
          auto const neighbor = lentCells.first;
          for (auto const &cell : lentCells.second)
          {
            assert(offsets.count(neighbor) == 1);
            auto const offset = offsets[neighbor];
            offsets[neighbor] += cell->GetNumberOfNodes();
            velocities.resize(cell->GetNumberOfNodes());
            std::fill(velocities.begin(), velocities.end(), LatticeVelocity { 0, 0, 0 });
            velocitiesOnMesh<Kernel, Stencil>(cell, latDat, velocities);
            for (auto const &vertex : util::enumerate(velocities))
            {
              sendVelocities.SetSend(neighbor, vertex.value, offset + vertex.index);
            }
          }
        }

        sendNodeCount.receive();
        sendVelocities.SetReceiveCounts(sendNodeCount.GetReceiveBuffer());
        sendVelocities.send();
      }

    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
