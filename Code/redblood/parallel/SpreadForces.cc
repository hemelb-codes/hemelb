//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/parallel/SpreadForces.h"
#include "util/Iterator.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      void SpreadForces::PostMessageLength(NodeDistributions const& distributions,
                                           CellContainer const &owned)
      {
        typedef NodeCharacterizer::Process2NodesMap::mapped_type NodeIndices;
        auto countNodesToSend =
            [this](int index, int, NodeIndices const & indices, CellContainer::value_type const&)
            {
              assert(int(sendNodeCount.GetSendBuffer().size()) > index);
              sendNodeCount.GetSendBuffer()[index] += indices.size();
            };

        sendNodeCount.GetSendBuffer().resize(sendNodeCount.GetCommunicator().GetNeighborsCount());
        std::fill(sendNodeCount.GetSendBuffer().begin(), sendNodeCount.GetSendBuffer().end(), 0);

        IterateOverMessageCells(distributions, owned, countNodesToSend);

        // Then send message
        sendNodeCount.send();
      }

      // Computes and caches forces
      LatticeEnergy SpreadForces::ComputeForces(CellContainer const &owned)
      {
        // Clear forces and make sure there are enough of them
        cellForces.clear();
        // compute forces for each
        LatticeEnergy energy(0);
        for (auto const &cell : owned)
        {
          // Create the map entry and allocate memory
          cellForces[cell->GetTag()].resize(cell->GetNumberOfNodes());

          energy += cell->Energy(cellForces[cell->GetTag()]);
        }
        return energy;
      }

      void SpreadForces::PostForcesAndNodes(NodeDistributions const &distributions,
                                            CellContainer const &owned)
      {
        sendPositions.SetSendCounts(sendNodeCount.GetSendBuffer());
        sendForces.SetSendCounts(sendNodeCount.GetSendBuffer());
        std::map<int, int> offsets;
        for (auto const neighbor : sendPositions.GetCommunicator().GetNeighbors())
        {
          offsets[neighbor] = 0;
        }
        typedef NodeCharacterizer::Process2NodesMap::mapped_type NodeIndices;
        auto provisionSendBuffer = [=](
            int, int neighbor, NodeIndices const & indices,
            CellContainer::value_type const &cell) mutable
        {
          auto const offset = offsets[neighbor];
          offsets[neighbor] += indices.size();
          for(auto const item: util::enumerate(indices))
          {
            auto const& node = cell->GetVertices()[item.value];
            sendPositions.SetSend(neighbor, node, offset + item.index);

            auto const& force = cellForces[cell->GetTag()][item.value];
            sendForces.SetSend(neighbor, force, offset + item.index);
          }
        };

        IterateOverMessageCells(distributions, owned, provisionSendBuffer);

        sendNodeCount.receive();
        sendPositions.SetReceiveCounts(sendNodeCount.GetReceiveBuffer());
        sendForces.SetReceiveCounts(sendNodeCount.GetReceiveBuffer());
        sendPositions.send();
        sendForces.send();
      }
    } // parallel
  } // redblood
} // hemelb
