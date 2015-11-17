//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <set>
#include <numeric>
#include "util/Iterator.h"
#include "net/MpiError.h"
#include "redblood/parallel/CellParallelization.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      void ExchangeCells::PostCellMessageLength(CellParallelization::NodeDistributions const &owned)
      {
        if(step % 3 != 0)
        {
          throw Exception() << "Out-of-order Exchange cell step called\n";
        }
        ++step;
        auto const neighbors = cellCount.GetCommunicator().GetNeighbors();
        cellCount.GetSendBuffer().resize(neighbors.size());
        totalNodeCount.GetSendBuffer().resize(neighbors.size());
        std::fill(cellCount.GetSendBuffer().begin(), cellCount.GetSendBuffer().end(), 0);
        std::fill(totalNodeCount.GetSendBuffer().begin(), totalNodeCount.GetSendBuffer().end(), 0);
        // Count the number of vertices and cells
        for(auto const & dist: owned)
        {
          for(auto item: util::enumerate(neighbors))
          {
            auto const nVertices = dist.second.CountNodes(item.value);
            if(nVertices > 0)
            {
              ++cellCount.GetSendBuffer()[item.index];
              totalNodeCount.GetSendBuffer()[item.index] += nVertices;
            }
          }
        }

        // Then send message
        cellCount.send();
        totalNodeCount.send();
      }

      void ExchangeCells::PostCells(
          CellParallelization::NodeDistributions const &owned, CellContainer const &)
      {

        // sets up nodeCount's send buffer
        SetupLocalNodeCount(owned);
        // nodeCount's receive buffer depends on the number of incoming cell from each neigbor
        cellCount.receive();
        nodeCount.SetReceiveCounts(cellCount.GetReceiveBuffer());
        nodeCount.send();

        // Now set up array to receive nodes
        totalNodeCount.receive();
      }

      void ExchangeCells::ReceiveCells()
      {
      }

      void ExchangeCells::SetupLocalNodeCount(CellParallelization::NodeDistributions const &owned)
      {
        auto const neighbors = cellCount.GetCommunicator().GetNeighbors();
        nodeCount.SetSendCounts(cellCount.GetSendBuffer());
        for(auto neighbor: neighbors)
        {
          int i(0);
          for(auto const & dist: owned)
          {
            auto const nVertices = dist.second.CountNodes(neighbor);
            if(nVertices > 0)
            {
              nodeCount.SetSend(neighbor, nVertices, i++);
            }
          }
        }
      }

    } // parallel
  } // redblood
}  // hemelb
