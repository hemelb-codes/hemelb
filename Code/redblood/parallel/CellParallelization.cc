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
        auto const neighbors = cellCount.GetCommunicator().GetNeighbors();

        // Sets up information about number of cells to send
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


        cellCount.receive();
        totalNodeCount.receive();


        nodeCount.SetReceiveCounts(cellCount.GetReceiveBuffer());
        nodeCount.send();
      }

      void ExchangeCells::ReceiveCells()
      {
      }

      std::vector<size_t> ExchangeCells::GetNodesPerCells(
              CellParallelization::NodeDistributions const &owned,
              std::vector<int> neighbors) const
      {
        std::vector<size_t> result;
        for(auto const & dist: owned)
        {
          for(auto item: util::enumerate(neighbors))
          {
            auto const nVertices = dist.second.CountNodes(item.value);
            if(nVertices > 0)
            {
              result.push_back(nVertices);
            }
          }
        }
        return result;
      }

    } // parallel
  } // redblood
}  // hemelb
