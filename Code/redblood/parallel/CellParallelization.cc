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
#include <algorithm>
#include "util/Iterator.h"
#include "net/MpiError.h"
#include "redblood/parallel/CellParallelization.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      namespace
      {
        std::vector<int> GetUUIDCounts(std::vector<int> const &cellCounts)
        {
          std::vector<int> result = cellCounts;
          auto const uuidSize = [](int a)
          {
            return a * sizeof(boost::uuids::uuid);
          };
          std::transform(result.begin(), result.end(), result.begin(), uuidSize);
          return result;
        }
      }

      void ExchangeCells::PostCellMessageLength(CellParallelization::NodeDistributions const &owned)
      {
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
          CellParallelization::NodeDistributions const &owned, CellContainer const & cells)
      {

        // sets up nodeCount's send buffer
        SetupLocalSendBuffers(owned, cells);
        // nodeCount's receive buffer depends on the number of incoming cell from each neigbor
        cellCount.receive();
        nodeCount.SetReceiveCounts(cellCount.GetReceiveBuffer());
        cellScales.SetReceiveCounts(cellCount.GetReceiveBuffer());
        cellUUIDs.SetReceiveCounts(GetUUIDCounts(cellCount.GetReceiveBuffer()));

        nodeCount.send();
        cellScales.send();
        cellUUIDs.send();

        totalNodeCount.receive();
      }

      void ExchangeCells::ReceiveCells()
      {
      }

      void ExchangeCells::SetupLocalSendBuffers(
          CellParallelization::NodeDistributions const &owned, CellContainer const &cells)
      {
        // Sets up size of messages to send to neighbors
        nodeCount.SetSendCounts(cellCount.GetSendBuffer());
        cellScales.SetSendCounts(cellCount.GetSendBuffer());

        cellUUIDs.SetSendCounts(GetUUIDCounts(cellCount.GetSendBuffer()));

        for(auto neighbor: nodeCount.GetCommunicator().GetNeighbors())
        {
          int i(0);
          for(auto const & dist: owned)
          {
            auto const nVertices = dist.second.CountNodes(neighbor);
            if(nVertices > 0)
            {
              auto byUUID = [&dist](CellContainer::const_reference cell)
              {
                return cell->GetTag() == dist.first;
              };
              auto const i_cell = std::find_if(cells.begin(), cells.end(), byUUID);
              assert(i_cell != cells.end());
              AddToLocalSendBuffers(neighbor, i++, nVertices, *i_cell);
            }
          }
        }
      }

      void ExchangeCells::AddToLocalSendBuffers(
          int neighbor, int nth, int nVertices, CellContainer::const_reference cell)
      {
        nodeCount.SetSend(neighbor, nVertices, nth);
        cellScales.SetSend(neighbor, cell->GetScale(), nth);
        // boost::uuids::uuid are pod structures making up an array of unsigned chars.
        auto const & uuid = cell->GetTag();
        for(size_t j(0); j < sizeof(boost::uuids::uuid); ++j)
        {
          cellUUIDs.SetSend(neighbor, *(uuid.begin() + j), nth * sizeof(boost::uuids::uuid) + j);
        }
      }
    } // parallel
  } // redblood
}  // hemelb
