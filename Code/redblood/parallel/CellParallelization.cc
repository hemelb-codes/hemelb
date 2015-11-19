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
#include <functional>
#include "util/Iterator.h"
#include "net/MpiError.h"
#include "redblood/parallel/CellParallelization.h"
#include "redblood/VertexBag.h"

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
        std::vector<int> GetTotalNodeCounts(std::vector<int> const &cellCounts)
        {
          std::vector<int> result = cellCounts;
          auto const uuidSize = [](int a)
          {
            return a * 3;
          };
          std::transform(result.begin(), result.end(), result.begin(), uuidSize);
          return result;
        }

        template<class ITER> ITER cellFromTag(boost::uuids::uuid const &tag, ITER first, ITER last)
        {
          auto byUUID = [&tag](CellContainer::const_reference cell)
          {
            return cell->GetTag() == tag;
          };
          return std::find_if(first, last, byUUID);
        }

        CellContainer::const_reference cellFromTag(
            boost::uuids::uuid const &tag, CellContainer const &cells)
        {
          auto const i_cell = cellFromTag(tag, cells.begin(), cells.end());
          assert(i_cell != cells.end());
          return *i_cell;
        }

        std::string getTemplateName(size_t index, std::vector<char> const &buffer)
        {
          // loop over n occurences
          auto first = buffer.begin();
          for(; index > 0 and first != buffer.end(); --index, ++first)
          {
            first = std::find(first, buffer.end(), '\0');
          }
          return std::string(first, std::find(first, buffer.end(), '\0'));
        }

        boost::uuids::uuid getTag(size_t index, std::vector<unsigned char> const &buffer)
        {
          auto const n = sizeof(boost::uuids::uuid);
          boost::uuids::uuid result;
          std::copy(buffer.begin() + index * n, buffer.begin() + (index + 1) * n, result.begin());
          return result;
        }
      }

      void ExchangeCells::PostCellMessageLength(
          NodeDistributions const &distributions, CellContainer const &owned)
      {
        auto const neighbors = cellCount.GetCommunicator().GetNeighbors();
        cellCount.GetSendBuffer().resize(neighbors.size());
        totalNodeCount.GetSendBuffer().resize(neighbors.size());
        nameLengths.GetSendBuffer().resize(neighbors.size());
        std::fill(cellCount.GetSendBuffer().begin(), cellCount.GetSendBuffer().end(), 0);
        std::fill(totalNodeCount.GetSendBuffer().begin(), totalNodeCount.GetSendBuffer().end(), 0);
        std::fill(nameLengths.GetSendBuffer().begin(), nameLengths.GetSendBuffer().end(), 0);
        // Count the number of vertices and cells
        for(auto const & dist: distributions)
        {
          auto const cell = cellFromTag(dist.first, owned);
          for(auto item: util::enumerate(neighbors))
          {
            auto const nVertices = dist.second.CountNodes(item.value);
            if(nVertices > 0)
            {
              ++cellCount.GetSendBuffer()[item.index];
              totalNodeCount.GetSendBuffer()[item.index] += nVertices;
              // Include null termination that will seperate one name from another
              nameLengths.GetSendBuffer()[item.index] += cell->GetTemplateName().size() + 1;
            }
          }
        }

        // Then send message
        cellCount.send();
        totalNodeCount.send();
        nameLengths.send();
      }

      void ExchangeCells::PostCells(
          NodeDistributions const &distributions, CellContainer const & owned,
          Ownership const &ownership)
      {
        auto const rankMap = net::MpiCommunicator::World().RankMap(cellCount.GetCommunicator());
        std::map<boost::uuids::uuid, proc_t> id;
        auto const getID = [&rankMap, &ownership](CellContainer::const_reference cell)
        {
          auto const worldIndex = ownership(cell);
          return std::make_pair(cell->GetTag(), rankMap.find(worldIndex)->second);
        };
        std::transform(owned.begin(), owned.end(), std::inserter(id, id.begin()), getID);
        PostCells(distributions, owned, id);
      }

      void ExchangeCells::PostCells(
          NodeDistributions const &distributions, CellContainer const & owned,
          std::map<boost::uuids::uuid, proc_t> const &ownership)
      {
        // sets up nodeCount's send buffer
        SetupLocalSendBuffers(distributions, owned, ownership);

        // nodeCount's receive buffer depends on the number of incoming cell from each neigbor
        cellCount.receive();
        nodeCount.SetReceiveCounts(cellCount.GetReceiveBuffer());
        cellScales.SetReceiveCounts(cellCount.GetReceiveBuffer());
        ownerIDs.SetReceiveCounts(cellCount.GetReceiveBuffer());
        cellUUIDs.SetReceiveCounts(GetUUIDCounts(cellCount.GetReceiveBuffer()));

        // Now can send messages that depend only on number of cells
        nodeCount.send();
        cellScales.send();
        ownerIDs.send();
        cellUUIDs.send();

        // Need total number of nodes to prepare buffers for positions
        totalNodeCount.receive();
        nodePositions.SetReceiveCounts(GetTotalNodeCounts(totalNodeCount.GetReceiveBuffer()));
        nodePositions.send();

        // Need template name lengths to prepare buffers for template names
        nameLengths.receive();
        templateNames.SetReceiveCounts(nameLengths.GetReceiveBuffer());
        templateNames.send();
      }

      void ExchangeCells::ReceiveCells(
              CellContainer &owned, std::map<proc_t, CellContainer> &lent,
              std::shared_ptr<TemplateCellContainer const> const &templateCells)
      {
        nodeCount.receive();
        cellScales.receive();
        ownerIDs.receive();
        cellUUIDs.receive();
        nodePositions.receive();
        templateNames.receive();

        for(size_t i(0); i < nodeCount.GetReceiveBuffer().size(); ++i)
        {
          auto const ownerID = ownerIDs.GetReceiveBuffer()[i];
          bool const isOwned = ownerID == ownerIDs.GetCommunicator().Rank();
          auto cell = isOwned ? RecreateOwnedCell(i, templateCells): RecreateLentCell(i);
          auto & container = isOwned ? owned: lent[ownerID];
          container.insert(cell);
        }
      }

      CellContainer::value_type ExchangeCells::RecreateLentCell(size_t index)
      {
        // first extract tag and template name
        auto const tag = getTag(index, cellUUIDs.GetReceiveBuffer());
        auto const templateName = getTemplateName(index, templateNames.GetReceiveBuffer());

        // create cell
        auto result = std::make_shared<VertexBag>(tag, templateName);
        // add nodes
        auto i_node = nodePositions.GetReceiveBuffer().cbegin();
        for(size_t i(0); i < index; ++index, i_node += nodeCount.GetReceiveBuffer()[i] * 3);
        for(size_t j(0); j < nodeCount.GetReceiveBuffer()[index]; ++j)
        {
          result->addVertex({*(i_node++), *(i_node++), *(i_node++)});
        }
        // and set the scale
        result->SetScale(cellScales.GetReceiveBuffer()[index]);
        return result;
      }

      CellContainer::value_type ExchangeCells::RecreateOwnedCell(
              size_t index, std::shared_ptr<TemplateCellContainer const> const &templateCells)
      {
        auto const templateName = getTemplateName(index, templateNames.GetReceiveBuffer());
        auto result = templateCells->find(templateName)->second->clone();
        result->SetTag(getTag(index, cellUUIDs.GetReceiveBuffer()));
        result->SetScale(cellScales.GetReceiveBuffer()[index]);
        auto i_node = nodePositions.GetReceiveBuffer().cbegin();
        for(size_t i(0); i < index; ++index, i_node += nodeCount.GetReceiveBuffer()[i] * 3);
        auto const Nnodes = nodeCount.GetReceiveBuffer()[index];
        assert(Nnodes == result->GetNumberOfNodes());
        for(size_t j(0); j < Nnodes; ++j)
        {
          result->GetVertices()[j].x = *(i_node++);
          result->GetVertices()[j].y = *(i_node++);
          result->GetVertices()[j].z = *(i_node++);
        }
        return std::move(result);
      }

      void ExchangeCells::SetupLocalSendBuffers(
          NodeDistributions const &distributions, CellContainer const &owned,
          std::map<boost::uuids::uuid, proc_t> const &ownership)
      {
        // Sets up size of messages to send to neighbors
        nodeCount.SetSendCounts(cellCount.GetSendBuffer());
        cellScales.SetSendCounts(cellCount.GetSendBuffer());
        ownerIDs.SetSendCounts(cellCount.GetSendBuffer());
        cellUUIDs.SetSendCounts(GetUUIDCounts(cellCount.GetSendBuffer()));
        nodePositions.SetSendCounts(GetTotalNodeCounts(totalNodeCount.GetSendBuffer()));
        templateNames.SetSendCounts(nameLengths.GetSendBuffer());
        templateNames.fillSend('\0');

        auto const thisRank = nodeCount.GetCommunicator().Rank();
        for(auto neighbor: nodeCount.GetCommunicator().GetNeighbors())
        {
          int i(0);
          for(auto const & dist: distributions)
          {
            auto const nVertices = dist.second.CountNodes(neighbor);
            if(nVertices > 0)
            {
              auto const cell = cellFromTag(dist.first, owned);
              auto const owner = ownership.find(dist.first)->second;
              owner == thisRank ?
                AddToLocalSendBuffers(neighbor, i++, nVertices, owner, cell, dist.second[neighbor]):
                AddToLocalSendBuffers(neighbor, i++, owner, cell);
            }
          }
        }
      }

      void ExchangeCells::AddToLocalSendBuffersAllButNodes(
          int neighbor, int nth, proc_t ownerID, CellContainer::const_reference cell)
      {
        // For now, onwership is left with sending process
        cellScales.SetSend(neighbor, cell->GetScale(), nth);
        ownerIDs.SetSend(neighbor, ownerID, nth);
        // boost::uuids::uuid are pod structures making up an array of unsigned chars.
        auto const & uuid = cell->GetTag();
        for(size_t j(0); j < sizeof(boost::uuids::uuid); ++j)
        {
          cellUUIDs.SetSend(neighbor, *(uuid.begin() + j), nth * sizeof(boost::uuids::uuid) + j);
        }
        // Copy template mesh name into message buffer
        // First, finds where out where we are in the stream.
        auto name_offset(0);
        for(int i(0); i < nth; ++i)
        {
          while(templateNames.GetSend(neighbor, ++name_offset) != '\0');
        }
        for(auto const item: util::enumerate(cell->GetTemplateName()))
        {
          templateNames.SetSend(neighbor, item.value, name_offset + item.index);
        }
      }

      void ExchangeCells::AddToLocalSendBuffers(
          int neighbor, int nth, proc_t ownerID, CellContainer::const_reference cell)
      {
        AddToLocalSendBuffersAllButNodes(neighbor, nth, ownerID, cell);
        nodeCount.SetSend(neighbor, cell->GetNumberOfNodes(), nth);
        size_t node_offset(0);
        for(int j(0); j < nth; ++j, node_offset += 3 * nodeCount.GetSend(neighbor, j));
        for(auto const item: util::enumerate(cell->GetVertices()))
        {
          nodePositions.SetSend(neighbor, item.value.x, node_offset + item.index * 3);
          nodePositions.SetSend(neighbor, item.value.y, node_offset + item.index * 3 + 1);
          nodePositions.SetSend(neighbor, item.value.z, node_offset + item.index * 3 + 2);
        }
      }

      void ExchangeCells::AddToLocalSendBuffers(
          int neighbor, int nth, int nVertices, proc_t ownerID,
          CellContainer::const_reference cell,
          NodeCharacterizer::Process2NodesMap::mapped_type const & indices)
      {
        AddToLocalSendBuffersAllButNodes(neighbor, nth, ownerID, cell);
        // For now, onwership is left with sending process
        nodeCount.SetSend(neighbor, nVertices, nth);
        size_t node_offset(0);
        for(int j(0); j < nth; ++j, node_offset += 3 * nodeCount.GetSend(neighbor, j));
        for(auto const item: util::enumerate(indices))
        {
          auto const& node = cell->GetVertices()[item.value];
          nodePositions.SetSend(neighbor, node.x, node_offset + item.index * 3);
          nodePositions.SetSend(neighbor, node.y, node_offset + item.index * 3 + 1);
          nodePositions.SetSend(neighbor, node.z, node_offset + item.index * 3 + 2);
        }
      }
    } // parallel
  } // redblood
}  // hemelb
