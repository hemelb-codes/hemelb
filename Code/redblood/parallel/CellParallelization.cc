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

        std::string getTemplateName(size_t index, std::vector<char> const &buffer)
        {
          // loop over n occurences
          auto first = buffer.begin();
          for (; index > 0 and first != buffer.end(); --index, ++first)
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

        std::map<boost::uuids::uuid, proc_t> getOwnership(net::MpiCommunicator const &graphComm,
                                                          net::MpiCommunicator const &world,
                                                          CellContainer const &owned,
                                                          ExchangeCells::Ownership const &ownership)
        {
          auto const rankMap = world.RankMap(graphComm);
          std::map<boost::uuids::uuid, proc_t> result;
          for (auto const &cell : owned)
          {
            result[cell->GetTag()] = rankMap.find(ownership(cell))->second;
          }
          return result;
        }
      }

      void ExchangeCells::PostCellMessageLength(NodeDistributions const& distributions,
                                                CellContainer const &owned,
                                                Ownership const & ownership)
      {
        PostCellMessageLength(distributions,
                              owned,
                              getOwnership(nodeCount.GetCommunicator(), simComm, owned, ownership));
      }

      void ExchangeCells::PostCellMessageLength(
          NodeDistributions const &distributions, CellContainer const &owned,
          std::map<boost::uuids::uuid, proc_t> const &ownership)
      {
        auto const neighbors = cellCount.GetCommunicator().GetNeighbors();
        cellCount.GetSendBuffer().resize(neighbors.size());
        totalNodeCount.GetSendBuffer().resize(neighbors.size());
        nameLengths.GetSendBuffer().resize(neighbors.size());
        std::fill(cellCount.GetSendBuffer().begin(), cellCount.GetSendBuffer().end(), 0);
        std::fill(totalNodeCount.GetSendBuffer().begin(), totalNodeCount.GetSendBuffer().end(), 0);
        std::fill(nameLengths.GetSendBuffer().begin(), nameLengths.GetSendBuffer().end(), 0);
        // Count the number of vertices and cells
        for (auto const & cell : owned)
        {
          auto const dist = distributions.find(cell->GetTag())->second;
          assert(distributions.count(cell->GetTag()) == 1);
          auto const newOwner = ownership.find(cell->GetTag())->second;
          for (auto item : util::enumerate(neighbors))
          {
            auto const nVertices = newOwner == item.value ?
              cell->GetNumberOfNodes() :
              dist.CountNodes(item.value);
            if (nVertices > 0)
            {
              assert(cellCount.GetSendBuffer().size() > item.index);
              assert(totalNodeCount.GetSendBuffer().size() > item.index);
              assert(nameLengths.GetSendBuffer().size() > item.index);
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

      void ExchangeCells::PostCells(NodeDistributions const &distributions,
                                    CellContainer const & owned, Ownership const &ownership)
      {
        PostCells(distributions,
                  owned,
                  getOwnership(nodeCount.GetCommunicator(), simComm, owned, ownership));
      }

      void ExchangeCells::PostCells(NodeDistributions const &distributions,
                                    CellContainer const & owned,
                                    std::map<boost::uuids::uuid, proc_t> const &ownership)
      {
        // sets up all main messages and the disowned cells (lent back to this process)
        disowned.clear();
        formelyOwned.clear();
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
        nodePositions.SetReceiveCounts(totalNodeCount.GetReceiveBuffer());
        nodePositions.send();

        // Need template name lengths to prepare buffers for template names
        nameLengths.receive();
        templateNames.SetReceiveCounts(nameLengths.GetReceiveBuffer());
        templateNames.send();
      }

      ExchangeCells::ChangedCells ExchangeCells::ReceiveCells(
          TemplateCellContainer const &templateCells)
      {
        nodeCount.receive();
        cellScales.receive();
        ownerIDs.receive();
        cellUUIDs.receive();
        nodePositions.receive();
        templateNames.receive();

        ChangedCells result;
        std::get<1>(result) = disowned;

        for (size_t i(0); i < nodeCount.GetReceiveBuffer().size(); ++i)
        {
          auto const ownerID = ownerIDs.GetReceiveBuffer()[i];
          bool const isOwned = ownerID == ownerIDs.GetCommunicator().Rank();
          auto cell = isOwned ?
            RecreateOwnedCell(i, templateCells) :
            RecreateLentCell(i);
          auto & container = isOwned ?
            std::get<0>(result) :
            std::get<2>(result)[ownerID];
          container.insert(cell);
        }

        // adds formely owned cells to lent cells
        for (auto const & item : formelyOwned)
        {
          for (auto const & cell : item.second)
          {
            std::get<2>(result)[item.first].insert(cell);
          }
        }
        return result;
      }

      CellContainer::value_type ExchangeCells::RecreateLentCell(size_t index)
      {
        // first extract tag and template name
        auto const tag = getTag(index, cellUUIDs.GetReceiveBuffer());
        auto const templateName = getTemplateName(index, templateNames.GetReceiveBuffer());

        // create cell
        auto result = std::make_shared<VertexBag>(tag, templateName);
        // add nodes
        auto const offset = std::accumulate(nodeCount.GetReceiveBuffer().cbegin(),
                                            nodeCount.GetReceiveBuffer().cbegin() + index,
                                            0);
        auto i_node = nodePositions.GetReceiveBuffer().cbegin() + offset;
        auto const Nnodes = nodeCount.GetReceiveBuffer()[index];
        for (size_t j(0); j < Nnodes; ++j)
        {
          result->addVertex(* (i_node++));
        }
        // and set the scale
        result->SetScale(cellScales.GetReceiveBuffer()[index]);
        return result;
      }

      CellContainer::value_type ExchangeCells::RecreateOwnedCell(
          size_t index, TemplateCellContainer const &templateCells)
      {
        auto const templateName = getTemplateName(index, templateNames.GetReceiveBuffer());
        assert(templateCells.count(templateName) == 1);
        auto result = templateCells.find(templateName)->second->clone();
        result->SetTag(getTag(index, cellUUIDs.GetReceiveBuffer()));
        assert(cellScales.GetReceiveBuffer().size() > index);
        result->SetScale(cellScales.GetReceiveBuffer()[index]);
        auto const offset = std::accumulate(nodeCount.GetReceiveBuffer().cbegin(),
                                            nodeCount.GetReceiveBuffer().cbegin() + index,
                                            0);
        auto i_node = nodePositions.GetReceiveBuffer().cbegin() + offset;
        auto const Nnodes = nodeCount.GetReceiveBuffer()[index];
        assert(site_t(Nnodes) == result->GetNumberOfNodes());
        std::copy(i_node, i_node + Nnodes, result->GetVertices().begin());
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
        nodePositions.SetSendCounts(totalNodeCount.GetSendBuffer());
        templateNames.SetSendCounts(nameLengths.GetSendBuffer());
        templateNames.fillSend('\0');

        auto const neighbors = nodeCount.GetCommunicator().GetNeighbors();
        auto const thisRank = nodeCount.GetCommunicator().Rank();
        std::map<int, int> nth;
        for (auto const &neighbor : neighbors)
        {
          nth[neighbor] = 0;
        }
        for (auto const &cell : owned)
        {
          auto const newOwner = ownership.find(cell->GetTag())->second;
          assert(distributions.count(cell->GetTag()) == 1);
          auto const & distribution = distributions.find(cell->GetTag())->second;
          for (auto const neighbor : neighbors)
          {
            auto const nVertices = distribution.CountNodes(newOwner == neighbor ?
              thisRank :
              neighbor);
            if (newOwner == neighbor and nVertices > 0)
            {
              AddDisownedToLocalSendBuffers(neighbor,
                                            nth[neighbor]++,
                                            cell,
                                            distribution[thisRank]);
            }
            else if (newOwner == neighbor)
            {
              AddDisownedToLocalSendBuffers(neighbor, nth[neighbor]++, cell);
            }
            else if (nVertices > 0)
            {
              AddOwnedToLocalSendBuffers(neighbor,
                                         nth[neighbor]++,
                                         nVertices,
                                         newOwner,
                                         cell,
                                         distribution[neighbor]);
            }
          }
        }
      }

      void ExchangeCells::AddToLocalSendBuffersAllButNodes(int neighbor, int nth, proc_t ownerID,
                                                           CellContainer::const_reference cell)
      {
        // For now, onwership is left with sending process
        cellScales.SetSend(neighbor, cell->GetScale(), nth);
        ownerIDs.SetSend(neighbor, ownerID, nth);
        // boost::uuids::uuid are pod structures making up an array of unsigned chars.
        auto const & uuid = cell->GetTag();
        for (size_t j(0); j < sizeof(boost::uuids::uuid); ++j)
        {
          cellUUIDs.SetSend(neighbor, * (uuid.begin() + j), nth * sizeof(boost::uuids::uuid) + j);
        }
        // Copy template mesh name into message buffer
        // First, finds where out where we are in the stream.
        auto name_offset(0);
        for (int i(0); i < nth; ++i)
        {
          while (templateNames.GetSend(neighbor, name_offset++) != '\0')
            ;
        }
        assert(nth == 0 or templateNames.GetSend(neighbor, name_offset - 1) == '\0');
        for (auto const item : util::enumerate(cell->GetTemplateName()))
        {
          templateNames.SetSend(neighbor, item.value, name_offset + item.index);
        }
      }

      void ExchangeCells::AddDisownedToLocalSendBuffers(int neighbor, int nth,
                                                        CellContainer::const_reference cell)
      {
        AddToLocalSendBuffersAllButNodes(neighbor, nth, neighbor, cell);
        nodeCount.SetSend(neighbor, cell->GetNumberOfNodes(), nth);
        size_t node_offset(0);
        for (int j(0); j < nth; ++j)
        {
          node_offset += nodeCount.GetSend(neighbor, j);
        }
        for (auto const item : util::enumerate(cell->GetVertices()))
        {
          nodePositions.SetSend(neighbor, item.value, node_offset + item.index);
        }
        // adds to disowned cells and to formely owned cells
        disowned.insert(cell);
      }

      void ExchangeCells::AddDisownedToLocalSendBuffers(
          int neighbor, int nth, CellContainer::const_reference cell,
          NodeCharacterizer::Process2NodesMap::mapped_type const & indices)
      {
        AddDisownedToLocalSendBuffers(neighbor, nth, cell);
        // create vertex bag if any nodes are lent to this object
        if (indices.size() > 0)
        {
          auto lentCell = std::make_shared<VertexBag>(cell->GetTag(), cell->GetTemplateName());
          lentCell->SetScale(cell->GetScale());
          // add nodes
          for (auto const index : indices)
          {
            lentCell->addVertex(cell->GetVertices()[index]);
          }
          if (formelyOwned.count(neighbor) == 1)
          {
            formelyOwned[neighbor] = CellContainer { lentCell };
          }
          else
          {
            formelyOwned[neighbor].emplace(std::move(lentCell));
          }
        }
      }

      void ExchangeCells::AddOwnedToLocalSendBuffers(
          int neighbor, int nth, int nVertices, proc_t ownerID, CellContainer::const_reference cell,
          NodeCharacterizer::Process2NodesMap::mapped_type const & indices)
      {
        AddToLocalSendBuffersAllButNodes(neighbor, nth, ownerID, cell);
        // For now, onwership is left with sending process
        nodeCount.SetSend(neighbor, nVertices, nth);
        size_t node_offset(0);
        for (int j(0); j < nth; ++j)
        {
          node_offset += nodeCount.GetSend(neighbor, j);
        }
        for (auto const item : util::enumerate(indices))
        {
          auto const& node = cell->GetVertices()[item.value];
          nodePositions.SetSend(neighbor, node, node_offset + item.index);
        }
      }

      void ExchangeCells::Update(CellContainer &owned, const ChangedCells &changes)
      {
        for (auto const& cell : std::get<0>(changes))
        {
          owned.insert(cell);
        }
        for (auto const& cell : std::get<1>(changes))
        {
          owned.erase(cell);
        }
      }

      void ExchangeCells::Update(NodeDistributions &distributions, const ChangedCells &changes,
                                 NodeCharacterizer::AssessNodeRange const &assessor)
      {
        for (auto const& cell : std::get<1>(changes))
        {
          distributions.erase(cell->GetTag());
        }
        for (auto const& cell : std::get<0>(changes))
        {
          if (distributions.count(cell->GetTag()) == 0)
          {
            distributions.emplace(std::piecewise_construct,
                                  std::forward_as_tuple(cell->GetTag()),
                                  std::forward_as_tuple(assessor, cell));
          }
        }
      }
    } // parallel
  } // redblood
} // hemelb
