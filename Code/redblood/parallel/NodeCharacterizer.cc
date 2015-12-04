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
#include "redblood/parallel/NodeCharacterizer.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      namespace details
      {
        template<class FUNCTION>
          NodeCharacterizer::Process2NodesMap meshMessenger(
              FUNCTION assessNodeRange, MeshData::Vertices const &vertices)
        {
          NodeCharacterizer::Process2NodesMap result;
          for(auto const vertex: util::enumerate(vertices))
          {
            for(auto const &process: assessNodeRange(vertex.value))
            {
              auto const i_proc = result.find(process);
              if(i_proc == result.end())
              {
                result[process] = {vertex.index};
              }
              else
              {
                result[process].insert(vertex.index);
              }
            }
          }
          return result;
        }
      }


      NodeCharacterizer::NodeCharacterizer(
              AssessNodeRange const& assessNodeRange, MeshData::Vertices const &vertices)
        : NodeCharacterizer(details::meshMessenger(assessNodeRange, vertices))
        {
        }

      NodeCharacterizer :: NodeCharacterizer(
          AssessNodeRange const & assessNodeRange, std::shared_ptr<CellBase const> cell)
        : NodeCharacterizer(assessNodeRange, cell->GetVertices())
        {
        }

      void NodeCharacterizer :: Reindex(
          AssessNodeRange const & assessor, MeshData::Vertices const & vertices)
      {
        affectedProcs = details::meshMessenger(assessor, vertices);
      }

      bool NodeCharacterizer::IsMidDomain(Index index) const
      {
        int found (0);
        for(auto const process: affectedProcs)
        {
          if(process.second.count(index) and ++found > 1)
          {
            return false;
          }
        }
        return true;
      }

      std::set<NodeCharacterizer::Process2NodesMap::key_type>
        NodeCharacterizer::AffectedProcs(Index index) const
      {
        std::set<NodeCharacterizer::Process2NodesMap::key_type> result;
        for(auto const process: affectedProcs)
        {
          if(process.second.count(index))
          {
            result.insert(process.first);
          }
        }
        return result;
      }

      void NodeCharacterizer::ReduceFrom(
          MeshData::Vertices &consolidated,
          Process2NodesMap::key_type node, MeshData::Vertices const& incoming) const
      {
        assert(affectedProcs.count(node) == 1);
        assert(affectedProcs.find(node)->second.size() == incoming.size());
        for(auto const && item: util::czip(affectedProcs.find(node)->second, incoming))
        {
          assert(std::get<0>(item) < consolidated.size());
          consolidated[std::get<0>(item)] += std::get<1>(item);
        }
      }

      void NodeCharacterizer::ReduceFrom(
          MeshData::Vertices &consolidated, MeshData::Vertices const& incoming) const
      {
        auto incoming_node = incoming.cbegin();
        for(auto const proc: affectedProcs)
        {
          for(auto const index: proc.second)
          {
            assert(incoming_node != incoming.end());
            consolidated[index] += *incoming_node;
            ++incoming_node;
          }
        }
      }

      void NodeCharacterizer::SpreadTo(
          std::vector<size_t> & sizes, MeshData::Vertices & outgoing,
          MeshData::Vertices const &vertices) const
      {
        for(auto const proc: affectedProcs)
        {
          sizes.push_back(proc.second.size());
          for(auto const index: proc.second)
          {
            assert(index < vertices.size());
            outgoing.push_back(vertices[index]);
          }
        }
      }
    } // parallel
  } // redblood
}  // hemelb
