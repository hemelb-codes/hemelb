// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
        NodeCharacterizer::Process2NodesMap meshMessenger(FUNCTION assessNodeRange,
                                                          MeshData::Vertices const &vertices)
        {
          NodeCharacterizer::Process2NodesMap result;
          for (auto const &vertex : util::enumerate(vertices))
          {
            for (auto const &process : assessNodeRange(vertex.value))
            {
              result[process].insert(vertex.index);
            }
          }
          return result;
        }
      }

      NodeCharacterizer::NodeCharacterizer(AssessNodeRange const& assessNodeRange,
                                           MeshData::Vertices const &vertices) :
          NodeCharacterizer(details::meshMessenger(assessNodeRange, vertices))
      {
      }

      NodeCharacterizer::NodeCharacterizer(AssessNodeRange const & assessNodeRange,
                                           std::shared_ptr<CellBase const> cell) :
          NodeCharacterizer(assessNodeRange, cell->GetVertices())
      {
      }

      void NodeCharacterizer::Reindex(AssessNodeRange const & assessor,
                                      MeshData::Vertices const & vertices)
      {
        affectedProcs = details::meshMessenger(std::cref(assessor), vertices);
      }

      bool NodeCharacterizer::IsMidDomain(Index index) const
      {
        int found(0);
        for (auto const &process : affectedProcs)
        {
          if (process.second.count(index) and ++found > 1)
          {
            return false;
          }
        }
        return true;
      }

      NodeCharacterizer::Process2NodesMap::value_type::second_type NodeCharacterizer::BoundaryIndices() const {
        if(affectedProcs.size() == 1)
        {
          return decltype(BoundaryIndices())();
        }
        auto i_first = affectedProcs.cbegin();
        auto current = i_first->second;
        for(++i_first; i_first != affectedProcs.cend(); ++i_first) {
          std::set<MeshData::Vertices::size_type> next;
          std::set_intersection(current.begin(), current.end(),
                                i_first->second.begin(), i_first->second.end(),
                                std::inserter(next, next.begin()));
          std::swap(next, current);
        }
        return current;
      }

      std::set<NodeCharacterizer::Process2NodesMap::key_type>
      NodeCharacterizer::AffectedProcs() const
      {
        std::set<Process2NodesMap::key_type> result;
        for(auto const &process: affectedProcs)
        {
          if(process.second.size() > 0)
          {
            result.insert(process.first);
          }
        }
        return result;
      }

      std::set<NodeCharacterizer::Process2NodesMap::key_type> NodeCharacterizer::AffectedProcs(
          Index index) const
      {
        std::set<NodeCharacterizer::Process2NodesMap::key_type> result;
        for (auto const &process : affectedProcs)
        {
          if (process.second.count(index))
          {
            result.insert(process.first);
          }
        }
        return result;
      }

      NodeCharacterizer::Process2NodesMap::key_type NodeCharacterizer::DominantAffectedProc() const
      {
        typedef NodeCharacterizer::Process2NodesMap::const_reference Arg;
        auto const max = std::max_element(affectedProcs.begin(), affectedProcs.end(), [](Arg a, Arg b) { return a.second.size() < b.second.size(); });
        if (affectedProcs.empty() or max->second.size() == 0)
        {
          return -1;
        }
        return max->first;
      }

      void NodeCharacterizer::ReduceFrom(MeshData::Vertices &consolidated,
                                         Process2NodesMap::key_type node,
                                         MeshData::Vertices const& incoming) const
      {
        assert(affectedProcs.count(node) == 1);
        assert(affectedProcs.find(node)->second.size() == incoming.size());
        for (auto const && item : util::czip(affectedProcs.find(node)->second, incoming))
        {
          assert(std::get<0>(item) < consolidated.size());
          consolidated[std::get<0>(item)] += std::get<1>(item);
        }
      }

      void NodeCharacterizer::ReduceFrom(MeshData::Vertices &consolidated,
                                         MeshData::Vertices const& incoming) const
      {
        auto incoming_node = incoming.cbegin();
        for (auto const proc : affectedProcs)
        {
          for (auto const index : proc.second)
          {
            assert(incoming_node != incoming.end());
            consolidated[index] += *incoming_node;
            ++incoming_node;
          }
        }
      }

      void NodeCharacterizer::SpreadTo(std::vector<size_t> & sizes, MeshData::Vertices & outgoing,
                                       MeshData::Vertices const &vertices) const
      {
        for (auto const proc : affectedProcs)
        {
          sizes.push_back(proc.second.size());
          for (auto const index : proc.second)
          {
            assert(index < vertices.size());
            outgoing.push_back(vertices[index]);
          }
        }
      }
    } // parallel
  } // redblood
} // hemelb
