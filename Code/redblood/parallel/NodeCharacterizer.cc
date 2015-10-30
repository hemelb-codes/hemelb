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
        void UpdateNodeCharacterization(
            std::set<proc_t> const &newProcs, std::set<proc_t> const &oldProcs,
            proc_t localRank, LatticePosition const &,
            std::map<proc_t, util::Packer> &)
        {
          bool const
            wasMidDomain = oldProcs.size() == 1,
            isMidDomain = newProcs.size() == 1,
            wasLocal = oldProcs.count(localRank),
            isLocal = newProcs.count(localRank);
          if(wasMidDomain == isMidDomain and wasLocal != isLocal)
          {
          }
        }


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

      // void MeshOwner::UpdateNodeCharacterization(
      //     AssessNodeRange const& assessNodeRange,
      //     proc_t const localRank,
      //     std::map<proc_t, util::Packer> &packers)
      // {
      //   MeshData::Vertices :: const_iterator vertex = cell->GetVertices().begin();
      //   MeshData::Vertices :: const_iterator const vertex_end = cell->GetVertices().end();
      //   ProcessorSets::iterator procs = affectedProcs.begin();
      //   assert(affectedProcs.size() == cell->GetVertices().size());
      //   for(; vertex != vertex_end; ++vertex, ++procs)
      //   {
      //     auto const newProcs = assessNodeRange(*vertex);
      //     parallel::details::UpdateNodeCharacterization(
      //         newProcs, *procs, localRank, *vertex, packers);
      //     *procs = std::move(newProcs);
      //   }
      // }
    } // parallel
  } // redblood
}  // hemelb
