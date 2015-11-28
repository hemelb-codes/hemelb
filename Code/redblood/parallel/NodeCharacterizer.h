//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_REDBLOODE_PARALLELIZATION_NODE_CHARACTERIZATION
#define HEMELB_REDBLOODE_PARALLELIZATION_NODE_CHARACTERIZATION

#include <set>
#include <memory>
#include <functional>

#include "util/Packer.h"
#include "redblood/Cell.h"
#include "redblood/stencil.h"
#include "redblood/Interpolation.h"
#include "redblood/VelocityInterpolation.h"
#include "geometry/LatticeData.h"
#include "Traits.h"


namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      namespace details
      {
        //! Maps process indices to a sequence of node indices
        typedef std::map<proc_t, std::set<MeshData::Vertices::size_type>> Process2NodesMap;
        //! Functions returning the set of affected procs for a given node
        typedef std::function<std::set<proc_t>(LatticePosition const&)> AssessNodeRange;
        //! Set of procs affected by this position
        //! \param[in] latDat will tell us which site belongs to which proc
        //! \param[in] iterator a  stencil iterator going over affected lattice points
        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(
            geometry::LatticeData const &latDat, InterpolationIterator<STENCIL> &&iterator);
        //! Set of procs affected by this position
        //! \param[in] latDat will tell us which site belongs to which proc
        //! \param[in] position for which to figure out affected processes
        //! \param[in] stencil giving interaction range
        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(
            geometry::LatticeData const &latDat, LatticePosition const &position);

        //! Simplest MPI range assessor loops over nodes + interpolation stencil
        template<class STENCIL>
          static AssessNodeRange AssessMPIFunction(geometry::LatticeData const &latDat)
          {
            return [&latDat](LatticePosition const &pos)
            {
              return details::positionAffectsProcs<STENCIL>(
                  latDat, InterpolationIterator<STENCIL>(pos));
            };
          }
      } /* details */

      class NodeCharacterizer
      {
        public:
          //! set of processors affected by a node
          typedef details::Process2NodesMap Process2NodesMap;
          //! Indices of the nodes
          typedef Process2NodesMap::value_type::second_type::size_type Index;
          //! A function to assess which processes a node may affect
          typedef details::AssessNodeRange AssessNodeRange;

          //! Constructs object from prior knowledge of how processors are affected
          NodeCharacterizer(Process2NodesMap const &affectedProcs)
            : affectedProcs(affectedProcs)
          {
          }
          //! Constructs object from prior knowledge of how processors are affected
          NodeCharacterizer(Process2NodesMap &&affectedProcs) : affectedProcs(affectedProcs)
          {
          }
          //! Constructs object using stencil as range object
          //! In practice, this will be the main constructor. Others are there for testing and
          //! convenience.
          NodeCharacterizer(
              AssessNodeRange const &assessNodeRange, MeshData::Vertices const &vertices);
          //! Constructs object using stencil as range object
          //! In practice, this will be the main constructor. Others are there for testing and
          //! convenience.
          NodeCharacterizer(
              AssessNodeRange const &assessNodeRange, std::shared_ptr<CellBase const> cell);

          //! Constructs object using a custom assessor function
          template<class ... ARGS>
          NodeCharacterizer(
              geometry::LatticeData const &latDat,
              std::shared_ptr<CellBase const> cell,
              Traits<ARGS...> const& )
          : NodeCharacterizer(
              details::AssessMPIFunction<typename Traits<ARGS...>::Stencil>(latDat), cell)
          {
          }

          //! Whether the node affects more than one processor
          bool IsMidDomain(Index index) const;
          //! Whether the node affects more than one processor
          bool IsBoundary(Index index) const
          {
            return not IsMidDomain(index);
          }
          //! \brief The processors affected by a given node
          std::set<Process2NodesMap::key_type> AffectedProcs(Index i) const;

          //! The nodes affected by a given proc
          Process2NodesMap::mapped_type const & operator[](Process2NodesMap::key_type proc) const
          {
            auto result = affectedProcs.find(proc);
            assert(result != affectedProcs.end());
            return result->second;
          }
          //! The number of nodes affected by a given proc
          Process2NodesMap::mapped_type::difference_type
          CountNodes(Process2NodesMap::key_type proc) const
          {
            auto result = affectedProcs.find(proc);
            return result == affectedProcs.end() ? 0: result->second.size();
          }

          //! Updates node characterization and return change in ownership
          void Reindex(AssessNodeRange const& assessNodeRange, MeshData::Vertices const &vertices);
          //! Reindex with normal mpi function
          template<class STENCIL>
          void Reindex(geometry::LatticeData const &latDat, std::shared_ptr<CellBase const> cell)
          {
            Reindex(details::AssessMPIFunction<STENCIL>(latDat), cell->GetVertices());
          }

          //! Consolidates result from another proc into an input array
          void ReduceFrom(
              MeshData::Vertices &consolidated,
              proc_t node, MeshData::Vertices const& incoming) const;
          //! \brief Reduce results from all nodes
          //! \details Assumes that the incomming array is aranged such that it contains first all
          //! the nodes from the lowest rank, then the next lowest, ...
          void ReduceFrom(
              MeshData::Vertices &consolidated, MeshData::Vertices const& incoming) const;
          //! \brief Spread data to different node
          //! \details Creates an array with all the nodes to send over to different processor. The
          //! output is tailored to fit MPI's gatherv.
          void SpreadTo(
              std::vector<size_t> & sizes, MeshData::Vertices &outgoing,
              MeshData::Vertices const& vertices) const;
          //! \brief Spread data to different node
          //! \details Creates an array with all the nodes to send over to different processor. The
          //! output is tailored to fit MPI's gatherv.
          std::pair<std::vector<size_t>, MeshData::Vertices>
          SpreadTo(MeshData::Vertices const& vertices) const
          {
            std::pair<std::vector<size_t>, MeshData::Vertices> result;
            SpreadTo(result.first, result.second, vertices);
            return result;
          }


        protected:
          //! Processes affected by a given processor
          Process2NodesMap affectedProcs;
      };

      namespace details
      {
        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(
            geometry::LatticeData const &latDat, InterpolationIterator<STENCIL> &&iterator)
        {
          std::set<proc_t> result;
          for(; iterator.IsValid(); ++iterator)
          {
            proc_t procid;
            site_t siteid;
            if (latDat.GetContiguousSiteId(*iterator, procid, siteid))
            {
              result.insert(procid);
            }
          }
          return result;
        }

        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(
            geometry::LatticeData const &latDat,
            LatticePosition const &position)
        {
          return positionAffectsProcs(latDat, interpolationIterator<STENCIL>(position));
        }
      } /* details */

    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
