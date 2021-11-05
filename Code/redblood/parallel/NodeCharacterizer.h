// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_REDBLOOD_PARALLEL_NODECHARACTERIZER_H
#define HEMELB_REDBLOOD_PARALLEL_NODECHARACTERIZER_H

#include <set>
#include <memory>
#include <functional>

#include "redblood/Cell.h"
#include "redblood/stencil.h"
#include "redblood/Interpolation.h"
#include "redblood/VelocityInterpolation.h"
#include "geometry/LatticeData.h"
#include "Traits.h"
#include "redblood/parallel/GraphBasedCommunication.h"

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
        //! \param[in] globalCoordsToProcMap will tell us which site belongs to which proc
        //! \param[in] iterator a  stencil iterator going over affected lattice points
        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(GlobalCoordsToProcMap const &globalCoordsToProcMap,
                                              InterpolationIterator<STENCIL> &&iterator);
        //! Set of procs affected by this position
        //! \param[in] globalCoordsToProcMap will tell us which site belongs to which proc
        //! \param[in] position for which to figure out affected processes
        //! \param[in] stencil giving interaction range
        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(GlobalCoordsToProcMap const &globalCoordsToProcMap,
                                              LatticePosition const &position);

        //! Simplest MPI range assessor loops over nodes + interpolation stencil
        template<class STENCIL>
        static AssessNodeRange AssessMPIFunction(GlobalCoordsToProcMap const &globalCoordsToProcMap)
        {
          return [&globalCoordsToProcMap](LatticePosition const &pos)
          {
            auto const& affectedProcs = details::positionAffectsProcs<STENCIL>(
                globalCoordsToProcMap, InterpolationIterator<STENCIL>(pos));
            // #652 No mesh vertex should be under the influence of no rank in a valid simulation.
            if (affectedProcs.size() == 0)
            {
              std::stringstream message;
              message << "Mesh vertex at " << pos << " is not affected by any flow subdomain.";
              log::Logger::Log<log::Error, log::OnePerCore>(message.str());

              throw std::exception();
            }
            return affectedProcs;
          };
        }
      } /* details */

      //! NodeCharacterizer computes which processes are affected by each node
      //! in a RBC mesh, keeps a map between process rank and node indices, and
      //! provides methods to query it. Mesh node and vertex refer to the same
      //! concept here.
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
          NodeCharacterizer(Process2NodesMap const &affectedProcs) :
              affectedProcs(affectedProcs)
          {
          }
          //! Constructs object from prior knowledge of how processors are affected
          NodeCharacterizer(Process2NodesMap &&affectedProcs) :
              affectedProcs(affectedProcs)
          {
          }
          //! Constructs object using stencil as range object
          //! In practice, this will be the main constructor. Others are there for testing and
          //! convenience.
          NodeCharacterizer(AssessNodeRange const &assessNodeRange,
                            MeshData::Vertices const &vertices);
          //! Constructs object using stencil as range object
          //! In practice, this will be the main constructor. Others are there for testing and
          //! convenience.
          NodeCharacterizer(AssessNodeRange const &assessNodeRange,
                            std::shared_ptr<CellBase const> cell);

          //! Constructs object using a custom assessor function
          template<class ... ARGS>
          NodeCharacterizer(GlobalCoordsToProcMap const &globalCoordsToProcMap,
                            std::shared_ptr<CellBase const> cell, Traits<ARGS...> const&) :
                  NodeCharacterizer(details::AssessMPIFunction<typename Traits<ARGS...>::Stencil>(globalCoordsToProcMap),
                                    cell)
          {
          }

          //! Whether the node affects only one processor
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
          Process2NodesMap::mapped_type::difference_type CountNodes(
              Process2NodesMap::key_type proc) const
          {
            auto result = affectedProcs.find(proc);
            return result == affectedProcs.end() ?
              0 :
              result->second.size();
          }

          //! Indices of nodes that affect more than one proc
          Process2NodesMap::value_type::second_type BoundaryIndices() const;
          //! Indices of the processors affected by any node in the mesh
          std::set<Process2NodesMap::key_type> AffectedProcs() const;

          //! The process affected by the largest number of nodes
          //! \return The rank of the process affected by the largest number of
          //! mesh nodes or -1 if this could not be determined
          Process2NodesMap::key_type DominantAffectedProc() const;

          //! Updates node characterization and return change in ownership
          void Reindex(AssessNodeRange const& assessNodeRange, MeshData::Vertices const &vertices);
          //! Reindex with normal mpi function
          template<class STENCIL>
          void Reindex(GlobalCoordsToProcMap const &globalCoordsToProcMap, std::shared_ptr<CellBase const> cell)
          {
            Reindex(details::AssessMPIFunction<STENCIL>(globalCoordsToProcMap), cell->GetVertices());
          }

          //! Consolidates result from another proc into an input array
          void ReduceFrom(MeshData::Vertices &consolidated, proc_t node,
                          MeshData::Vertices const& incoming) const;
          //! \brief Reduce results from all nodes
          //! \details Assumes that the incomming array is aranged such that it contains first all
          //! the nodes from the lowest rank, then the next lowest, ...
          void ReduceFrom(MeshData::Vertices &consolidated,
                          MeshData::Vertices const& incoming) const;
          //! \brief Spread data to different node
          //! \details Creates an array with all the nodes to send over to different processor. The
          //! output is tailored to fit MPI's gatherv.
          void SpreadTo(std::vector<size_t> & sizes, MeshData::Vertices &outgoing,
                        MeshData::Vertices const& vertices) const;
          //! \brief Spread data to different node
          //! \details Creates an array with all the nodes to send over to different processor. The
          //! output is tailored to fit MPI's gatherv.
          std::pair<std::vector<size_t>, MeshData::Vertices> SpreadTo(
              MeshData::Vertices const& vertices) const
          {
            std::pair<std::vector<size_t>, MeshData::Vertices> result;
            SpreadTo(result.first, result.second, vertices);
            return result;
          }

        protected:
          //! Nodes affected by a given processor
          Process2NodesMap affectedProcs;
      };

      namespace details
      {
        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(GlobalCoordsToProcMap const &globalCoordsToProcMap,
                                              InterpolationIterator<STENCIL> &&iterator)
        {
          std::set<proc_t> result;
          for (; iterator.IsValid(); ++iterator)
          {
            auto const& id = globalCoordsToProcMap.find(*iterator);
            //! @todo #668 Some unit tests throw the warning below. Requires further investigation.
            if(id == globalCoordsToProcMap.end())
            {
              std::stringstream message;
              message << "No owner recorded for lattice site " << *iterator << ". Not a problem iff outside flow domain.";
              log::Logger::Log<log::Debug, log::OnePerCore>(message.str());
            }
            else
            {
              result.insert(id->second);
            }
          }

          return result;
        }

        template<class STENCIL>
        std::set<proc_t> positionAffectsProcs(GlobalCoordsToProcMap const &globalCoordsToProcMap,
                                              LatticePosition const &position)
        {
          return positionAffectsProcs(globalCoordsToProcMap, interpolationIterator<STENCIL>(position));
        }
      } /* details */

    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
