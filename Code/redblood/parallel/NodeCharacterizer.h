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
#include "redblood/cell.h"
#include "redblood/stencil.h"
#include "redblood/interpolation.h"
#include "redblood/VelocityInterpolation.h"
#include "geometry/latticeData.h"
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
          : NodeCharacterizer(AssessMPIFunction<typename Traits<ARGS...>::Stencil>(latDat), cell)
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

          Process2NodesMap::mapped_type const & operator[](Process2NodesMap::key_type proc) const
          {
            auto result = affectedProcs.find(proc);
            assert(result != affectedProcs.end());
            return result->second;
          }

          //! Updates node characterization and return change in ownership
          void Reindex(AssessNodeRange const& assessNodeRange, MeshData::Vertices const &vertices);
          //! Reindex with normal mpi function
          template<class STENCIL>
          void Reindex(geometry::LatticeData const &latDat, std::shared_ptr<CellBase const> cell)
          {
            Reindex(AssessMPIFunction<STENCIL>(latDat), cell->GetVertices());
          }

          //! Consolidates result from another proc into an input array
          // void ReduceFrom(
          //     MeshData::Vertices const &consolidated,
          //     proc_t node, MeshData::Vertices const& incoming) const;

        protected:
          template<class STENCIL>
            static AssessNodeRange AssessMPIFunction(geometry::LatticeData const &latDat)
            {
              return std::bind(
                details::positionAffectsProcs<STENCIL>,
                latDat, std::placeholders::_1
              );
            }
          //! Processes affected by a given processor
          Process2NodesMap affectedProcs;
      };

      // //! An object for nodes owned by this process
      // class MeshOwner : public NodeCharacterizer
      // {
      //   public:
      //     MeshOwner(
      //         AssessNodeRange const& assessNodeRange,
      //         std::shared_ptr<CellBase> cell)
      //       : NodeCharacterizer(assessNodeRange, cell), cell(cell)
      //     {
      //     }
      //     //! Constructs object using a custom assessor function
      //     template<class ... ARGS>
      //     MeshOwner(
      //         geometry::LatticeData const &latDat,
      //         std::shared_ptr<CellBase> cell,
      //         redblood::stencil::types stencil=redblood::stencil::types::FOUR_POINT)
      //       : NodeCharacterizer(latDat, cell, Traits<ARGS...>()), cell(cell)
      //     {
      //     }
      //
      //     //! Apply functor to mid-domain nodes
      //     template<class APPLY> void ApplyToMidDomainNodes(APPLY apply, proc_t proc) const;
      //     //! Update mid-domain positions by computing velocities
      //     template <typename KERNEL> void UpdateLocalMidDomainPositions(
      //         geometry::LatticeData const &latDat,
      //         stencil::types stencil = stencil::types::FOUR_POINT) const;
      //
      //     //! Recomputes node locations and packs diffs for other procs
      //     void UpdateNodeCharacterization(
      //         AssessNodeRange const& assessNodeRange,
      //         proc_t const localRank,
      //         std::map<proc_t, util::Packer> &packer);
      //
      //   protected:
      //     using NodeCharacterizer::affectedProcs;
      //     //! Reference to vertices characterized by this object
      //     std::shared_ptr<CellBase> cell;
      // };

      // //! \brief An object for nodes owned by another processor
      // //! \details Internally, there are two kind of nodes: mid-domain nodes and boundary nodes.
      // //! The former do not interact with sites owned by other procs. The latter do. The interaction
      // //! range is determined by the size of the stencil for velocity interpolation and force
      // //! spreading. (Though a larger stencil could also be used).
      // class DistributedNodes
      // {
      //   public:
      //     //! type of the vector of nodes for each proc
      //     typedef std::vector<LatticePosition> Nodes;
      //     //! constructor
      //     DistributedNodes()
      //     {
      //     }
      //     //! Explicit construction from a map of proc -> vector of nodes
      //     DistributedNodes(std::map<proc_t, Nodes> && mid, std::map<proc_t, Nodes> && bound)
      //       : midDomain(mid), boundary(bound)
      //     {
      //     }
      //     //! destructor
      //     virtual ~DistributedNodes();
      //
      //     //! Apply functor to mid-domain nodes
      //     template<class KERNEL> void UpdateMidDomainPositions(
      //         geometry::LatticeData const &latDat,
      //         stencil::types stencil = stencil::types::FOUR_POINT);
      //     //! Send new node positions to owner
      //     void PackMidDomainNodes(util::Packer &packer, proc_t proc) const
      //     {
      //       auto const nodes = midDomain.find(proc);
      //       assert(nodes != midDomain.end());
      //       packer << nodes->second;
      //     }
      //     //! Send partial boundary velocity to owner
      //     // void ComputeAndPackBoundaryDomainNodes(util::Packer &packer, proc_t proc) const;
      //
      //   protected:
      //     //! MidDomain nodes and the processor they are owned by
      //     std::map<proc_t, Nodes> midDomain;
      //     //! Boundary nodes and the processor they are owned by
      //     std::map<proc_t, Nodes> boundary;
      // };

      // template<class APPLY> void MeshOwner::ApplyToMidDomainNodes(APPLY apply, proc_t proc) const
      // {
      //   auto i_procs = affectedProcs.cbegin();
      //   auto const i_procs_end = affectedProcs.cend();
      //   auto i_vertex = cell->GetVertices().begin();
      //   for(; i_procs != i_procs_end; ++i_procs, ++i_vertex)
      //   {
      //     if(i_procs->count(proc) == 1 and i_procs->size() == 1)
      //     {
      //       apply(*i_vertex);
      //     }
      //   }
      // }

      // template<class APPLY>
      //   void DistributedNodes::ApplyToMidDomainNodes(APPLY apply, proc_t proc) const
      // {
      //   for(; i_procs != i_procs_end; ++i_procs, ++i_vertex)
      //   {
      //     if(i_procs->count(proc) == 1 and i_procs->size() == 1)
      //     {
      //       apply(*i_vertex);
      //     }
      //   }
      // }

      // template<class KERNEL> void MeshOwner::UpdateLocalMidDomainPositions(
      //     geometry::LatticeData const &latDat, stencil::types stencil) const
      // {
      //   auto velocity = [&latDat, &stencil](LatticePosition &position)
      //   {
      //     position += interpolateVelocity<KERNEL>(latDat, position, stencil);
      //   };
      //   ApplyToMidDomainNodes(velocity, latDat.GetLocalRank());
      // }
      //
      // template<class KERNEL> void DistributedNodes::UpdateLocalMidDomainPositions(
      //     geometry::LatticeData const &latDat, stencil::types stencil)
      // {
      //   assert(midDomain.count(proc) == 1);
      //   for(auto & position: midDomain[proc])
      //   {
      //     position += interpolateVelocity<KERNEL>(latDat, position, stencil);
      //   }
      // }

      // template<class KERNEL> void DistributedNodes::ComputeAndPackBoundaryDomainNodes(
      //     util::Packer & packer, proc_t proc,
      //     geometry::LatticeData const &latDat, stencil::types stencil) const
      // {
      //   assert(boundary.count(proc) == 1);
      //   for(auto & position: boundary[proc])
      //   {
      //     packer << interpolateVelocity<KERNEL>(latDat, position, stencil);
      //   }
      // }

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
