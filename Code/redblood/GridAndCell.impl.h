//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

namespace details
{
  //! Iterates over vertices of a mesh and the nearby nodes of a grid
  //! The functor argument is called with the current vertex index, the
  //! global site index triplet, and the associated interpolation weight.
  template<class FUNCTOR, class STENCIL>
  void spreadForce2Grid(std::vector<LatticePosition> const &positions, FUNCTOR functor)
  {
    for (auto const item : util::enumerate(positions))
    {
      InterpolationIterator<STENCIL> spreader { item.value };

      for (; spreader; ++spreader)
      {
        functor(item.index, *spreader, spreader.weight());
      }
    }
  }

  //! Iterates over vertices of a mesh and the nearby nodes of a grid
  //! The functor argument is called with the current vertex index, the
  //! global site index triplet, and the associated interpolation weight.
  template<class FUNCTOR, class STENCIL>
  void spreadForce2Grid(std::shared_ptr<CellBase const> cell, FUNCTOR functor)
  {
    spreadForce2Grid<FUNCTOR, STENCIL>(cell->GetVertices(), functor);
  }

  class SpreadForces
  {
    public:
      SpreadForces(std::vector<LatticePosition> const &forces, geometry::LatticeData &latticeData) :
          latticeData(latticeData), i_force(forces.begin())
      {
      }

      void operator()(size_t vertex, LatticeVector const &site, Dimensionless weight)
      {
        proc_t procid;
        site_t siteid;

        if (latticeData.GetContiguousSiteId(site, procid, siteid))
        {
          assert (procid == latticeData.GetCommunicator().Rank());
          auto siteOb = latticeData.GetSite(site);
          siteOb.AddToForce(* (i_force + vertex) * weight);
        }
      }

    protected:
      geometry::LatticeData &latticeData;
      std::vector<LatticeForceVector>::const_iterator const i_force;
  };

} // namespace details::anonymous
