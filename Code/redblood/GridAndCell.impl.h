// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_GRIDANDCELL_IMPL_H
#define HEMELB_REDBLOOD_GRIDANDCELL_IMPL_H

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

#endif
