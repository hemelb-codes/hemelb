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
  namespace
  {
    //! Iterates over vertices of a mesh and the nearby nodes of a grid
    //! The functor argument is called with the current vertex index, the
    //! global site index triplet, and the associated interpolation weight.
    template<class FUNCTOR, class STENCIL>
    void spreadForce2Grid(std::shared_ptr<CellBase const> cell, FUNCTOR functor)
    {
      typedef MeshData::Vertices::const_iterator const_iterator;
      // Spread them onto lattice
      const_iterator i_vertex = cell->GetVertices().begin();
      const_iterator const i_end = cell->GetVertices().end();

      for (size_t i(0); i_vertex != i_end; ++i_vertex, ++i)
      {
        InterpolationIterator<STENCIL> spreader = interpolationIterator<STENCIL>(*i_vertex);

        for (; spreader; ++spreader)
        {
          functor(i, *spreader, spreader.weight());
        }
      }
    }

    class SpreadForces
    {
      public:
        SpreadForces(std::vector<LatticePosition> const &forces, geometry::LatticeData &latticeData) :
            latticeData(latticeData), forces(forces)
        {
        }

        void operator()(size_t vertex, LatticeVector const &site, Dimensionless weight)
        {
          proc_t procid;
          site_t siteid;

          if (latticeData.GetContiguousSiteId(site, procid, siteid))
            geometry::Site < geometry::LatticeData
                > (latticeData.GetSite(site)).AddToForce(forces[vertex] * weight);
        }

      protected:
        geometry::LatticeData &latticeData;
        std::vector<LatticeForceVector> const &forces;
    };
  }
} // namespace details::anonymous
