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

    // Loops over nodes, computes velocity and does something
    // struct + member makes up for lack of partial function
    // specialization in c++ pre 11
    template<class KERNEL> struct VelocityNodeLoop
    {
      VelocityNodeLoop(
        stencil::types stencil,
        CellBase const &cell,
        geometry::LatticeData const &latDat
      ) : stencil(stencil), cell(cell), latticeData(latDat) {}
      // Loop and does something
      template<class FUNCTOR> void loop(FUNCTOR apply)
      {
        typedef MeshData::Vertices::const_iterator const_iterator;
        const_iterator i_current = cell.GetVertices().begin();
        const_iterator const i_end = cell.GetVertices().end();

        for(; i_current != i_end; ++i_current)
        {
          PhysicalVelocity const velocity
            = interpolateVelocity<KERNEL>(latticeData, *i_current, stencil);
          apply(velocity);
        }
      }

      stencil::types const stencil;
      CellBase const &cell;
      geometry::LatticeData const &latticeData;
    };

    //! Updates an assignable iterator of some kind
    template<class ITERATOR> struct TransformIterator
    {
      ITERATOR iterator;
      TransformIterator(ITERATOR iterator) : iterator(iterator) {}
      void operator()(typename ITERATOR::value_type const & value)
      {
        *(iterator++) = value;
      }
    };

    //! Updates an assignable iterator of some kind
    template<class ITERATOR>
    TransformIterator<ITERATOR> transform_iterator(ITERATOR iterator)
    {
      return TransformIterator<ITERATOR>(iterator);
    }

    //! Iterates over vertices of a mesh and the nearby nodes of a grid
    //! The functor argument is called with the current vertex index, the
    //! global site index triplet, and the associated interpolation weight.
    template<class FUNCTOR> void spreadForce2Grid(
      CellBase const &cell,
      FUNCTOR functor,
      stencil::types stencil
    )
    {
      typedef MeshData::Vertices::const_iterator const_iterator;
      // Spread them onto lattice
      const_iterator i_vertex = cell.GetVertices().begin();
      const_iterator const i_end = cell.GetVertices().end();

      for(size_t i(0); i_vertex != i_end; ++i_vertex, ++i)
      {
        InterpolationIterator spreader
          = interpolationIterator(*i_vertex, stencil);

        for(; spreader; ++spreader)
        {
          functor(i, *spreader, spreader.weight());
        }
      }
    }

    class SpreadForces
    {
      public:
        SpreadForces(
          std::vector<LatticePosition> const &forces,
          geometry::LatticeData &latticeData
        ) : latticeData(latticeData), forces(forces) {}

        void operator()(
          size_t vertex, LatticeVector const &site,
          Dimensionless weight)
        {
          proc_t procid;
          site_t siteid;

          if(latticeData.GetContiguousSiteId(site, procid, siteid))
            geometry::Site<geometry::LatticeData>(latticeData.GetSite(site))
            .AddToForce(forces[vertex] * weight);
        }

      protected:
        geometry::LatticeData &latticeData;
        std::vector<LatticeForceVector> const & forces;
    };

    template<class LATTICE>
    class SpreadForcesAndWallForces : public SpreadForces
    {
      public:
        SpreadForcesAndWallForces(
          Cell const &cell,
          std::vector<LatticePosition> const &forces,
          geometry::LatticeData &latticeData
        ) : SpreadForces(forces, latticeData), cell(cell) {}
        void operator()(
          size_t vertexIn, LatticeVector const &siteIn,
          Dimensionless weight)
        {
          proc_t procid;
          site_t siteid;

          if(not latticeData.GetContiguousSiteId(siteIn, procid, siteid))
          {
            return;
          }

          geometry::Site<geometry::LatticeData> site(latticeData.GetSite(siteid));
          site.AddToForce(forces[vertexIn] * weight);
          LatticePosition const vertex(cell.GetVertices()[vertexIn]);

          for(size_t i(1); i < LATTICE::NUMVECTORS; ++i)
          {
            PhysicalDistance const distance = site.GetWallDistance<LATTICE>(i);

            if(not site.HasWall(i))
            {
              continue;
            }

            // Direction of streaming from wall to this site
            LatticePosition const direction = LatticePosition(
                                                LATTICE::CX[i], LATTICE::CY[i], LATTICE::CZ[i]
                                              );
            LatticePosition const wallnode = LatticePosition(siteIn)
                                             + direction.GetNormalised() * distance;
            site.AddToForce(cell.nodeWall(vertex, wallnode) * weight);
          }
        }
      protected:
        Cell const &cell;
    };

  }
} // namespace details::anonymous

