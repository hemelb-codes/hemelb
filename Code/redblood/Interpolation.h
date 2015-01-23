//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_INTERPOLATION_H
#define HEMELB_REDBLOOD_INTERPOLATION_H

#include "units.h"
#include "Exception.h"
#include "geometry/LatticeData.h"
#include "redblood/stencil.h"

#include <cassert>
#include <boost/shared_array.hpp>

namespace hemelb
{
  namespace redblood
  {

    //! Iterates over points on the lattice
    class IndexIterator;
    //! \brief Interpolates onto an off-lattice point
    //! \details Given an off-lattice point and a stencil, iterates over the
    //! points on the lattice which have non-zero weight.
    class InterpolationIterator;

    //! Creates an interpolator for a given stencil
    template<class STENCIL>
    InterpolationIterator interpolationIterator(LatticePosition const &in);
    //! Creates an interpolator for a given stencil
    InterpolationIterator interpolationIterator(
      LatticePosition const &where, stencil::types stencil);

    class IndexIterator
    {

      public:
        IndexIterator(LatticeVector const &center, LatticeCoordinate width)
          : min(center - LatticeVector::Ones() * width),
            max(center + LatticeVector::Ones() * width), current(min) {}
        IndexIterator(LatticeVector const &min, LatticeVector const &max)
          : min(min), max(max), current(min) {}

        //! Returns current position
        LatticeVector const & operator*() const
        {
          return current;
        }

        //! Increments to next position on lattice
        void operator++();

        //! True if iterator is still valid
        bool isValid() const
        {
          return current[0] <= max[0];
        }
        //! True if iterator is still valid
        operator bool() const
        {
          return isValid();
        }

        //! computes local index
        site_t ContiguousSiteId(geometry::LatticeData const &latDat) const
        {
          return latDat.GetContiguousSiteId(current);
        }
        // //! computes local index
        // geometry::Site<geometry::LatticeData const>
        //   Site(geometry::LatticeData const &latDat) const {
        //     return latDat.GetSite(current);
        // }
        // //! computes local index
        // geometry::Site<geometry::LatticeData>
        //   Site(geometry::LatticeData const &latDat) const {
        //     return latDat.GetSite(current);
        // }

      protected:
        //! Minimum indices
        LatticeVector const min;
        //! Maximum indices
        LatticeVector const max;
        //! Current point
        LatticeVector current;
    };

    class InterpolationIterator : public IndexIterator
    {
      public:
        //! Lattice
        template<class STENCIL>
        InterpolationIterator(
          LatticePosition const &node, STENCIL const & stencil);

        //! Returns weight for current point
        Dimensionless weight() const
        {
          assert(xWeight && yWeight && zWeight);
          assert(current[0] >= min[0]);
          assert(current[0] <= max[0]);
          assert(current[1] >= min[1]);
          assert(current[1] <= max[1]);
          assert(current[2] >= min[2]);
          assert(current[2] <= max[2]);
          return xWeight[current[0] - min[0]]
                 * yWeight[current[1] - min[1]]
                 * zWeight[current[2] - min[2]];
        }
        //! Weights for each direction
        util::Vector3D<Dimensionless> weights() const
        {
          assert(xWeight && yWeight && zWeight);
          assert(current[0] >= min[0]);
          assert(current[0] <= max[0]);
          assert(current[1] >= min[1]);
          assert(current[1] <= max[1]);
          assert(current[2] >= min[2]);
          assert(current[2] <= max[2]);
          return util::Vector3D<Dimensionless>(
                   xWeight[current[0] - min[0]],
                   yWeight[current[1] - min[1]],
                   zWeight[current[2] - min[2]]
                 );
        }


      protected:
        //! Weight alongst x direction;
        boost::shared_array<Dimensionless> xWeight;
        //! Weight alongst y direction;
        boost::shared_array<Dimensionless> yWeight;
        //! Weight alongst z direction;
        boost::shared_array<Dimensionless> zWeight;

        static LatticeVector minimumPosition(LatticePosition const &node,
                                             size_t range);
        static LatticeVector maximumPosition(LatticePosition const &node,
                                             size_t range);
    };

    template<class STENCIL>
    InterpolationIterator :: InterpolationIterator(
      LatticePosition const &node, STENCIL const & stencil)
      : IndexIterator(
          minimumPosition(node, STENCIL::range),
          maximumPosition(node, STENCIL::range)
        ), xWeight(new Dimensionless[STENCIL::range]),
        yWeight(new Dimensionless[STENCIL::range]),
        zWeight(new Dimensionless[STENCIL::range])
    {
      for(LatticeVector::value_type i(0); i < STENCIL::range; ++i)
      {
        xWeight[i] = stencil(node[0] - Dimensionless(min[0] + i));
        yWeight[i] = stencil(node[1] - Dimensionless(min[1] + i));
        zWeight[i] = stencil(node[2] - Dimensionless(min[2] + i));
      }
    }

    template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
                                 InterpolationIterator interpolator)
    {
      PhysicalVelocity result(0, 0, 0);

      for(; interpolator; ++interpolator)
      {
        result += gridfunc(*interpolator) * interpolator.weight();
      }

      return result;
    }
    template<class GRID_FUNCTION, class STENCIL>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
                                 LatticePosition const &pos, STENCIL stencil)
    {
      return interpolate(gridfunc, OffLatticeIterator(pos, stencil));
    }
    template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
                                 LatticePosition const &pos, stencil::types stencil)
    {
      return interpolate(gridfunc, interpolationIterator(pos, stencil));
    }
    template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
                                 Dimensionless const &x, Dimensionless const &y,
                                 Dimensionless const &z, stencil::types stencil)
    {
      return interpolate(gridfunc, LatticePosition(x, y, z), stencil);
    }

    // Creates an interpolator for a given stencil
    template<class STENCIL>
    InterpolationIterator interpolationIterator(LatticePosition const &in)
    {
      return InterpolationIterator(in, STENCIL());
    }
  }
} // hemelb::redblood
#endif
