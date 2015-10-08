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
#include <array>

namespace hemelb
{
  namespace redblood
  {

    namespace
    {
      inline int minimumPosImpl(Dimensionless x, size_t range)
      {
        return static_cast<int>(std::floor(x - 0.5 * Dimensionless(range)) + 1);
      }
      inline int maximumPosImpl(Dimensionless x, size_t range)
      {
        return static_cast<int>(std::floor(x + 0.5 * Dimensionless(range)));
      }
    }

    //! Iterates over points on the lattice
    class IndexIterator;
    //! \brief Interpolates onto an off-lattice point
    //! \details Given an off-lattice point and a stencil, iterates over the
    //! points on the lattice which have non-zero weight.
    template <class Stencil>
    class InterpolationIterator;

    //! Creates an interpolator for a given stencil
    template <class Stencil>
    InterpolationIterator<Stencil> interpolationIterator(LatticePosition const &in);

    //! \brief Iterates over lattice points surrounding a given coordinate
    //! \details Given a coordinate :math:`\Omega`, this object iterates over integer coordinates
    //! that are less than a given distance from it. The distance is an L1 norm (iteration is over a
    //! cube)
    class IndexIterator
    {
      public:
        IndexIterator(LatticeVector const &center, LatticeCoordinate width) :
            min(center - LatticeVector::Ones() * width),
                max(center + LatticeVector::Ones() * width), current(min)
        {
        }
        IndexIterator(LatticeVector const &min, LatticeVector const &max) :
            min(min), max(max), current(min)
        {
        }

        //! Returns current position
        LatticeVector const &operator*() const
        {
          return current;
        }

        //! Increments to next position on lattice
        void operator++();

        //! True if iterator is still valid
        bool IsValid() const
        {
          return current[0] <= max[0];
        }
        //! True if iterator is still valid
        operator bool() const
        {
          return IsValid();
        }

        //! computes local index
        site_t ContiguousSiteId(geometry::LatticeData const &latDat) const
        {
          return latDat.GetContiguousSiteId(current);
        }

      protected:
        //! Minimum indices
        LatticeVector min;
        //! Maximum indices
        LatticeVector max;
        //! Current point
        LatticeVector current;
    };

    //! \brief Iterates over lattice points close to a given coordinate
    //! \details Adds knowledge of weights for iteration
    template <class Stencil>
    class InterpolationIterator : public IndexIterator
    {
      public:
        //! Lattice
        InterpolationIterator(LatticePosition const &node);

        //! Returns weight for current point
        Dimensionless weight() const
        {
          assert(current[0] >= min[0]);
          assert(current[0] <= max[0]);
          assert(current[1] >= min[1]);
          assert(current[1] <= max[1]);
          assert(current[2] >= min[2]);
          assert(current[2] <= max[2]);
          return xWeight[current[0] - min[0]] * yWeight[current[1] - min[1]]
              * zWeight[current[2] - min[2]];
        }
        //! Weights for each direction
        util::Vector3D<Dimensionless> weights() const
        {
          assert(current[0] >= min[0]);
          assert(current[0] <= max[0]);
          assert(current[1] >= min[1]);
          assert(current[1] <= max[1]);
          assert(current[2] >= min[2]);
          assert(current[2] <= max[2]);
          return util::Vector3D<Dimensionless>(xWeight[current[0] - min[0]],
                                               yWeight[current[1] - min[1]],
                                               zWeight[current[2] - min[2]]);
        }

      protected:
        //! Weight alongst x direction;
        std::array<Dimensionless, Stencil::GetRange()> xWeight;
        //! Weight alongst y direction;
        std::array<Dimensionless, Stencil::GetRange()> yWeight;
        //! Weight alongst z direction;
        std::array<Dimensionless, Stencil::GetRange()> zWeight;

        static LatticeVector minimumPosition(LatticePosition const &node, size_t range);
        static LatticeVector maximumPosition(LatticePosition const &node, size_t range);
    };

    template <class Stencil>
    LatticeVector InterpolationIterator<Stencil>::minimumPosition(LatticePosition const &node, size_t range)
    {
      return LatticeVector(minimumPosImpl(node.x, range),
                           minimumPosImpl(node.y, range),
                           minimumPosImpl(node.z, range));
    }
    template <class Stencil>
    LatticeVector InterpolationIterator<Stencil>::maximumPosition(LatticePosition const &node, size_t range)
    {
      return LatticeVector(maximumPosImpl(node.x, range),
                           maximumPosImpl(node.y, range),
                           maximumPosImpl(node.z, range));
    }

    template <class Stencil>
    InterpolationIterator<Stencil>::InterpolationIterator(LatticePosition const &node) :
        IndexIterator(
            minimumPosition(node, Stencil::GetRange()),
            maximumPosition(node, Stencil::GetRange())
        )
    {
      for (LatticeVector::value_type i(0); i < LatticeVector::value_type(Stencil::GetRange()); ++i)
      {
        xWeight[i] = Stencil::stencil(node[0] - Dimensionless(min[0] + i));
        yWeight[i] = Stencil::stencil(node[1] - Dimensionless(min[1] + i));
        zWeight[i] = Stencil::stencil(node[2] - Dimensionless(min[2] + i));
      }
    }

    template<class GRID_FUNCTION, class Stencil>
    LatticeVelocity interpolate(GRID_FUNCTION const &gridfunc, InterpolationIterator<Stencil> interpolator)
    {
      LatticeVelocity result(0, 0, 0);

      for (; interpolator; ++interpolator)
      {
        result += gridfunc(*interpolator) * interpolator.weight();
      }

      return result;
    }
    template<class GRID_FUNCTION, class Stencil>
    LatticeVelocity interpolate(GRID_FUNCTION const &gridfunc, LatticePosition const &pos)
    {
      return interpolate(gridfunc, interpolationIterator<Stencil>(pos));
    }
    template<class GRID_FUNCTION, class Stencil>
    LatticeVelocity interpolate(GRID_FUNCTION const &gridfunc, Dimensionless const &x,
                                 Dimensionless const &y, Dimensionless const &z)
    {
      return interpolate<GRID_FUNCTION, Stencil>(gridfunc, LatticePosition(x, y, z));
    }

    // Creates an interpolator for a given stencil
    template <class Stencil>
    InterpolationIterator<Stencil> interpolationIterator(LatticePosition const &in)
    {
      return InterpolationIterator<Stencil>(in);
    }
  }
} // hemelb::redblood
#endif
