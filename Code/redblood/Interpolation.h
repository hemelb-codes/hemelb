// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
      inline int minimumPosImpl(Dimensionless const x, size_t const range)
      {
        return static_cast<int>(std::floor(x - 0.5 * Dimensionless(range)) + 1);
      }
      inline int maximumPosImpl(Dimensionless const x, size_t const range)
      {
        return minimumPosImpl(x, range) + range - 1;
      }
    }

    //! Iterates over points on the lattice
    class IndexIterator;
    //! \brief Interpolates onto an off-lattice point
    //! \details Given an off-lattice point and a stencil, iterates over the
    //! points on the lattice which have non-zero weight.
    template<class STENCIL>
    class InterpolationIterator;

    //! Creates an interpolator for a given stencil
    template<class STENCIL>
    InterpolationIterator<STENCIL> interpolationIterator(LatticePosition const &in);

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
        //! Returns current position
        LatticeVector const *operator->() const
        {
          return &current;
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
    template<class STENCIL>
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
          assert(current[0] - min[0] < LatticeCoordinate(xWeight.size()));
          assert(current[1] - min[1] < LatticeCoordinate(yWeight.size()));
          assert(current[2] - min[2] < LatticeCoordinate(zWeight.size()));
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
          assert(current[0] - min[0] < xWeight.size());
          assert(current[1] - min[1] < yWeight.size());
          assert(current[2] - min[2] < zWeight.size());
          return util::Vector3D<Dimensionless>(xWeight[current[0] - min[0]],
                                               yWeight[current[1] - min[1]],
                                               zWeight[current[2] - min[2]]);
        }

      protected:
        //! Weight alongst x direction;
        std::array<Dimensionless, STENCIL::GetRange()> xWeight;
        //! Weight alongst y direction;
        std::array<Dimensionless, STENCIL::GetRange()> yWeight;
        //! Weight alongst z direction;
        std::array<Dimensionless, STENCIL::GetRange()> zWeight;

        static LatticeVector minimumPosition(LatticePosition const &node, size_t range);
        static LatticeVector maximumPosition(LatticePosition const &node, size_t range);
    };

    template<class STENCIL>
    LatticeVector InterpolationIterator<STENCIL>::minimumPosition(LatticePosition const &node,
                                                                  size_t range)
    {
      return LatticeVector(minimumPosImpl(node.x, range),
                           minimumPosImpl(node.y, range),
                           minimumPosImpl(node.z, range));
    }
    template<class STENCIL>
    LatticeVector InterpolationIterator<STENCIL>::maximumPosition(LatticePosition const &node,
                                                                  size_t range)
    {
      return LatticeVector(maximumPosImpl(node.x, range),
                           maximumPosImpl(node.y, range),
                           maximumPosImpl(node.z, range));
    }

    template<class STENCIL>
    InterpolationIterator<STENCIL>::InterpolationIterator(LatticePosition const &node) :
            IndexIterator(minimumPosition(node, STENCIL::GetRange()),
                          maximumPosition(node, STENCIL::GetRange()))
    {
      // Lattice interpolation range contains the specified off-lattice position
      assert(node[0] >= min[0]);
      assert(node[0] < max[0]);
      assert(node[1] >= min[1]);
      assert(node[1] < max[1]);
      assert(node[2] >= min[2]);
      assert(node[2] < max[2]);

      // The size of the interpolation range is consistent with the stencil size
      assert(max[0] - min[0] + 1 == STENCIL::GetRange());
      assert(max[1] - min[1] + 1 == STENCIL::GetRange());
      assert(max[2] - min[2] + 1 == STENCIL::GetRange());

      for (LatticeVector::value_type i(0); i < LatticeVector::value_type(STENCIL::GetRange()); ++i)
      {
        xWeight[i] = STENCIL::stencil(node[0] - Dimensionless(min[0] + i));
        yWeight[i] = STENCIL::stencil(node[1] - Dimensionless(min[1] + i));
        zWeight[i] = STENCIL::stencil(node[2] - Dimensionless(min[2] + i));
      }
    }

    template<class GRID_FUNCTION, class STENCIL>
    LatticeVelocity interpolate(GRID_FUNCTION const &gridfunc,
                                InterpolationIterator<STENCIL> interpolator)
    {
      LatticeVelocity result(0, 0, 0);

      for (; interpolator; ++interpolator)
      {
        result += gridfunc(*interpolator) * interpolator.weight();
      }

      return result;
    }
    template<class GRID_FUNCTION, class STENCIL>
    LatticeVelocity interpolate(GRID_FUNCTION const &gridfunc, LatticePosition const &pos)
    {
      return interpolate(gridfunc, interpolationIterator<STENCIL>(pos));
    }
    template<class GRID_FUNCTION, class STENCIL>
    LatticeVelocity interpolate(GRID_FUNCTION const &gridfunc, Dimensionless const &x,
                                Dimensionless const &y, Dimensionless const &z)
    {
      return interpolate<GRID_FUNCTION, STENCIL>(gridfunc, LatticePosition(x, y, z));
    }

    // Creates an interpolator for a given stencil
    template<class STENCIL>
    InterpolationIterator<STENCIL> interpolationIterator(LatticePosition const &in)
    {
      return InterpolationIterator<STENCIL>(in);
    }
  }
} // hemelb::redblood
#endif
