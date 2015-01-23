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

namespace hemelb { namespace redblood {

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

  class IndexIterator {

    public:
      IndexIterator(LatticeVector const &center, LatticeCoordinate width)
          : min_(center - LatticeVector::Ones() * width),
            max_(center + LatticeVector::Ones() * width), current_(min_) {}
      IndexIterator(LatticeVector const &min, LatticeVector const &max)
          : min_(min), max_(max), current_(min) {}

      //! Returns current position
      LatticeVector const & operator*() const { return current_; }

      //! Increments to next position on lattice
      void operator++();

      //! True if iterator is still valid
      bool isValid() const { return current_[0] <= max_[0]; }
      //! True if iterator is still valid
      operator bool() const { return isValid(); }

      //! computes local index
      site_t ContiguousSiteId(geometry::LatticeData const &latDat) const{
        return latDat.GetContiguousSiteId(current_);
      }
      // //! computes local index
      // geometry::Site<geometry::LatticeData const>
      //   Site(geometry::LatticeData const &latDat) const {
      //     return latDat.GetSite(current_);
      // }
      // //! computes local index
      // geometry::Site<geometry::LatticeData>
      //   Site(geometry::LatticeData const &latDat) const {
      //     return latDat.GetSite(current_);
      // }

    protected:
      //! Minimum indices
      LatticeVector const min_;
      //! Maximum indices
      LatticeVector const max_;
      //! Current point
      LatticeVector current_;
  };

  class InterpolationIterator : public IndexIterator {
    public:
      //! Lattice
      template<class STENCIL>
        InterpolationIterator(
            LatticePosition const &node, STENCIL const & stencil);

      //! Returns weight for current point
      Dimensionless weight() const {
        assert(xWeight_ && yWeight_ && zWeight_);
        assert(current_[0] >= min_[0]);
        assert(current_[0] <= max_[0]);
        assert(current_[1] >= min_[1]);
        assert(current_[1] <= max_[1]);
        assert(current_[2] >= min_[2]);
        assert(current_[2] <= max_[2]);
        return xWeight_[current_[0] - min_[0]]
          * yWeight_[current_[1] - min_[1]]
          * zWeight_[current_[2] - min_[2]];
      }
      //! Weights for each direction
      util::Vector3D<Dimensionless> weights() const {
        assert(xWeight_ && yWeight_ && zWeight_);
        assert(current_[0] >= min_[0]);
        assert(current_[0] <= max_[0]);
        assert(current_[1] >= min_[1]);
        assert(current_[1] <= max_[1]);
        assert(current_[2] >= min_[2]);
        assert(current_[2] <= max_[2]);
        return util::Vector3D<Dimensionless>(
            xWeight_[current_[0] - min_[0]],
            yWeight_[current_[1] - min_[1]],
            zWeight_[current_[2] - min_[2]]
        );
      }


    protected:
      //! Weight alongst x direction;
      boost::shared_array<Dimensionless> xWeight_;
      //! Weight alongst y direction;
      boost::shared_array<Dimensionless> yWeight_;
      //! Weight alongst z direction;
      boost::shared_array<Dimensionless> zWeight_;

      static LatticeVector minimumPosition_(LatticePosition const &node,
              size_t range);
      static LatticeVector maximumPosition_(LatticePosition const &node,
              size_t range);
  };

  template<class STENCIL>
    InterpolationIterator :: InterpolationIterator(
         LatticePosition const &node, STENCIL const & stencil)
      : IndexIterator(
          minimumPosition_(node, STENCIL::range),
          maximumPosition_(node, STENCIL::range)
        ), xWeight_(new Dimensionless[STENCIL::range]),
           yWeight_(new Dimensionless[STENCIL::range]),
           zWeight_(new Dimensionless[STENCIL::range]) {
      for(LatticeVector::value_type i(0); i < STENCIL::range; ++i) {
        xWeight_[i] = stencil(node[0] - Dimensionless(min_[0] + i));
        yWeight_[i] = stencil(node[1] - Dimensionless(min_[1] + i));
        zWeight_[i] = stencil(node[2] - Dimensionless(min_[2] + i));
      }
    }

  template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
        InterpolationIterator interpolator) {
      PhysicalVelocity result(0, 0, 0);
      for(; interpolator; ++interpolator)
        result += gridfunc(*interpolator) * interpolator.weight();
      return result;
    }
  template<class GRID_FUNCTION, class STENCIL>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
        LatticePosition const &pos, STENCIL stencil) {
      return interpolate(gridfunc, OffLatticeIterator(pos, stencil));
    }
  template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
        LatticePosition const &pos, stencil::types stencil) {
      return interpolate(gridfunc, interpolationIterator(pos, stencil));
    }
  template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &gridfunc,
        Dimensionless const &x, Dimensionless const &y,
        Dimensionless const &z, stencil::types stencil) {
      return interpolate(gridfunc, LatticePosition(x, y, z), stencil);
  }

  // Creates an interpolator for a given stencil
  template<class STENCIL>
    InterpolationIterator interpolationIterator(LatticePosition const &in) {
      return InterpolationIterator(in, STENCIL());
    }
}} // hemelb::redblood
#endif
