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
    InterpolationIterator interpolationIterator(LatticePosition const &_in);
  //! Creates an interpolator for a given stencil
  InterpolationIterator interpolationIterator(
      LatticePosition const &_where, stencil::types _stencil);

  class IndexIterator {

    public:
      IndexIterator(LatticeVector const &_center, LatticeCoordinate _width)
          : min_(_center - LatticeVector::Ones() * _width),
            max_(_center + LatticeVector::Ones() * _width), current_(min_) {}
      IndexIterator(LatticeVector const &_min, LatticeVector const &_max)
          : min_(_min), max_(_max), current_(_min) {}

      //! Returns current position
      LatticeVector const & operator*() const { return current_; }

      //! Increments to next position on lattice
      void operator++();

      //! True if iterator is still valid
      bool isValid() const { return current_[0] <= max_[0]; }
      //! True if iterator is still valid
      operator bool() const { return isValid(); }

      //! computes local index
      site_t ContiguousSiteId(geometry::LatticeData const &_latDat) const{
        return _latDat.GetContiguousSiteId(current_);
      }
      // //! computes local index
      // geometry::Site<geometry::LatticeData const>
      //   Site(geometry::LatticeData const &_latDat) const {
      //     return _latDat.GetSite(current_);
      // }
      // //! computes local index
      // geometry::Site<geometry::LatticeData>
      //   Site(geometry::LatticeData const &_latDat) const {
      //     return _latDat.GetSite(current_);
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
            LatticePosition const &_node, STENCIL const & _stencil);

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

      static LatticeVector minimumPosition_(LatticePosition const &_node,
              size_t _range);
      static LatticeVector maximumPosition_(LatticePosition const &_node,
              size_t _range);
  };

  template<class STENCIL>
    InterpolationIterator :: InterpolationIterator(
         LatticePosition const &_node, STENCIL const & _stencil)
      : IndexIterator(
          minimumPosition_(_node, STENCIL::range),
          maximumPosition_(_node, STENCIL::range)
        ), xWeight_(new Dimensionless[STENCIL::range]),
           yWeight_(new Dimensionless[STENCIL::range]),
           zWeight_(new Dimensionless[STENCIL::range]) {
      for(LatticeVector::value_type i(0); i < STENCIL::range; ++i) {
        xWeight_[i] = _stencil(_node[0] - Dimensionless(min_[0] + i));
        yWeight_[i] = _stencil(_node[1] - Dimensionless(min_[1] + i));
        zWeight_[i] = _stencil(_node[2] - Dimensionless(min_[2] + i));
      }
    }

  template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &_gridfunc,
        InterpolationIterator _interpolator) {
      PhysicalVelocity result(0, 0, 0);
      for(; _interpolator; ++_interpolator)
        result += _gridfunc(*_interpolator) * _interpolator.weight();
      return result;
    }
  template<class GRID_FUNCTION, class STENCIL>
    PhysicalVelocity interpolate(GRID_FUNCTION const &_gridfunc,
        LatticePosition const &_pos, STENCIL _stencil) {
      return interpolate(_gridfunc, OffLatticeIterator(_pos, _stencil));
    }
  template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &_gridfunc,
        LatticePosition const &_pos, stencil::types _stencil) {
      return interpolate(_gridfunc, interpolationIterator(_pos, _stencil));
    }
  template<class GRID_FUNCTION>
    PhysicalVelocity interpolate(GRID_FUNCTION const &_gridfunc,
        Dimensionless const &_x, Dimensionless const &_y,
        Dimensionless const &_z, stencil::types _stencil) {
      return interpolate(_gridfunc, LatticePosition(_x, _y, _z), _stencil);
  }

  // Creates an interpolator for a given stencil
  template<class STENCIL>
    InterpolationIterator interpolationIterator(LatticePosition const &_in) {
      return InterpolationIterator(_in, STENCIL());
    }
}} // hemelb::redblood
#endif
