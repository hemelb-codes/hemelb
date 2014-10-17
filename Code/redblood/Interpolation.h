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

namespace hemelb { namespace redblood {

  //! Iterates over points on the lattice
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

    protected:
      //! Minimum indices
      LatticeVector const min_;
      //! Maximum indices
      LatticeVector const max_;
      //! Current point
      LatticeVector current_;
  };

  //! \brief Interpolates onto an off-lattice point
  //! \details Given an off-lattice point and a stencil, iterates over the
  //! points on the lattice which have non-zero weight.
  class OffLatticeInterpolator : public IndexIterator {
    public:
      //! Lattice
      template<class STENCIL>
        OffLatticeInterpolator(
            LatticePosition const &_node, STENCIL const & _stencil);

      //! Returns weight for current point
      Dimensionless weight() const {
        return xWeight_[current_[0] - min_[0]]
          * yWeight_[current_[1] - min_[1]]
          * zWeight_[current_[2] - min_[2]];
      }


    protected:
      //! Weight alongst x direction;
      std::vector<Dimensionless> xWeight_;
      //! Weight alongst y direction;
      std::vector<Dimensionless> yWeight_;
      //! Weight alongst z direction;
      std::vector<Dimensionless> zWeight_;

      static LatticeVector minimumPosition_(LatticePosition const &_node,
              size_t _range);
      static LatticeVector maximumPosition_(LatticePosition const &_node,
              size_t _range);
  };

  template<class STENCIL>
    OffLatticeInterpolator :: OffLatticeInterpolator(
         LatticePosition const &_node, STENCIL const & _stencil)
      : IndexIterator(
          minimumPosition_(_node, STENCIL::range),
          maximumPosition_(_node, STENCIL::range)
        ), xWeight_(STENCIL::range), yWeight_(STENCIL::range),
        zWeight_(STENCIL::range) {
      for(size_t i(0); i < STENCIL::range; ++i) {
        xWeight_[i] = _stencil(_node[0] - Dimensionless(min_[0] + i));
        yWeight_[i] = _stencil(_node[1] - Dimensionless(min_[1] + i));
        zWeight_[i] = _stencil(_node[2] - Dimensionless(min_[2] + i));
      }
    }

}} // hemelb::redblood
#endif
