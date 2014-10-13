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
  class RegionIterator {

    public:
      RegionIterator(LatticeVector const &_center, LatticeCoordinate _width)
          : min_(_center - LatticeVector::Ones() * _width),
            max_(_center + LatticeVector::Ones() * _width), current_(min_) {}
      RegionIterator(LatticeVector const &_min, LatticeVector const &_max)
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

}} // hemelb::redblood
#endif
