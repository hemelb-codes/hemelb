//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/Interpolation.h"
#include <cmath>

namespace hemelb { namespace redblood {

  void IndexIterator :: operator++() {
#     ifndef NDEBUG
        if(not isValid())
          throw Exception() << "Cannot increment invalid iterator\n";
#     endif
      if(current_[2] < max_[2]) {
        ++current_[2];
        return;
      }
      current_[2] = min_[2];
      if(current_[1] < max_[1]) {
        ++current_[1];
        return;
      }
      current_[1] = min_[1];
      ++current_[0];
  }

  namespace {
    int minimumPosImpl(Dimensionless _x, size_t _range) {
        return static_cast<int>(
          std::floor(_x - 0.5 * Dimensionless(_range)) + 1
        );
    }
    int maximumPosImpl(Dimensionless _x, size_t _range) {
        return static_cast<int>(
          std::floor(_x + 0.5 * Dimensionless(_range))
        );
    }
  }

  LatticeVector InterpolationIterator :: minimumPosition_(
      LatticePosition const &_node,
      size_t _range) {
    return LatticeVector(
        minimumPosImpl(_node.x, _range),
        minimumPosImpl(_node.y, _range),
        minimumPosImpl(_node.z, _range)
    );
  }
  LatticeVector InterpolationIterator :: maximumPosition_(
      LatticePosition const &_node,
      size_t _range) {
    return LatticeVector(
        maximumPosImpl(_node.x, _range),
        maximumPosImpl(_node.y, _range),
        maximumPosImpl(_node.z, _range)
    );
  }

  InterpolationIterator interpolationIterator(
      LatticePosition const &_where, stencil::types _stencil) {
#   define HEMELB_STENCIL_MACRO(NAME, Name)                     \
        case stencil::NAME:                                     \
          return interpolationIterator<stencil::Name>(_where)
      switch(_stencil) {
        HEMELB_STENCIL_MACRO(FOUR_POINT, FourPoint);
        HEMELB_STENCIL_MACRO(COSINE_APPROX, CosineApprox);
        HEMELB_STENCIL_MACRO(THREE_POINT, ThreePoint);
        HEMELB_STENCIL_MACRO(TWO_POINT, TwoPoint);
      }
#   undef HEMELB_STENCIL_MACRO
  }


}}
