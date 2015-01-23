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
    int minimumPosImpl(Dimensionless x, size_t range) {
        return static_cast<int>(
          std::floor(x - 0.5 * Dimensionless(range)) + 1
        );
    }
    int maximumPosImpl(Dimensionless x, size_t range) {
        return static_cast<int>(
          std::floor(x + 0.5 * Dimensionless(range))
        );
    }
  }

  LatticeVector InterpolationIterator :: minimumPosition_(
      LatticePosition const &node,
      size_t range) {
    return LatticeVector(
        minimumPosImpl(node.x, range),
        minimumPosImpl(node.y, range),
        minimumPosImpl(node.z, range)
    );
  }
  LatticeVector InterpolationIterator :: maximumPosition_(
      LatticePosition const &node,
      size_t range) {
    return LatticeVector(
        maximumPosImpl(node.x, range),
        maximumPosImpl(node.y, range),
        maximumPosImpl(node.z, range)
    );
  }

  InterpolationIterator interpolationIterator(
      LatticePosition const &where, stencil::types stencil) {
#   define HEMELB_STENCIL_MACRO(NAME, Name)                     \
        case stencil::NAME:                                     \
          return interpolationIterator<stencil::Name>(where)
      switch(stencil) {
        HEMELB_STENCIL_MACRO(FOUR_POINT, FourPoint);
        HEMELB_STENCIL_MACRO(COSINE_APPROX, CosineApprox);
        HEMELB_STENCIL_MACRO(THREE_POINT, ThreePoint);
        HEMELB_STENCIL_MACRO(TWO_POINT, TwoPoint);
      }
#   undef HEMELB_STENCIL_MACRO
  }


}}
