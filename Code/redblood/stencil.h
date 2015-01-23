//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_STENCIL_H
#define HEMELB_REDBLOOD_STENCIL_H

#include <boost/core/enable_if.hpp>
#include <cmath>
#include "units.h"
#include "constants.h"

namespace hemelb { namespace redblood {

  namespace stencil {
    //! Constant to name stencils without referring to type
    //! Useful to create factories
    enum types {
      FOUR_POINT,
      COSINE_APPROX,
      THREE_POINT,
      TWO_POINT
    };
    //! Four point stencil
    struct FourPoint;
    //! Approximation to the four-point stencil
    struct CosineApprox;
    //! Three-point stencil
    struct ThreePoint;
    //! Four point stencil
    struct TwoPoint;

    // Four point stencil
    inline Dimensionless fourPoint(Dimensionless const _x) {
      Dimensionless xAbs(std::abs(_x));
      if(xAbs < 1)
        return 1./8. * (3. - 2*xAbs + std::sqrt(1. + 4.*xAbs - 4.*xAbs*xAbs));
      else if(xAbs < 2)
        return 1./8. * (5. - 2*xAbs - std::sqrt(-7. + 12.*xAbs - 4.*xAbs*xAbs));
      else
        return 0.;
    }

    // Approximation to the four-point stencil
    inline Dimensionless cosineApprox(Dimensionless const _x) {
      return std::abs(_x) < 2 ? 0.25 * (1. + std::cos(PI * _x * 0.5)): 0.;
    }

    // Three-point stencil
    inline Dimensionless threePoint(Dimensionless const _x) {
      Dimensionless xAbs(std::abs(_x));
      if(xAbs < 0.5)
        return 1./3. * (1 + std::sqrt(1. - 3.*xAbs*xAbs));
      else if(xAbs < 1.5)
        return 1./6. * (5. - 3.*xAbs - std::sqrt(-2. + 6.*xAbs - 3.*xAbs*xAbs));
      else
        return 0.;
    }

    // Two-point stencil
    inline Dimensionless twoPoint(Dimensionless const _x) {
      Dimensionless xAbs(std::abs(_x));
      return xAbs < 1 ? 1. - xAbs: 0;
    }

#   define HEMELB_STENCIL_MACRO(NAME, STENCIL)                        \
      struct NAME {                                                   \
        NAME() {}                                                     \
        static Dimensionless stencil(Dimensionless _x) {              \
          return STENCIL(_x);                                         \
        }                                                             \
        static Dimensionless stencil(LatticePosition const &_x) {     \
          return STENCIL(_x.x) * STENCIL(_x.y) * STENCIL(_x.z);       \
        }                                                             \
        Dimensionless operator()(Dimensionless _x) const {            \
          return NAME::stencil(_x);                                   \
        }                                                             \
        Dimensionless operator()(LatticePosition const &_x) {         \
          return NAME::stencil(_x);                                   \
        }                                                             \
        static const size_t range;                                    \
      }
      HEMELB_STENCIL_MACRO(FourPoint, fourPoint);
      HEMELB_STENCIL_MACRO(CosineApprox, cosineApprox);
      HEMELB_STENCIL_MACRO(ThreePoint, threePoint);
      HEMELB_STENCIL_MACRO(TwoPoint, twoPoint);
#   undef HEMELB_STENCIL_MACRO
  }

}}

#endif
