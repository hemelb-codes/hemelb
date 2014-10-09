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
    // Four point stencil
    inline Dimensionless fourPoint(Dimensionless const _x) {
      Dimensionless x(std::abs(_x));
      if(x < 1)
        return 1./8. * (3. - 2*x + std::sqrt(1. + 4.*x - 4.*x*x));
      else if(x < 2)
        return 1./8. * (5. - 2*x - std::sqrt(-7. + 12.*x - 4.*x*x));
      else
        return 0.;
    }

    // Approximation to the four-point stencil
    inline Dimensionless cosineApprox(Dimensionless const _x) {
      return std::abs(_x) < 2 ? 0.25 * (1. + std::cos(PI * _x * 0.5)): 0.;
    }

    // Three-point stencil
    inline Dimensionless threePoint(Dimensionless const _x) {
      Dimensionless x(std::abs(_x));
      if(x < 0.5)
        return 1./3. * (1 + std::sqrt(1. - 3.*x*x));
      else if(x < 1.5)
        return 1./6. * (5. - 3.*x - std::sqrt(-2. + 6.*x - 3.*x*x));
      else
        return 0.;
    }

    // Two-point stencil
    inline Dimensionless twoPoint(Dimensionless const _x) {
      Dimensionless x(std::abs(_x));
      return x < 1 ? 1. - x: 0;
    }

  }

}}

#endif
