// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_STENCIL_H
#define HEMELB_REDBLOOD_STENCIL_H

#include <cmath>
#include <type_traits>
#include "units.h"
#include "constants.h"

namespace hemelb
{
  namespace redblood
  {
    namespace stencil
    {
      //! Constant to name stencils without referring to type
      //! Useful to create factories
      enum class types
        : unsigned int
        {
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
      inline Dimensionless fourPoint(Dimensionless const x)
      {
        Dimensionless xAbs(std::abs(x));

        if (xAbs < 1)
        {
          return 1. / 8. * (3. - 2 * xAbs + std::sqrt(1. + 4. * xAbs - 4. * xAbs * xAbs));
        }
        else if (xAbs < 2)
        {
          return 1. / 8. * (5. - 2 * xAbs - std::sqrt(-7. + 12. * xAbs - 4. * xAbs * xAbs));
        }
        else
        {
          return 0.;
        }
      }

      // Approximation to the four-point stencil
      inline Dimensionless cosineApprox(Dimensionless const x)
      {
        return std::abs(x) < 2 ?
          0.25 * (1. + std::cos(PI * x * 0.5)) :
          0.;
      }

      // Three-point stencil
      inline Dimensionless threePoint(Dimensionless const x)
      {
        Dimensionless xAbs(std::abs(x));

        if (xAbs < 0.5)
        {
          return 1. / 3. * (1 + std::sqrt(1. - 3. * xAbs * xAbs));
        }
        else if (xAbs < 1.5)
        {
          return 1. / 6. * (5. - 3. * xAbs - std::sqrt(-2. + 6. * xAbs - 3. * xAbs * xAbs));
        }
        else
        {
          return 0.;
        }
      }

      // Two-point stencil
      inline Dimensionless twoPoint(Dimensionless const x)
      {
        Dimensionless xAbs(std::abs(x));
        return xAbs < 1 ?
          1. - xAbs :
          0;
      }

#define HEMELB_STENCIL_MACRO(NAME, STENCIL, RANGE)           \
  struct NAME                                                \
  {                                                          \
    static Dimensionless stencil(Dimensionless x)            \
    {                                                        \
      return STENCIL(x);                                     \
    }                                                        \
    static Dimensionless stencil(LatticePosition const &x)   \
    {                                                        \
      return STENCIL(x.x) * STENCIL(x.y) * STENCIL(x.z);     \
    }                                                        \
    static constexpr size_t GetRange()                       \
    {                                                        \
      return RANGE;                                          \
    }                                                        \
  };                                                         \
  static_assert(std::is_pod<NAME>::value, "Can be a struct")
      HEMELB_STENCIL_MACRO(FourPoint, fourPoint, 4);
      HEMELB_STENCIL_MACRO(CosineApprox, cosineApprox, 4);
      HEMELB_STENCIL_MACRO(ThreePoint, threePoint, 3);
      HEMELB_STENCIL_MACRO(TwoPoint, twoPoint, 2);
#undef HEMELB_STENCIL_MACRO
    }
  }
}

#endif
